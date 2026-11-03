// SPDX-License-Identifier: GPL-3.0-or-later
// meosar-rx: Cospas-Sarsat 406 MHz beacon bursts relayed on the MEOSAR L-band downlink
// (Galileo SAR 1544.05-1544.15 MHz, GLONASS-K 1544.85-1544.95 MHz).
//
// The input is channelised as it streams into one 125 kS/s stream per SAR band. Every 8 s
// block of a band stream (plus 0.3 s before and 0.9 s after) goes through the FGB chain:
// coherent CW detector -> coherent demodulator -> Chase-II BCH -> relay combining. One JSON
// record per candidate on stdout (and UDP), same fields as the Python prototype.
#include <getopt.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "detect.h"
#include "dsp.h"
#include "fgb.h"
#include "fgb_demod.h"
#include "input.h"
#include "out.h"
#include "sgb_rx.h"

#define D1 20
#define PRE 0.3
#define POST 1.3   // a whole SGB burst (1.0 s) must fit after the block end

typedef struct {
	const char *name;
	double lo, hi, mid;
	chan_stream_t *cs;
	cf_t *buf;            // band samples with absolute indices [base, base+len)
	long long base;
	size_t len, cap;
} band_t;

static volatile sig_atomic_t stop;
static void on_sig(int s) { (void)s; stop = 1; }

static void usage(void) {
	fprintf(stderr,
	        "usage: meosar-rx -i INPUT [options]\n"
	        "  -i INPUT         recording (.sigmf-meta/.sigmf-data or raw with -F), or airspy[:SERIAL],\n"
	        "                   hackrf[:SERIAL], rtl_tcp:HOST:PORT, spyserver:HOST:PORT\n"
	        "  -F FMT           raw file format: cu8, cs8, cs16, cf32\n"
	        "  -c HZ            centre frequency (live input; default 1544.5e6)\n"
	        "  -s RATE          sample rate (live input; default 2.5e6)\n"
	        "  -g GAIN -b       gain, bias-tee (antenna power) for live input\n"
	        "  --start S --dur S  part of a recording to process\n"
	        "  --bands LIST     galileo,glonass (default), quiet = noise-only test band\n"
	        "  --thr DB         detector threshold (13)    --L N  Chase order (10, 0 = hard)\n"
	        "  --no-combine     no relay combining         --no-llr  leave soft bits out of records\n"
	        "  --no-sgb         FGB only                   --sgb-threads N (4)  --sgb-span HZ (8000)\n"
	        "  --sgb-thr Z (19) acquisition threshold      --sgb-prez Z (7) preamble threshold\n"
	        "  -u HOST:PORT     UDP output (repeatable)    -q  nothing on stdout\n");
}

static const char *tier(int a) { return a == 2 ? "single" : a == 1 ? "weak" : ""; }

static void emit(const fgb_rec_t *r, const fgb_rec_t *all, double t0_ms, int with_llr) {
	jbuf_t b = {0};
	char hex[40] = "";
	if (r->frame.nbits) fgb_hex(r->frame.bits, r->frame.nbits, hex);
	jb_printf(&b, "{\"t\": %.3f, \"band\": \"%s\"", r->t_frame, r->band);
	if (r->combined) {
		jb_printf(&b, ", \"combined\": %d, \"relays\": [", r->combined);
		for (int k = 0; k < r->combined; k++) {
			const fgb_rec_t *q = &all[r->relay_idx[k]];
			char qh[40] = "";
			if (q->accept && q->frame.nbits) fgb_hex(q->frame.bits, q->frame.nbits, qh);
			jb_printf(&b, "%s{\"f_hz\": %.1f, \"ebn0\": %.2f, \"dt_ms\": %.2f, \"hex\": \"%s\"}", k ? ", " : "",
			          q->f_hz, q->ebn0, 1e3 * (q->t_cw - r->t_cw), qh);
		}
		jb_raw(&b, "]");
	} else {
		jb_printf(&b, ", \"det_snr\": %.1f", r->det_snr);
	}
	jb_printf(&b, ", \"t_cw\": %.6f, \"f_hz\": %.1f", r->t_cw, r->f_hz);
	if (!r->combined)
		jb_printf(&b, ", \"fdot\": %.2f, \"cn0\": %.2f, \"selftest_sync\": %s, \"pol\": %d", r->fdot, r->cn0,
		          r->selftest_sync ? "true" : "false", r->pol);
	jb_printf(&b, ", \"ebn0\": %.2f, \"sync\": %.2f, \"ok\": %s, \"hex\": \"%s\", \"long\": %s", r->ebn0, r->sync,
	          r->frame.ok ? "true" : "false", hex, r->frame.long_frame ? "true" : "false");
	jb_printf(&b, ", \"cost\": %.2f, \"flips\": [%d, %d], \"reasons\": ", isfinite(r->frame.bch1_cost) ? r->frame.bch1_cost + r->frame.bch2_cost : -1.0,
	          r->frame.bch1_flips, r->frame.bch2_flips);
	jb_str(&b, r->frame.reasons);
	jb_printf(&b, ", \"accept\": \"%s\"", tier(r->accept));
	if (t0_ms > 0) {
		double ms = t0_ms + r->t_cw * 1000;
		time_t s = (time_t)(ms / 1000);
		struct tm tm;
		gmtime_r(&s, &tm);
		char ts[40];
		strftime(ts, sizeof ts, "%Y-%m-%dT%H:%M:%S", &tm);
		jb_printf(&b, ", \"utc\": \"%s.%03d+00:00\"", ts, (int)fmod(ms, 1000));
	}
	if (with_llr && !r->combined) {
		jb_raw(&b, ", \"llr\": [");
		for (int i = 0; i < FGB_BITS; i++) jb_printf(&b, "%s%.1f", i ? ", " : "", r->llr[i]);
		jb_raw(&b, "]");
	}
	jb_raw(&b, "}");
	out_emit(b.s);
	jb_free(&b);
}

static void emit_sgb(const sgb_rec_t *r, double t0_ms) {
	jbuf_t b = {0};
	char hex[64] = "";
	if (r->ok || r->reasons[0]) sgb_hex(r->bits, hex);
	jb_printf(&b, "{\"kind\": \"sgb\", \"t\": %.6f, \"band\": \"%s\", \"f_hz\": %.1f, \"relays\": [", r->t, r->band, r->f_hz);
	for (int k = 0; k < r->nrelay; k++)
		jb_printf(&b, "%s{\"f_hz\": %.1f, \"dt_ms\": %.2f, \"ebn0\": %.2f, \"pre_z\": %.1f, \"alone\": %s}", k ? ", " : "",
		          r->rf[k], r->rdt_ms[k], r->reb[k], r->rpz[k], r->ralone[k] ? "true" : "false");
	jb_printf(&b, "], \"ok\": %s, \"how\": \"", r->ok ? "true" : "false");
	if (r->hard) jb_printf(&b, "hard %d errors", r->nerr);
	else if (r->flips >= 0) jb_printf(&b, "chase %d flips", r->flips);
	else jb_raw(&b, "fail");
	jb_printf(&b, "\", \"hex\": \"%s\", \"accept\": \"%s\", \"reasons\": ", hex, tier(r->accept));
	jb_str(&b, r->reasons);
	if (t0_ms > 0) {
		double ms = t0_ms + r->t * 1000;
		time_t sec = (time_t)(ms / 1000);
		struct tm tm;
		gmtime_r(&sec, &tm);
		char ts[40];
		strftime(ts, sizeof ts, "%Y-%m-%dT%H:%M:%S", &tm);
		jb_printf(&b, ", \"utc\": \"%s.%03d+00:00\"", ts, (int)fmod(ms, 1000));
	}
	jb_raw(&b, "}");
	out_emit(b.s);
	jb_free(&b);
}

int main(int argc, char **argv) {
	const char *inspec = NULL, *fmt = NULL, *gain = NULL, *bands_s = "galileo,glonass";
	double center = 1544.5e6, rate = 2.5e6, start = 0, dur = -1, thr = 13;
	int L = 10, combine = 1, bias = 0, quiet = 0, with_llr = 1, do_sgb = 1, sgb_threads = 4;
	double sgb_span = 8000, sgb_thr = 19, sgb_prez = 7;
	static struct option lo[] = {
		{"start", 1, 0, 1}, {"dur", 1, 0, 2}, {"bands", 1, 0, 3}, {"thr", 1, 0, 4}, {"L", 1, 0, 5},
		{"no-combine", 0, 0, 6}, {"no-llr", 0, 0, 7}, {"no-sgb", 0, 0, 8}, {"sgb-threads", 1, 0, 9},
		{"sgb-span", 1, 0, 10}, {"sgb-thr", 1, 0, 11}, {"sgb-prez", 1, 0, 12}, {"help", 0, 0, 'h'}, {0, 0, 0, 0}};
	int o;
	while ((o = getopt_long(argc, argv, "i:F:c:s:g:bu:qh", lo, NULL)) != -1) {
		switch (o) {
		case 'i': inspec = optarg; break;
		case 'F': fmt = optarg; break;
		case 'c': center = atof(optarg); break;
		case 's': rate = atof(optarg); break;
		case 'g': gain = optarg; break;
		case 'b': bias = 1; break;
		case 'u': if (out_add_udp(optarg)) return 1; break;
		case 'q': quiet = 1; break;
		case 1: start = atof(optarg); break;
		case 2: dur = atof(optarg); break;
		case 3: bands_s = optarg; break;
		case 4: thr = atof(optarg); break;
		case 5: L = atoi(optarg); break;
		case 6: combine = 0; break;
		case 7: with_llr = 0; break;
		case 8: do_sgb = 0; break;
		case 9: sgb_threads = atoi(optarg); break;
		case 10: sgb_span = atof(optarg); break;
		case 11: sgb_thr = atof(optarg); break;
		case 12: sgb_prez = atof(optarg); break;
		default: usage(); return o == 'h' ? 0 : 1;
		}
	}
	if (!inspec) { usage(); return 1; }
	fgb_init();
	input_t *in;
	if (input_is_live(inspec)) {
		live_opts_t lo2 = {.center_hz = center, .rate = rate, .half_span = 0.5e6, .gain = gain, .bias = bias};
		in = input_open_live(inspec, &lo2);
	} else {
		in = input_open_file(inspec, fmt ? iq_format_parse(fmt) : -1);
	}
	if (!in) return 1;
	if (!in->rate) in->rate = rate;
	if (!in->center_hz) in->center_hz = center;
	double fs = in->rate, fc = in->center_hz, f1 = fs / D1;
	out_init(!quiet);
	signal(SIGINT, on_sig);
	signal(SIGTERM, on_sig);

	band_t bands[3];
	int nb = 0;
	char tmp[64];
	snprintf(tmp, sizeof tmp, "%s", bands_s);
	for (char *tok = strtok(tmp, ","); tok && nb < 3; tok = strtok(NULL, ",")) {
		band_t *b = &bands[nb];
		memset(b, 0, sizeof *b);
		if (!strcmp(tok, "galileo")) { b->name = "galileo"; b->lo = 1544.04e6; b->hi = 1544.16e6; }
		else if (!strcmp(tok, "glonass")) { b->name = "glonass"; b->lo = 1544.84e6; b->hi = 1544.96e6; }
		else if (!strcmp(tok, "quiet")) { b->name = "quiet"; b->lo = 1543.62e6; b->hi = 1543.88e6; }
		else { fprintf(stderr, "unknown band %s\n", tok); return 1; }
		if (b->lo < fc - fs / 2 || b->hi > fc + fs / 2) { fprintf(stderr, "band %s outside the input\n", tok); continue; }
		b->mid = (b->lo + b->hi) / 2;
		b->cs = chan_stream_new(fs, b->mid - fc, D1, (b->hi - b->lo) * 1.05);
		nb++;
	}
	fprintf(stderr, "meosar-rx: %s, %.0f S/s at %.4f MHz, %d band(s), threshold %.1f dB, Chase L=%d%s%s\n",
	        in->name, fs, fc / 1e6, nb, thr, L, combine ? ", relay combining" : "", do_sgb ? ", SGB" : "");
	sgb_acq_t *acq = do_sgb ? sgb_acq_new(f1, sgb_span, sgb_threads) : NULL;
	sgb_rec_t *srecs = malloc(sizeof(sgb_rec_t) * 256);
	long long n_sgb = 0;

	// skip to --start
	const int CH = 65536;
	float *iq = malloc(sizeof(float) * 2 * CH);
	cf_t *xin = malloc(sizeof(cf_t) * CH), *yout = malloc(sizeof(cf_t) * (CH / D1 + 16));
	long long skip = (long long)(start * fs), done_in = 0;
	while (skip > 0 && !stop) {
		int k = in->read(in, iq, skip > CH ? CH : (int)skip);
		if (k <= 0) break;
		skip -= k;
	}
	double t_block = start, blk = 8.0, t_end = dur > 0 ? start + dur : INFINITY;
	long long need_end = (long long)((t_block + blk + POST - start) * f1);   // band index (rel. to start)
	int cap = 4096;
	fgb_rec_t *recs = malloc(sizeof(fgb_rec_t) * cap);
	cand_t *cands = malloc(sizeof(cand_t) * 4096);
	long long n_cand = 0, n_acc = 0;
	clock_t tic = clock();
	int eof = 0;
	while (!stop && t_block < t_end - 1.0) {
		// fill until every band covers the block end
		while (!eof && !stop && (nb == 0 || bands[0].base + (long long)bands[0].len < need_end)) {
			long long want = isfinite(t_end) ? (long long)((t_end - start) * fs) - done_in : CH;
			if (want <= 0) { eof = 1; break; }
			int k = in->read(in, iq, want < CH ? (int)want : CH);
			if (k <= 0) { eof = 1; break; }
			done_in += k;
			for (int i = 0; i < k; i++) xin[i] = iq[2 * i] + iq[2 * i + 1] * I;
			for (int bi = 0; bi < nb; bi++) {
				band_t *b = &bands[bi];
				size_t w = chan_stream_push(b->cs, xin, (size_t)k, yout, CH / D1 + 16);
				if (b->len + w > b->cap) { b->cap = (b->len + w) * 2; b->buf = realloc(b->buf, sizeof(cf_t) * b->cap); }
				memcpy(b->buf + b->len, yout, sizeof(cf_t) * w);
				b->len += w;
			}
		}
		double t_stop = t_block + blk;
		if (eof) t_stop = start + (double)(bands[0].base + (long long)bands[0].len) / f1 - POST;
		if (t_stop > t_end) t_stop = t_end;
		if (t_stop <= t_block) break;
		int nrec = 0;
		for (int bi = 0; bi < nb; bi++) {
			band_t *b = &bands[bi];
			double t_y0 = start + (double)b->base / f1;
			int nc = fgb_detect(b->buf, b->len, f1, b->mid, b->lo, b->hi, t_y0, thr, cands, 4096);
			int first = nrec;
			for (int c = 0; c < nc && nrec < cap - 64; c++) {
				if (cands[c].t < t_block || cands[c].t >= t_stop) continue;   // owned by another block
				fgb_rec_t *r = &recs[nrec];
				if (!fgb_process_candidate(b->buf, b->len, f1, b->mid, t_y0, cands[c].t, cands[c].f, L, r)) continue;
				r->det_snr = cands[c].snr_db;
				r->band = b->name;
				nrec++;
			}
			n_cand += nrec - first;
			if (combine) {
				// combine within this band only
				int nn = fgb_combine_relays(recs + first, nrec - first, cap - first, L);
				for (int k = nrec - first; k < nn; k++)
					for (int j = 0; j < recs[first + k].combined; j++) recs[first + k].relay_idx[j] += 0;
				nrec = first + nn;
			}
			for (int k = first; k < nrec; k++) {
				emit(&recs[k], recs + first, in->start_ms, with_llr);
				n_acc += recs[k].accept > 0;
			}
			if (acq) {
				int ns = sgb_process_block(acq, b->buf, b->len, f1, b->mid, t_y0, t_block, t_stop, sgb_thr, sgb_prez, srecs, 256);
				for (int k = 0; k < ns; k++) {
					srecs[k].band = b->name;
					emit_sgb(&srecs[k], in->start_ms);
					n_sgb += srecs[k].accept > 0;
				}
			}
		}
		// advance and trim
		t_block = t_stop;
		need_end = (long long)((t_block + blk + POST - start) * f1);
		long long keep_from = (long long)((t_block - PRE - start) * f1);
		for (int bi = 0; bi < nb; bi++) {
			band_t *b = &bands[bi];
			if (keep_from > b->base) {
				size_t drop = (size_t)(keep_from - b->base);
				if (drop > b->len) drop = b->len;
				memmove(b->buf, b->buf + drop, sizeof(cf_t) * (b->len - drop));
				b->len -= drop;
				b->base += (long long)drop;
			}
		}
		double cpu = (double)(clock() - tic) / CLOCKS_PER_SEC;
		fprintf(stderr, "[%7.1f s] FGB candidates %lld accepted %lld, SGB accepted %lld  cpu %.2fx real time%s\n", t_block, n_cand, n_acc, n_sgb,
		        cpu > 0 ? (t_block - start) / cpu : 0, input_dropped ? "  (dropped samples)" : "");
		if (eof) break;
	}
	in->close(in);
	out_close();
	return 0;
}
