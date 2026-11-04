// SPDX-License-Identifier: GPL-3.0-or-later
#include "fgb_demod.h"

#include <fftw3.h>
#define RE(a, i) (((float *)(a))[2 * (i)])
#define IM(a, i) (((float *)(a))[2 * (i) + 1])
#include <math.h>
#include <stdlib.h>
#include <string.h>

static const uint8_t KNOWN_N[24] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 0,0,0,1,0,1,1,1,1};
static const uint8_t KNOWN_S[24] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 0,1,1,0,1,0,0,0,0};

// ---- FFT helper with a small plan cache ----
typedef struct { int n; fftwf_plan p; fftwf_complex *in, *out; } fft_t;
static fft_t fft_cache[8];
static fft_t *fft_get(int n) {
	for (int i = 0; i < 8; i++)
		if (fft_cache[i].n == n) return &fft_cache[i];
	for (int i = 0; i < 8; i++)
		if (!fft_cache[i].n) {
			fft_t *f = &fft_cache[i];
			f->n = n;
			f->in = fftwf_malloc(sizeof(fftwf_complex) * n);
			f->out = fftwf_malloc(sizeof(fftwf_complex) * n);
			f->p = fftwf_plan_dft_1d(n, f->in, f->out, FFTW_FORWARD, FFTW_ESTIMATE);
			return f;
		}
	return NULL;
}

static int next_pow2(int n) { int p = 1; while (p < n) p <<= 1; return p; }

// Manchester waveform of nb bits at spb samples/bit into m[0..n)
static void manchester(const uint8_t *bits, int nb, double spb, int n, float *m) {
	for (int i = 0; i < n; i++) {
		double t = i / spb;
		int k = (int)t;
		if (k > nb - 1) k = nb - 1;
		float b = bits[k] ? 1.f : -1.f;
		m[i] = (t - k) < 0.5 ? b : -b;
	}
}

static double score_fd(const cf_t *x, int n, double fs, double f, double d) {
	double re = 0, im = 0;
	for (int i = 0; i < n; i++) {
		double t = i / fs, ph = -2 * M_PI * (f * t + 0.5 * d * t * t);
		double c = cos(ph), s = sin(ph);
		re += crealf(x[i]) * c - cimagf(x[i]) * s;
		im += crealf(x[i]) * s + cimagf(x[i]) * c;
	}
	return hypot(re, im);
}

static void estimate_carrier(const cf_t *x, int n, double fs, double *f_out, double *fd_out) {
	int nfft = next_pow2(n * 4);
	fft_t *F = fft_get(nfft);
	double bestp = -1, f0 = 0, fd0 = 0;
	for (double fd = -12.0; fd <= 12.01; fd += 1.5) {
		for (int i = 0; i < nfft; i++) {
			if (i < n) {
				double t = i / fs, ph = -M_PI * fd * t * t;
				cf_t v = x[i] * (float)cos(ph) + x[i] * (float)sin(ph) * I;
				RE(F->in, i) = crealf(v); IM(F->in, i) = cimagf(v);
			} else RE(F->in, i) = IM(F->in, i) = 0;
		}
		fftwf_execute(F->p);
		for (int k = 0; k < nfft; k++) {
			double fr = (k < nfft / 2 ? k : k - nfft) * fs / nfft;
			if (fabs(fr) > 80) continue;
			double p = RE(F->out, k) * RE(F->out, k) + IM(F->out, k) * IM(F->out, k);
			if (p > bestp) { bestp = p; f0 = fr; fd0 = fd; }
		}
	}
	double sf = fs / nfft, sd = 0.75, cur = score_fd(x, n, fs, f0, fd0);
	for (int it = 0; it < 25; it++) {
		int moved = 0;
		double cand[4][2] = {{sf, 0}, {-sf, 0}, {0, sd}, {0, -sd}};
		for (int c = 0; c < 4; c++) {
			double s = score_fd(x, n, fs, f0 + cand[c][0], fd0 + cand[c][1]);
			if (s > cur) { cur = s; f0 += cand[c][0]; fd0 += cand[c][1]; moved = 1; }
		}
		if (!moved) { sf /= 2; sd /= 2; if (sf < 0.02) break; }
	}
	*f_out = f0;
	*fd_out = fd0;
}

static void bit_stats(const cf_t *z, int n, double spb, int nstart, int pol, double *r) {
	for (int b = 0; b < FGB_BITS; b++) {
		double a = nstart + b * spb;
		int h0 = (int)lround(a), h1 = (int)lround(a + spb / 2), h2 = (int)lround(a + spb);
		r[b] = 0;
		if (h2 > n) continue;
		double s = 0;
		for (int i = h0; i < h1; i++) s += cimagf(z[i]);
		for (int i = h1; i < h2; i++) s -= cimagf(z[i]);
		r[b] = s * pol;
	}
}

// centred moving average with reflect boundary (as scipy.ndimage.uniform_filter1d)
static void movavg(const double *x, int n, int w, double *y) {
	int h = w / 2;
	// prefix sums over the reflect-padded signal: index j in [-h, n-1+w-1-h]
	int lo = -h, hi = n - 1 + (w - 1 - h), len = hi - lo + 1;
	double *ps = malloc(sizeof(double) * (len + 1));
	ps[0] = 0;
	for (int j = lo; j <= hi; j++) {
		int k = j;
		while (k < 0 || k >= n) k = k < 0 ? -k - 1 : 2 * n - k - 1;
		ps[j - lo + 1] = ps[j - lo] + x[k];
	}
	for (int i = 0; i < n; i++) {
		int a = i - h - lo;            // window [i-h, i-h+w-1]
		y[i] = (ps[a + w] - ps[a]) / w;
	}
	free(ps);
}

void fgb_demodulate(const cf_t *x, size_t nn, double fs, double t_cw_guess, fgb_demod_t *d) {
	memset(d, 0, sizeof *d);
	int n = (int)nn;
	const double search_s = 0.08, smooth_s = 0.030;
	double spb = fs / FGB_BIT_RATE;
	int i0 = (int)((t_cw_guess - search_s) * fs); if (i0 < 0) i0 = 0;
	int i1 = (int)((t_cw_guess + FGB_CW_SECONDS + FGB_BITS / FGB_BIT_RATE + search_s) * fs); if (i1 > n) i1 = n;
	double f0, fd;
	estimate_carrier(x + i0, i1 - i0, fs, &f0, &fd);
	d->f_hz = f0; d->fdot = fd;
	cf_t *z = malloc(sizeof(cf_t) * n);
	double sre = 0, sim = 0;
	for (int i = 0; i < n; i++) {
		double t = (i - i0) / fs, ph = -2 * M_PI * (f0 * t + 0.5 * fd * t * t);
		z[i] = x[i] * (float)cos(ph) + x[i] * (float)sin(ph) * I;
		if (i >= i0 && i < i1) { sre += crealf(z[i]); sim += cimagf(z[i]); }
	}
	double a0 = -atan2(sim, sre);
	cf_t rot0 = (float)cos(a0) + (float)sin(a0) * I;
	for (int i = 0; i < n; i++) z[i] *= rot0;

	// timing on the 24 known bits, both sync words and polarities
	int nref = (int)lround(24 * spb);
	float *refn = malloc(sizeof(float) * nref), *refs = malloc(sizeof(float) * nref);
	manchester(KNOWN_N, 24, spb, nref, refn);
	manchester(KNOWN_S, 24, spb, nref, refs);
	int nd = (int)lround(FGB_BITS * spb);
	int lo = (int)((t_cw_guess + FGB_CW_SECONDS - search_s) * fs); if (lo < 0) lo = 0;
	int hi = (int)((t_cw_guess + FGB_CW_SECONDS + search_s) * fs);
	if (hi > n - nref - nd / 2) hi = n - nref - nd / 2;
	if (hi <= lo) { free(z); free(refn); free(refs); return; }
	double best = -1, bestc = 0; int bestk = 0, st = 0;
	for (int k = lo; k <= hi; k++) {
		double cn = 0, cs = 0;
		for (int i = 0; i < nref; i++) { cn += cimagf(z[k + i]) * refn[i]; cs += cimagf(z[k + i]) * refs[i]; }
		if (fabs(cn) > best) { best = fabs(cn); bestc = cn; bestk = k; st = 0; }
		if (fabs(cs) > best) { best = fabs(cs); bestc = cs; bestk = k; st = 1; }
	}
	free(refn); free(refs);
	int pol = bestc > 0 ? 1 : -1, n_data = bestk;
	d->polarity = pol; d->selftest_sync = st;
	const uint8_t *known = st ? KNOWN_S : KNOWN_N;
	int nt = nd; if (n_data + nt > n) nt = n - n_data;

	// decision-directed phase smoothing
	int ncw0 = n_data - (int)(FGB_CW_SECONDS * fs); if (ncw0 < 0) ncw0 = 0;
	int w = (int)(smooth_s * fs); if (w < 3) w = 3; w |= 1;
	double r[FGB_BITS];
	float *m = calloc((size_t)n, sizeof(float)), *mm = malloc(sizeof(float) * nd);
	double *cre = malloc(sizeof(double) * n), *cim = malloc(sizeof(double) * n);
	double *sre_ = malloc(sizeof(double) * n), *sim_ = malloc(sizeof(double) * n);
	int last = n_data + nt - 1;
	for (int it = 0; it < 3; it++) {
		bit_stats(z, n, spb, n_data, pol, r);
		uint8_t hard[FGB_BITS];
		for (int b = 0; b < FGB_BITS; b++) hard[b] = r[b] > 0;
		memcpy(hard, known, 24);
		manchester(hard, FGB_BITS, spb, nd, mm);
		memset(m, 0, sizeof(float) * n);
		for (int i = 0; i < nt; i++) m[n_data + i] = mm[i] * pol;
		for (int i = 0; i < n; i++) {
			if (i >= ncw0 && i <= last) {
				double ph = -FGB_MOD_INDEX * m[i];
				cf_t c = z[i] * ((float)cos(ph) + (float)sin(ph) * I);
				cre[i] = crealf(c); cim[i] = cimagf(c);
			} else cre[i] = cim[i] = 0;
		}
		movavg(cre, n, w, sre_);
		movavg(cim, n, w, sim_);
		double ph_lo = atan2(sim_[ncw0], sre_[ncw0]), ph_hi = atan2(sim_[last], sre_[last]);
		for (int i = 0; i < n; i++) {
			double ph = i < ncw0 ? ph_lo : i > last ? ph_hi : atan2(sim_[i], sre_[i]);
			z[i] *= (float)cos(-ph) + (float)sin(-ph) * I;
		}
	}
	bit_stats(z, n, spb, n_data, pol, r);
	// amplitude from the CW, noise from before it
	int c0 = ncw0 + (int)(0.01 * fs), c1 = n_data - (int)(0.01 * fs);
	double A = 0;
	if (c1 - c0 > 10) { for (int i = c0; i < c1; i++) A += crealf(z[i]); A /= (c1 - c0); }
	int nz = ncw0 - (int)(0.02 * fs);
	if (nz < 4 * spb) { free(z); free(m); free(mm); free(cre); free(cim); free(sre_); free(sim_); return; }
	int nsb = (int)((nz - spb) / spb);
	double sig2 = 0;
	for (int k = 0; k < nsb; k++) {
		int a = (int)(k * spb), b = (int)(k * spb + spb / 2), e = (int)((k + 1) * spb);
		double s = 0;
		for (int i = a; i < b; i++) s += cimagf(x[i]);
		for (int i = b; i < e; i++) s -= cimagf(x[i]);
		sig2 += s * s;
	}
	sig2 /= nsb;
	double mu = A * sin(FGB_MOD_INDEX) * spb;
	for (int b = 0; b < FGB_BITS; b++) d->llr[b] = 2 * mu * r[b] / sig2;
	d->ebn0_db = 10 * log10(fmax(mu * mu / (2 * sig2), 1e-12));
	d->cn0_dbhz = d->ebn0_db + 10 * log10(FGB_BIT_RATE) - 10 * log10(sin(FGB_MOD_INDEX) * sin(FGB_MOD_INDEX));
	d->t_data = n_data / fs;
	d->t_cw = d->t_data - FGB_CW_SECONDS;
	double ks = 0;
	for (int b = 0; b < 24; b++) ks += r[b] * (2 * known[b] - 1);
	d->sync_score = ks / 24 / sqrt(sig2);
	d->ok = 1;
	free(z); free(m); free(mm); free(cre); free(cim); free(sre_); free(sim_);
}

#define D3 6

int fgb_process_candidate(const cf_t *y, size_t n, double f1, double mid, double t_y0,
                          double t_frame, double f_abs, int L, fgb_rec_t *r) {
	memset(r, 0, sizeof *r);
	r->t_frame = t_frame; r->f_det = f_abs;
	double t_cw_guess = t_frame + 0.01, t_win = t_cw_guess - 0.25;
	long i0 = lround((t_win - t_y0) * f1), len = (long)(1.05 * f1);
	if (i0 < 0 || i0 + len > (long)n) return 0;
	const cf_t *seg = y + i0;
	double f2 = f1 / D3;
	size_t m = ((size_t)len + D3 - 1) / D3;
	cf_t *xb = malloc(sizeof(cf_t) * m);
	channelise(seg, (size_t)len, f_abs - mid, f1, D3, 8000, 0, xb);
	// coarse frequency within one detector bin, 140 ms after the guess
	int a = (int)(0.25 * f2), b = (int)(0.39 * f2), N = 1 << 15;
	fft_t *F = fft_get(N);
	for (int i = 0; i < N; i++) {
		if (i < b - a) { RE(F->in, i) = crealf(xb[a + i]); IM(F->in, i) = cimagf(xb[a + i]); }
		else RE(F->in, i) = IM(F->in, i) = 0;
	}
	fftwf_execute(F->p);
	double bp = -1, ff = 0;
	for (int k = 0; k < N; k++) {
		double fr = (k < N / 2 ? k : k - N) * f2 / N;
		if (fabs(fr) >= 40) continue;
		double p = RE(F->out, k) * RE(F->out, k) + IM(F->out, k) * IM(F->out, k);
		if (p > bp) { bp = p; ff = fr; }
	}
	channelise(seg, (size_t)len, f_abs - mid + ff, f1, D3, 8000, 0, xb);
	fgb_demod_t d;
	fgb_demodulate(xb, m, f2, 0.25, &d);
	free(xb);
	if (!d.ok) return 0;
	fgb_decode_soft(d.llr, L, &r->frame);
	r->valid = 1;
	r->t_cw = t_win + d.t_cw;
	r->f_hz = f_abs + ff + d.f_hz;
	r->fdot = d.fdot; r->ebn0 = d.ebn0_db; r->cn0 = d.cn0_dbhz; r->sync = d.sync_score;
	r->selftest_sync = d.selftest_sync; r->pol = d.polarity;
	memcpy(r->llr, d.llr, sizeof d.llr);
	r->accept = fgb_accept(r);
	return 1;
}

int fgb_accept(const fgb_rec_t *r) {
	if (!r->frame.ok) return 0;
	if (r->ebn0 >= 1.5 && r->sync >= 1.3 && r->frame.bch1_flips <= 3) return 2;   // ID from BCH-1; PDF-2 flagged separately
	if (r->ebn0 >= 0.0 && r->sync >= 1.0) return 1;
	return 0;
}

static int by_tcw(const void *a, const void *b) {
	double d = (*(fgb_rec_t *const *)a)->t_cw - (*(fgb_rec_t *const *)b)->t_cw;
	return (d > 0) - (d < 0);
}

int fgb_combine_relays(fgb_rec_t *recs, int n, int cap, int L) {
	fgb_rec_t **pool = malloc(sizeof(*pool) * (n + 1));
	int np = 0;
	for (int i = 0; i < n; i++)
		if (recs[i].valid && !recs[i].combined && recs[i].sync >= 0.8) pool[np++] = &recs[i];
	qsort(pool, (size_t)np, sizeof(*pool), by_tcw);
	char *used = calloc((size_t)np + 1, 1);
	int total = n;
	for (int i = 0; i < np && total < cap; i++) {
		if (used[i]) continue;
		int grp[8], ng = 0;
		grp[ng++] = i;
		for (int j = i + 1; j < np && ng < 8; j++) {
			if (pool[j]->t_cw - pool[i]->t_cw > 0.040) break;
			if (used[j]) continue;
			int same = 0;
			for (int k = 0; k < ng; k++)
				if (fabs(pool[j]->f_hz - pool[grp[k]]->f_hz) < 20) same = 1;
			if (same) { used[j] = 1; continue; }
			grp[ng++] = j;
		}
		if (ng < 2) continue;
		fgb_rec_t *c = &recs[total];
		memset(c, 0, sizeof *c);
		double lin = 0, sy = 0;
		for (int k = 0; k < ng; k++) {
			fgb_rec_t *q = pool[grp[k]];
			used[grp[k]] = 1;
			for (int b = 0; b < FGB_BITS; b++) c->llr[b] += q->llr[b];
			lin += pow(10, q->ebn0 / 10);
			sy += q->sync * q->sync;
			c->relay_idx[k] = (int)(q - recs);
		}
		fgb_rec_t *q0 = pool[grp[0]];
		c->valid = 1; c->combined = ng; c->band = q0->band;
		c->t_frame = q0->t_frame; c->t_cw = q0->t_cw; c->f_hz = q0->f_hz;
		c->ebn0 = 10 * log10(lin); c->sync = sqrt(sy);
		fgb_decode_soft(c->llr, L, &c->frame);
		c->accept = fgb_accept(c);
		total++;
	}
	free(pool); free(used);
	return total;
}
