// SPDX-License-Identifier: GPL-3.0-or-later
#include "fgb_repeat.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define SEG0 24                 // BCH-1 codeword: bits 24..105 (0-based)
#define SEGN 82
#define Z_KNOWN 6.0
#define Z_PAIR 5.0              // a pair links this strongly: two bursts may decode together
#define Z_LINK 4.0              // weaker links count only in groups of three bursts or more
#define Z_MEMBER 3.0
#define LODDS_MIN 6.9          // ln 1000: a match, or a cluster decode, against every other known codeword
#define LLR_CLIP 6.0            // per bit, in log-likelihood ratios: on strict decodes, bits with |LLR|
                                // 6-9 are wrong 0.3% of the time (worth ~5.8, not 6-9)
#define SAME_BURST_S 1.5        // closer than this: the same transmission
#define MAX_KNOWN 2048
#define MAX_WEAK 1024
#define MAX_POOL 4096
#define KNOWN_FORGET_S 86400.0

typedef struct {
	char id[16];
	float s[SEGN];              // codeword as +-1
	fgb_frame_t frame;
	double t_last;
} known_t;

typedef struct {
	char id[16];
	fgb_frame_t frame;
	double t;
} weak_t;

typedef struct {
	double llr[FGB_BITS];
	float v[SEGN];              // codeword part, normalised to unit length
	double t_cw, f_hz, ebn0;
	const char *band;
} unit_t;

struct fgb_repeat {
	double window;
	int L;
	known_t *known; int nk;
	weak_t *weak; int nw;
	unit_t *pool; int np;
	double t_now;
};

fgb_repeat_t *fgb_repeat_new(double window_s, int L) {
	fgb_repeat_t *R = calloc(1, sizeof *R);
	R->window = window_s; R->L = L;
	R->known = calloc(MAX_KNOWN, sizeof(known_t));
	R->weak = calloc(MAX_WEAK, sizeof(weak_t));
	R->pool = calloc(MAX_POOL, sizeof(unit_t));
	return R;
}

void fgb_repeat_free(fgb_repeat_t *R) {
	if (!R) return;
	free(R->known); free(R->weak); free(R->pool); free(R);
}

static void id_of(const fgb_frame_t *f, char id[16]) {
	char hex[40];
	fgb_hex(f->bits, f->nbits, hex);
	memcpy(id, hex, 15);
	id[15] = 0;
}

// Correlation of soft values with a +-1 codeword, normalised: ~N(0,1) when unrelated.
static double zscore(const double *llr, const float *s) {
	double num = 0, den = 0;
	for (int i = 0; i < SEGN; i++) { num += llr[SEG0 + i] * s[i]; den += llr[SEG0 + i] * llr[SEG0 + i]; }
	return den > 0 ? num / sqrt(den) : 0;
}

static void codeword(const fgb_frame_t *f, float *s) {
	for (int i = 0; i < SEGN; i++) s[i] = f->bits[SEG0 + i] ? 1.0f : -1.0f;
}

static void known_add(fgb_repeat_t *R, const fgb_frame_t *f, double t) {
	char id[16];
	id_of(f, id);
	for (int i = 0; i < R->nk; i++)
		if (!strcmp(R->known[i].id, id)) { R->known[i].frame = *f; codeword(f, R->known[i].s); R->known[i].t_last = t; return; }
	int slot = R->nk;
	if (R->nk == MAX_KNOWN) {   // full: replace the beacon heard longest ago
		slot = 0;
		for (int i = 1; i < R->nk; i++) if (R->known[i].t_last < R->known[slot].t_last) slot = i;
	} else {
		R->nk++;
	}
	known_t *k = &R->known[slot];
	memcpy(k->id, id, 16);
	k->frame = *f; codeword(f, k->s); k->t_last = t;
}

// A weak decode confirms a beacon when the same ID was also decoded in a different burst.
static void weak_add(fgb_repeat_t *R, const fgb_frame_t *f, double t) {
	char id[16];
	id_of(f, id);
	for (int i = 0; i < R->nw; i++) {
		if (strcmp(R->weak[i].id, id)) continue;
		if (fabs(t - R->weak[i].t) > SAME_BURST_S) {
			known_add(R, f, t);
			R->weak[i] = R->weak[--R->nw];
		}
		return;
	}
	if (R->nw == MAX_WEAK) memmove(R->weak, R->weak + 1, sizeof(weak_t) * (--R->nw));
	weak_t *w = &R->weak[R->nw++];
	memcpy(w->id, id, 16); w->frame = *f; w->t = t;
}

static void prune(fgb_repeat_t *R) {
	int j = 0;
	for (int i = 0; i < R->np; i++)
		if (R->t_now - R->pool[i].t_cw <= R->window) { if (j != i) R->pool[j] = R->pool[i]; j++; }
	R->np = j;
	j = 0;
	for (int i = 0; i < R->nw; i++)
		if (R->t_now - R->weak[i].t <= R->window) { if (j != i) R->weak[j] = R->weak[i]; j++; }
	R->nw = j;
	j = 0;
	for (int i = 0; i < R->nk; i++)
		if (R->t_now - R->known[i].t_last <= KNOWN_FORGET_S) { if (j != i) R->known[j] = R->known[i]; j++; }
	R->nk = j;
}

static void unit_from(const fgb_rec_t *r, unit_t *u) {
	memcpy(u->llr, r->llr, sizeof u->llr);
	double n = 0;
	for (int i = 0; i < SEGN; i++) n += r->llr[SEG0 + i] * r->llr[SEG0 + i];
	n = n > 0 ? sqrt(n) : 1;
	for (int i = 0; i < SEGN; i++) u->v[i] = (float)(r->llr[SEG0 + i] / n);
	u->t_cw = r->t_cw; u->f_hz = r->f_hz; u->ebn0 = r->ebn0; u->band = r->band;
}

static double pair_z(const unit_t *a, const unit_t *b) {
	double d = 0;
	for (int i = 0; i < SEGN; i++) d += (double)a->v[i] * b->v[i];
	return d * sqrt((double)SEGN);
}

// Log-likelihood ratio of codeword s against the likeliest other known codeword, from calibrated
// LLRs (2 mu r / sigma^2): sum over the bits where they differ of llr * s, each bit clipped at
// LLR_CLIP (the demodulator is overconfident beyond it). +1e9 when there is none.
// Unlike z, which is capped at sqrt(82) (a perfect burst of one beacon still scores 7.5 against a
// sibling 7 bits away), it keeps growing with the evidence. Beacons of one family (consecutive
// serials) have codewords only 7-12 bits apart.
static double lodds_vs_known(const fgb_repeat_t *R, const double *llr, const float *s) {
	double worst = 1e9;
	for (int k = 0; k < R->nk; k++) {
		const float *o = R->known[k].s;
		double lo = 0;
		int d = 0;
		for (int i = 0; i < SEGN; i++)
			if (o[i] != s[i]) { lo += fmax(-LLR_CLIP, fmin(LLR_CLIP, llr[SEG0 + i])) * s[i]; d++; }
		if (d && lo < worst) worst = lo;
	}
	return worst;
}

// Try to decode the new unit u together with earlier unmatched units. On success fills *out,
// removes the members from the pool and returns 1.
static int try_cluster(fgb_repeat_t *R, const unit_t *u, fgb_rep_t *out) {
	int idx[FGB_REP_MAXMEM], ni = 0, strong = 0;
	for (int i = 0; i < R->np && ni < FGB_REP_MAXMEM - 1; i++) {
		if (fabs(R->pool[i].t_cw - u->t_cw) <= SAME_BURST_S) continue;
		double z = pair_z(u, &R->pool[i]);
		if (z >= Z_LINK) { idx[ni++] = i; strong |= z >= Z_PAIR; }
	}
	if (!ni) return 0;
	// one level of expansion: units that match a member strongly, though not u itself
	for (int a = 0, na = ni; a < na && ni < FGB_REP_MAXMEM - 1; a++)
		for (int i = 0; i < R->np && ni < FGB_REP_MAXMEM - 1; i++) {
			int have = 0;
			for (int k = 0; k < ni; k++) if (idx[k] == i) have = 1;
			if (!have && fabs(R->pool[i].t_cw - u->t_cw) > SAME_BURST_S && pair_z(&R->pool[idx[a]], &R->pool[i]) >= Z_PAIR)
				idx[ni++] = i;
		}
	fgb_frame_t f;
	float s[SEGN];
	double cl_lodds = 0;
	for (;;) {
		// Noise pairs reach Z_LINK now and then (p ~ 3e-5 each); two of them on one burst are rare.
		if (!ni || (!strong && ni < 2)) return 0;
		double sum[FGB_BITS];
		memcpy(sum, u->llr, sizeof sum);
		for (int k = 0; k < ni; k++)
			for (int b = 0; b < FGB_BITS; b++) sum[b] += R->pool[idx[k]].llr[b];
		fgb_decode_soft(sum, R->L, &f);
		if (!f.ok || f.bch1_flips > 3 || f.pdf2 == 2) return 0;   // PDF-2 doubtful (orbitography has none)
		codeword(&f, s);
		// Every member must match the decoded codeword on its own, and prefer it to every other
		// known codeword (log-odds >= 0; a weak burst cannot do more than that against a sibling).
		// Members that prefer another beacon are left out and the rest tried again; if u itself
		// does, nothing stands.
		if (zscore(u->llr, s) < Z_MEMBER || lodds_vs_known(R, u->llr, s) < 0) return 0;
		int nk = 0;
		for (int k = 0; k < ni; k++) {
			const double *l = R->pool[idx[k]].llr;
			double z = zscore(l, s);
			if (z >= Z_MEMBER && lodds_vs_known(R, l, s) >= 0) idx[nk++] = idx[k];
		}
		if (nk == ni) {
			// The cluster as a whole must favour the decoded codeword over every other known
			// beacon's by odds of 1000:1: a family's bursts summed together decode to one
			// sibling or another.
			cl_lodds = lodds_vs_known(R, sum, s);
			if (cl_lodds < LODDS_MIN) return 0;
			break;
		}
		ni = nk;
		strong = 0;
		for (int k = 0; k < ni; k++) strong |= pair_z(u, &R->pool[idx[k]]) >= Z_PAIR;
	}
	memset(out, 0, sizeof *out);
	out->kind = 2; out->frame = f; out->lodds = fmin(cl_lodds, 999);
	const unit_t *all[FGB_REP_MAXMEM];
	int nall = 0;
	for (int k = 0; k < ni; k++) all[nall++] = &R->pool[idx[k]];
	all[nall++] = u;
	double tdist[FGB_REP_MAXMEM];
	int nd = 0;
	for (int k = 0; k < nall; k++) {
		const unit_t *m = all[k];
		out->mem[out->nmem++] = (fgb_rep_member_t){m->t_cw, m->f_hz, m->ebn0, zscore(m->llr, s),
		                                            fmin(lodds_vs_known(R, m->llr, s), 999), m->band};
		int seen = 0;
		for (int j = 0; j < nd; j++) if (fabs(tdist[j] - m->t_cw) <= SAME_BURST_S) seen = 1;
		if (!seen) tdist[nd++] = m->t_cw;
	}
	out->bursts = nd;
	// remove the members from the pool (highest index first)
	for (int a = 0; a < ni; a++)
		for (int b = a + 1; b < ni; b++)
			if (idx[b] > idx[a]) { int t = idx[a]; idx[a] = idx[b]; idx[b] = t; }
	for (int k = 0; k < ni; k++) R->pool[idx[k]] = R->pool[--R->np];
	return 1;
}

int fgb_repeat_block(fgb_repeat_t *R, const fgb_rec_t *recs, int n, fgb_rep_t *out, int cap) {
	// Which records are relays summed into a combined record (those are covered by it).
	char *inside = calloc((size_t)n + 1, 1);
	for (int i = 0; i < n; i++)
		if (recs[i].combined)
			for (int k = 0; k < recs[i].combined; k++)
				if (recs[i].relay_idx[k] >= 0 && recs[i].relay_idx[k] < n) inside[recs[i].relay_idx[k]] = 1;
	// A burst counts as decoded when its combined record or any of its relays was accepted.
	char *decoded = calloc((size_t)n + 1, 1);
	for (int i = 0; i < n; i++) {
		if (!recs[i].valid) continue;
		if (recs[i].t_cw > R->t_now) R->t_now = recs[i].t_cw;
		int acc = recs[i].accept;
		if (recs[i].combined)
			for (int k = 0; k < recs[i].combined; k++) acc |= recs[recs[i].relay_idx[k]].accept;
		decoded[i] = acc > 0;
		// teach the known table first, so later bursts in this block can match
		if (recs[i].accept == 2) known_add(R, &recs[i].frame, recs[i].t_cw);
		else if (recs[i].accept == 1) weak_add(R, &recs[i].frame, recs[i].t_cw);
	}
	prune(R);
	int nout = 0;
	for (int i = 0; i < n && nout < cap; i++) {
		const fgb_rec_t *r = &recs[i];
		if (!r->valid || inside[i] || decoded[i]) continue;
		// 1. a known beacon? The best correlation, if it is clearly likelier than every other known
		// beacon (two beacons with consecutive serials correlate with each other's bursts almost
		// as well; z2 records the runner-up's correlation).
		int best = -1;
		double bz = -1e9, z2 = -1e9, lo = 0;
		for (int k = 0; k < R->nk; k++) {
			double z = zscore(r->llr, R->known[k].s);
			if (z > bz) { z2 = bz; bz = z; best = k; }
			else if (z > z2) z2 = z;
		}
		if (best >= 0 && bz >= Z_KNOWN) lo = lodds_vs_known(R, r->llr, R->known[best].s);
		// Too weak, or ambiguous between similar codewords: not attributed; it may still cluster.
		if (bz < Z_KNOWN || lo < LODDS_MIN) best = -1;
		if (best >= 0) {
			fgb_rep_t *o = &out[nout++];
			memset(o, 0, sizeof *o);
			o->kind = 1; o->frame = R->known[best].frame; o->z = bz; o->z2 = z2 > -1e8 ? z2 : 0; o->nmem = 1; o->bursts = 1;
			o->lodds = fmin(lo, 999);
			o->mem[0] = (fgb_rep_member_t){r->t_cw, r->f_hz, r->ebn0, bz, o->lodds, r->band};
			R->known[best].t_last = r->t_cw;
			continue;
		}
		// 2. with earlier unmatched bursts?
		if (r->sync < 0.8) continue;   // not burst-like enough to keep
		unit_t u;
		unit_from(r, &u);
		if (try_cluster(R, &u, &out[nout])) {
			known_add(R, &out[nout].frame, r->t_cw);
			nout++;
			continue;
		}
		if (R->np == MAX_POOL) memmove(R->pool, R->pool + 1, sizeof(unit_t) * (--R->np));
		R->pool[R->np++] = u;
	}
	free(inside); free(decoded);
	return nout;
}
