// SPDX-License-Identifier: GPL-3.0-or-later
#include "fgb.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const uint8_t FGB_SYNC_NORMAL[9] = {0, 0, 0, 1, 0, 1, 1, 1, 1};
const uint8_t FGB_SYNC_SELFTEST[9] = {0, 1, 1, 0, 1, 0, 0, 0, 0};

typedef struct {
	int n, r;              // code length, parity bits
	uint32_t g;            // generator
	uint32_t pos_rem[82];  // remainder of a 1 at position i (i = 0 highest power)
	uint32_t *table;       // syndrome -> packed error positions, 0xFFFFFFFF = none
} bch_t;

static bch_t B1 = {82, 21, 0x26D9E3, {0}, NULL};   // 1001101101100111100011
static bch_t B2 = {38, 12, 0x1539, {0}, NULL};     // 1010100111001

#define PACK(n, a, b, c) ((uint32_t)(n) | (uint32_t)(a) << 2 | (uint32_t)(b) << 9 | (uint32_t)(c) << 16)

static void bch_build(bch_t *c) {
	uint32_t r = 1, top = 1u << c->r;
	uint32_t rem_pow[82];
	for (int k = 0; k < c->n; k++) {   // x^k mod g
		rem_pow[k] = r;
		r <<= 1;
		if (r & top) r ^= c->g;
	}
	for (int i = 0; i < c->n; i++) c->pos_rem[i] = rem_pow[c->n - 1 - i];
	size_t sz = (size_t)1 << c->r;
	c->table = malloc(sz * sizeof(uint32_t));
	memset(c->table, 0xFF, sz * sizeof(uint32_t));
	c->table[0] = PACK(0, 0, 0, 0);
	int t = c->r == 21 ? 3 : 2;
	for (int a = 0; a < c->n; a++) {
		uint32_t s = c->pos_rem[a];
		if (c->table[s] == 0xFFFFFFFFu) c->table[s] = PACK(1, a, 0, 0);
	}
	for (int a = 0; a < c->n; a++)
		for (int b = a + 1; b < c->n; b++) {
			uint32_t s = c->pos_rem[a] ^ c->pos_rem[b];
			if (c->table[s] == 0xFFFFFFFFu) c->table[s] = PACK(2, a, b, 0);
		}
	if (t >= 3)
		for (int a = 0; a < c->n; a++)
			for (int b = a + 1; b < c->n; b++)
				for (int d = b + 1; d < c->n; d++) {
					uint32_t s = c->pos_rem[a] ^ c->pos_rem[b] ^ c->pos_rem[d];
					if (c->table[s] == 0xFFFFFFFFu) c->table[s] = PACK(3, a, b, d);
				}
}

void fgb_init(void) {
	if (B1.table) return;
	bch_build(&B1);
	bch_build(&B2);
}

// Chase-II: returns 1 and the codeword (in out) closest to the soft input, cost and flips.
static int bch_chase(const bch_t *c, const double *llr, int L, uint8_t *out, double *cost_out, int *flips_out) {
	int n = c->n;
	uint8_t hard[82];
	double rel[82];
	uint32_t s0 = 0;
	for (int i = 0; i < n; i++) {
		hard[i] = llr[i] > 0;
		rel[i] = fabs(llr[i]);
		if (hard[i]) s0 ^= c->pos_rem[i];
	}
	if (L > n) L = n;
	int weak[16];
	// L least reliable positions (simple selection)
	uint8_t taken[82] = {0};
	for (int j = 0; j < L; j++) {
		int best = -1;
		for (int i = 0; i < n; i++)
			if (!taken[i] && (best < 0 || rel[i] < rel[best])) best = i;
		weak[j] = best;
		taken[best] = 1;
	}
	double best_cost = INFINITY;
	uint8_t best_flip[82];
	int found = 0;
	for (uint32_t mask = 0; mask < (1u << L); mask++) {
		uint32_t s = s0;
		for (int j = 0; j < L; j++)
			if (mask >> j & 1) s ^= c->pos_rem[weak[j]];
		uint32_t e = c->table[s];
		if (e == 0xFFFFFFFFu) continue;
		uint8_t flip[82] = {0};
		for (int j = 0; j < L; j++)
			if (mask >> j & 1) flip[weak[j]] ^= 1;
		int ne = e & 3;
		if (ne >= 1) flip[(e >> 2) & 0x7F] ^= 1;
		if (ne >= 2) flip[(e >> 9) & 0x7F] ^= 1;
		if (ne >= 3) flip[(e >> 16) & 0x7F] ^= 1;
		double cost = 0;
		for (int i = 0; i < n; i++)
			if (flip[i]) cost += rel[i];
		if (cost < best_cost) {
			best_cost = cost;
			memcpy(best_flip, flip, (size_t)n);
			found = 1;
		}
	}
	if (!found) return 0;
	int nf = 0;
	for (int i = 0; i < n; i++) {
		out[i] = hard[i] ^ best_flip[i];
		nf += best_flip[i];
	}
	*cost_out = best_cost;
	*flips_out = nf;
	return 1;
}

int fgb_country(const uint8_t *bits) {
	int v = 0;
	for (int i = 26; i < 36; i++) v = v << 1 | bits[i];
	return v;
}

int fgb_is_orbitography(const uint8_t *bits) {
	return bits[25] && !bits[36] && !bits[37] && !bits[38];
}

int fgb_mid_valid(int mid) {
	for (int i = 0; i < FGB_NMIDS; i++)
		if (FGB_MIDS[i] == mid) return 1;
	return 0;
}

static void add_reason(fgb_frame_t *f, const char *r) {
	size_t l = strlen(f->reasons);
	snprintf(f->reasons + l, sizeof f->reasons - l, "%s%s", l ? "; " : "", r);
}

void fgb_decode_soft(const double *llr, int L, fgb_frame_t *f) {
	memset(f, 0, sizeof *f);
	f->bch1_cost = INFINITY;
	f->bch1_flips = -1;
	uint8_t hard[FGB_BITS];
	for (int i = 0; i < FGB_BITS; i++) hard[i] = llr[i] > 0;
	int dn = 0, ds = 0;
	for (int i = 0; i < 9; i++) {
		dn += hard[15 + i] != FGB_SYNC_NORMAL[i];
		ds += hard[15 + i] != FGB_SYNC_SELFTEST[i];
	}
	f->selftest = ds < dn;
	char tmp[64];
	int dmin = dn < ds ? dn : ds;
	uint8_t c1[82];
	int have1 = bch_chase(&B1, llr + 24, L, c1, &f->bch1_cost, &f->bch1_flips);
	// Frame sync: at most 1 of 9 bits wrong, or 2 when BCH-1 found a strong codeword (<= 3 flips,
	// within its guaranteed correction), which noise rarely produces.
	if (dmin > 2 || (dmin == 2 && !(have1 && f->bch1_flips <= 3))) {
		snprintf(tmp, sizeof tmp, "sync errors %d", dmin);
		add_reason(f, tmp);
	}
	if (!have1) {
		add_reason(f, "bch1 no codeword");
		return;
	}
	memcpy(f->bits, hard, 24);
	memcpy(f->bits + 24, c1, 82);
	f->long_frame = c1[0];
	if (f->long_frame) {
		f->nbits = 144;
		if (fgb_is_orbitography(f->bits)) {
			memcpy(f->bits + 106, hard + 106, 38);   // not BCH-2 protected
		} else {
			// BCH-1 alone carries the beacon ID; a doubtful PDF-2 is kept and flagged, not fatal.
			uint8_t c2[38];
			if (bch_chase(&B2, llr + 106, L < 8 ? L : 8, c2, &f->bch2_cost, &f->bch2_flips)) {
				memcpy(f->bits + 106, c2, 38);
				f->pdf2 = f->bch2_flips <= 2 ? 1 : 2;
			} else {
				memcpy(f->bits + 106, hard + 106, 38);
				f->bch2_flips = -1;
				f->bch2_cost = 0;
				f->pdf2 = 2;
			}
		}
	} else {
		f->nbits = 112;
		memcpy(f->bits + 106, hard + 106, 6);
	}
	int mid = fgb_country(f->bits);
	if (!fgb_mid_valid(mid)) {
		snprintf(tmp, sizeof tmp, "country %d not assigned", mid);
		add_reason(f, tmp);
	}
	f->ok = f->reasons[0] == 0;
}

void fgb_hex(const uint8_t *bits, int nbits, char *out) {
	int n = nbits - 24, k = 0;
	for (int i = 0; i < n; i += 4) {
		int v = 0;
		for (int j = 0; j < 4; j++) v = v << 1 | (i + j < n ? bits[24 + i + j] : 0);
		out[k++] = "0123456789ABCDEF"[v];
	}
	out[k] = 0;
}

static void bch_encode(const bch_t *c, const uint8_t *msg, uint8_t *par) {
	int k = c->n - c->r;
	uint32_t rem = 0, top = 1u << c->r;
	for (int i = 0; i < k + c->r; i++) {
		rem = rem << 1 | (i < k ? msg[i] : 0);
		if (rem & top) rem ^= c->g;
	}
	for (int i = 0; i < c->r; i++) par[i] = rem >> (c->r - 1 - i) & 1;
}

int fgb_build(const char *hex, uint8_t *bits) {
	uint8_t p[120] = {0};
	int np = 0;
	for (const char *h = hex; *h && np < 120; h++) {
		int v = (*h >= 'a') ? *h - 'a' + 10 : (*h >= 'A') ? *h - 'A' + 10 : *h - '0';
		for (int j = 3; j >= 0; j--) p[np++] = v >> j & 1;
	}
	memset(bits, 1, 15);
	memcpy(bits + 15, FGB_SYNC_NORMAL, 9);
	memcpy(bits + 24, p, 61);
	bch_encode(&B1, bits + 24, bits + 85);
	if (p[0]) {
		memcpy(bits + 106, p + 82, 26);
		bch_encode(&B2, bits + 106, bits + 132);
		return 144;
	}
	memcpy(bits + 106, p + 82, 6);
	return 112;
}
