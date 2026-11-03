// SPDX-License-Identifier: GPL-3.0-or-later
#include "sgb.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static const uint32_t SEED_I[2] = {0x000001, 0x52C9F0}, SEED_Q[2] = {0x1AC1FC, 0x3CE928};
static float *chips[2][2];
static uint8_t gexp[512], glog[256];
static uint8_t synpos[250][12];   // syndrome contribution of a 1 at position i

static void gen(uint32_t s, float *out) {
	for (int i = 0; i < SGB_NCHIPS; i++) {
		out[i] = (s & 1) ? -1.f : 1.f;
		uint32_t fb = (s ^ (s >> 18)) & 1;
		s = (s >> 1) | (fb << 22);
	}
}

void sgb_init(void) {
	if (chips[0][0]) return;
	for (int m = 0; m < 2; m++) {
		chips[m][0] = malloc(sizeof(float) * SGB_NCHIPS);
		chips[m][1] = malloc(sizeof(float) * SGB_NCHIPS);
		gen(SEED_I[m], chips[m][0]);
		gen(SEED_Q[m], chips[m][1]);
	}
	int x = 1;
	for (int i = 0; i < 255; i++) {
		gexp[i] = (uint8_t)x;
		glog[x] = (uint8_t)i;
		x <<= 1;
		if (x & 0x100) x ^= 0x11D;
	}
	for (int i = 255; i < 512; i++) gexp[i] = gexp[i - 255];
	for (int i = 0; i < 250; i++)
		for (int j = 1; j <= 12; j++) synpos[i][j - 1] = gexp[((249 - i) * j) % 255];
}

const float *sgb_chips_i(int mode) { return chips[mode][0]; }
const float *sgb_chips_q(int mode) { return chips[mode][1]; }

static uint8_t gmul(uint8_t a, uint8_t b) { return (a && b) ? gexp[glog[a] + glog[b]] : 0; }
static uint8_t gdiv(uint8_t a, uint8_t b) { return a ? gexp[(glog[a] - glog[b] + 255) % 255] : 0; }

void sgb_parity(const uint8_t *info, uint8_t *par) {
	const uint64_t G = 0x1C7EB85DF3C97ULL;
	uint64_t rem = 0;
	for (int i = 0; i < 202 + 48; i++) {
		rem = rem << 1 | (i < 202 ? info[i] : 0);
		if (rem >> 48 & 1) rem ^= G;
	}
	for (int i = 0; i < 48; i++) par[i] = rem >> (47 - i) & 1;
}

static void syndromes(const uint8_t *cw, uint8_t *s) {
	memset(s, 0, 12);
	for (int i = 0; i < 250; i++)
		if (cw[i])
			for (int j = 0; j < 12; j++) s[j] ^= synpos[i][j];
}

int sgb_syndromes_zero(const uint8_t *cw) {
	uint8_t s[12];
	syndromes(cw, s);
	for (int j = 0; j < 12; j++)
		if (s[j]) return 0;
	return 1;
}

static int bm_from_syn(const uint8_t *s, uint8_t *cw) {
	int any = 0;
	for (int j = 0; j < 12; j++) any |= s[j];
	if (!any) return 0;
	uint8_t C[13] = {1}, B[13] = {1}, T[13];
	int L = 0, m = 1;
	uint8_t b = 1;
	for (int n = 0; n < 12; n++) {
		uint8_t d = s[n];
		for (int i = 1; i <= L; i++) d ^= gmul(C[i], s[n - i]);
		if (d == 0) { m++; continue; }
		uint8_t coef = gdiv(d, b);
		if (2 * L <= n) {
			memcpy(T, C, 13);
			for (int i = m; i < 13; i++) C[i] ^= gmul(coef, B[i - m]);
			L = n + 1 - L;
			memcpy(B, T, 13);
			b = d;
			m = 1;
		} else {
			for (int i = m; i < 13; i++) C[i] ^= gmul(coef, B[i - m]);
			m++;
		}
	}
	if (L > 6) return -1;
	int pos[6], np = 0;
	for (int i = 0; i < 250 && np <= L; i++) {
		int e = 249 - i;
		uint8_t x = gexp[(255 - e) % 255], v = 0, xp = 1;
		for (int k = 0; k <= L; k++) { v ^= gmul(C[k], xp); xp = gmul(xp, x); }
		if (v == 0) { if (np < 6) pos[np] = i; np++; }
	}
	if (np != L) return -1;
	for (int k = 0; k < np; k++) cw[pos[k]] ^= 1;
	return L;
}

int sgb_bm_decode(uint8_t *cw) {
	uint8_t s[12];
	syndromes(cw, s);
	return bm_from_syn(s, cw);
}

int sgb_chase(const double *llr, int L, uint8_t *out, double *cost_out) {
	uint8_t hard[250], s0[12];
	double rel[250];
	for (int i = 0; i < 250; i++) { hard[i] = llr[i] > 0; rel[i] = fabs(llr[i]); }
	syndromes(hard, s0);
	int weak[16];
	uint8_t taken[250] = {0};
	if (L > 16) L = 16;
	for (int j = 0; j < L; j++) {
		int bi = -1;
		for (int i = 0; i < 250; i++)
			if (!taken[i] && (bi < 0 || rel[i] < rel[bi])) bi = i;
		weak[j] = bi;
		taken[bi] = 1;
	}
	double best = INFINITY;
	int found = 0, bestf = -1;
	uint8_t test[250];
	for (uint32_t mask = 0; mask < (1u << L); mask++) {
		uint8_t s[12];
		memcpy(s, s0, 12);
		memcpy(test, hard, 250);
		for (int j = 0; j < L; j++)
			if (mask >> j & 1) {
				test[weak[j]] ^= 1;
				for (int k = 0; k < 12; k++) s[k] ^= synpos[weak[j]][k];
			}
		if (bm_from_syn(s, test) < 0) continue;
		double c = 0;
		int nf = 0;
		for (int i = 0; i < 250; i++)
			if (test[i] != hard[i]) { c += rel[i]; nf++; }
		if (c < best) { best = c; memcpy(out, test, 250); bestf = nf; found = 1; }
	}
	if (!found) return -1;
	*cost_out = best;
	return bestf;
}

void sgb_hex(const uint8_t *bits, char *out) {
	uint8_t b[252] = {0, 0};
	memcpy(b + 2, bits, 250);
	for (int i = 0; i < 63; i++) {
		int v = b[4 * i] << 3 | b[4 * i + 1] << 2 | b[4 * i + 2] << 1 | b[4 * i + 3];
		out[i] = "0123456789ABCDEF"[v];
	}
	out[63] = 0;
}
