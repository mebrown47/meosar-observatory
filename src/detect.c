// SPDX-License-Identifier: GPL-3.0-or-later
#include "detect.h"

#include <fftw3.h>
#define RE(a, i) (((float *)(a))[2 * (i)])
#define IM(a, i) (((float *)(a))[2 * (i) + 1])
#include <math.h>
#include <stdlib.h>
#include <string.h>

static float kth(float *a, int n, int k) {   // quickselect, destroys a
	int lo = 0, hi = n - 1;
	while (lo < hi) {
		float p = a[(lo + hi) / 2];
		int i = lo, j = hi;
		while (i <= j) {
			while (a[i] < p) i++;
			while (a[j] > p) j--;
			if (i <= j) { float t = a[i]; a[i] = a[j]; a[j] = t; i++; j--; }
		}
		if (k <= j) hi = j;
		else if (k >= i) lo = i;
		else return a[k];
	}
	return a[k];
}

static int by_snr_desc(const void *a, const void *b) {
	double d = ((const cand_t *)b)->snr_db - ((const cand_t *)a)->snr_db;
	return (d > 0) - (d < 0);
}
static int by_time(const void *a, const void *b) {
	double d = ((const cand_t *)a)->t - ((const cand_t *)b)->t;
	return (d > 0) - (d < 0);
}

int fgb_detect(const cf_t *y, size_t n, double f1, double mid, double lo, double hi,
               double t_y0, double thr_db, cand_t *out, int max) {
	const int N = DET_NFFT;
	if (n < (size_t)N) return 0;
	int nfr = (int)((n - (size_t)N) / DET_HOP) + 1;
	if (nfr < 8) return 0;
	// bins inside [lo, hi): fftshifted index s -> frequency (s - N/2) * f1/N + mid
	int s0 = -1, s1 = -1;
	for (int s = 0; s < N; s++) {
		double f = (s - N / 2) * f1 / N + mid;
		if (f >= lo && f < hi) { if (s0 < 0) s0 = s; s1 = s; }
	}
	int nb = s1 - s0 + 1;
	// One plan per thread: the bands are searched in parallel.
	static _Thread_local fftwf_plan plan;
	static _Thread_local fftwf_complex *in, *outf;
	static _Thread_local float win[DET_NFFT];
	if (!plan) {
		in = fftwf_malloc(sizeof(fftwf_complex) * N);
		outf = fftwf_malloc(sizeof(fftwf_complex) * N);
		fftw_plan_lock();
		plan = fftwf_plan_dft_1d(N, in, outf, FFTW_FORWARD, FFTW_MEASURE);
		fftw_plan_unlock();
		for (int i = 0; i < N; i++) win[i] = (float)(0.5 - 0.5 * cos(2 * M_PI * i / (N - 1)));
	}
	float *P = malloc(sizeof(float) * (size_t)nfr * nb);
	for (int fr = 0; fr < nfr; fr++) {
		const cf_t *x = y + (size_t)fr * DET_HOP;
		for (int i = 0; i < N; i++) {
			RE(in, i) = crealf(x[i]) * win[i];
			IM(in, i) = cimagf(x[i]) * win[i];
		}
		fftwf_execute(plan);
		for (int b = 0; b < nb; b++) {
			int s = s0 + b, k = (s - N / 2 + N) % N;   // unshift
			P[(size_t)fr * nb + b] = RE(outf, k) * RE(outf, k) + IM(outf, k) * IM(outf, k);
		}
	}
	// per-bin noise mean from the median over frames
	float *noise = malloc(sizeof(float) * nb), *col = malloc(sizeof(float) * nfr);
	for (int b = 0; b < nb; b++) {
		for (int fr = 0; fr < nfr; fr++) col[fr] = P[(size_t)fr * nb + b];
		float med = kth(col, nfr, nfr / 2);
		noise[b] = med / (float)M_LN2;
	}
	float thr = (float)pow(10, thr_db / 10);
	int nc = 0, cap = 4096;
	cand_t *c = malloc(sizeof(cand_t) * cap);
	for (int fr = 0; fr < nfr; fr++)
		for (int b = 0; b < nb; b++) {
			float r = P[(size_t)fr * nb + b] / noise[b];
			if (!(r > thr)) continue;
			int lm = 1;
			for (int i = fr - 1; i <= fr + 1 && lm; i++)
				for (int j = b - 2; j <= b + 2; j++) {
					if (i < 0 || i >= nfr || j < 0 || j >= nb) continue;
					if (P[(size_t)i * nb + j] / noise[j] > r) { lm = 0; break; }
				}
			if (!lm) continue;
			if (nc == cap) { cap *= 2; c = realloc(c, sizeof(cand_t) * cap); }
			c[nc].t = t_y0 + (double)fr * DET_HOP / f1;
			c[nc].f = (s0 + b - N / 2) * f1 / N + mid;
			c[nc].snr_db = 10 * log10(r);
			nc++;
		}
	free(P); free(noise); free(col);
	// merge: strongest first, drop those within 40 Hz / 0.4 s of a kept one
	qsort(c, (size_t)nc, sizeof(cand_t), by_snr_desc);
	int nk = 0;
	for (int i = 0; i < nc; i++) {
		int dup = 0;
		for (int k = 0; k < nk && !dup; k++)
			if (fabs(c[i].f - out[k].f) < 40 && fabs(c[i].t - out[k].t) < 0.4) dup = 1;
		if (!dup && nk < max) out[nk++] = c[i];
	}
	free(c);
	qsort(out, (size_t)nk, sizeof(cand_t), by_time);
	return nk;
}
