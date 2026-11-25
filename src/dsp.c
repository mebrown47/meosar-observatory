// SPDX-License-Identifier: GPL-3.0-or-later
#include "dsp.h"

#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

static pthread_mutex_t fftw_plan_mu = PTHREAD_MUTEX_INITIALIZER;
void fftw_plan_lock(void) { pthread_mutex_lock(&fftw_plan_mu); }
void fftw_plan_unlock(void) { pthread_mutex_unlock(&fftw_plan_mu); }

static double bessel_i0(double x) {
	double s = 1, t = 1;
	for (int k = 1; k < 60; k++) {
		t *= (x / (2 * k)) * (x / (2 * k));
		s += t;
		if (t < 1e-12 * s) break;
	}
	return s;
}

float *fir_lowpass(int ntaps, double cutoff_hz, double fs, double beta) {
	float *h = malloc(sizeof(float) * (size_t)ntaps);
	double fc = cutoff_hz / (fs / 2), m = (ntaps - 1) / 2.0, sum = 0;
	double *w = malloc(sizeof(double) * (size_t)ntaps);
	for (int i = 0; i < ntaps; i++) {
		double x = i - m;
		double sinc = x == 0 ? 1.0 : sin(M_PI * fc * x) / (M_PI * fc * x);
		double r = 2.0 * i / (ntaps - 1) - 1.0;
		double kw = bessel_i0(beta * sqrt(fmax(0.0, 1 - r * r))) / bessel_i0(beta);
		w[i] = fc * sinc * kw;
		sum += w[i];
	}
	for (int i = 0; i < ntaps; i++) h[i] = (float)(w[i] / sum);
	free(w);
	return h;
}

size_t channelise(const cf_t *x, size_t n, double f_hz, double fs, int D, double bw_hz,
                  long long n0, cf_t *out) {
	int ntaps = 8 * D + 1, c = (ntaps - 1) / 2;
	float *h = fir_lowpass(ntaps, bw_hz / 2, fs, 8.0);
	cf_t *m = malloc(sizeof(cf_t) * n);
	double w = -2 * M_PI * f_hz / fs;
	// mixer: exact phase per sample (double), referenced to n0
	for (size_t i = 0; i < n; i++) {
		double ph = fmod(w * (double)((long long)i + n0), 2 * M_PI);
		m[i] = x[i] * (float)cos(ph) + x[i] * (float)sin(ph) * I;
	}
	size_t nout = (n + (size_t)D - 1) / (size_t)D;
	for (size_t k = 0; k < nout; k++) {
		long long base = (long long)k * D + c;
		float re = 0, im = 0;
		int j0 = 0, j1 = ntaps - 1;
		if (base - j1 < 0) j1 = (int)base;
		if (base - j0 >= (long long)n) j0 = (int)(base - (long long)n + 1);
		for (int j = j0; j <= j1; j++) {
			cf_t v = m[base - j];
			re += h[j] * crealf(v);
			im += h[j] * cimagf(v);
		}
		out[k] = re + im * I;
	}
	free(m);
	free(h);
	return nout;
}

chan_stream_t *chan_stream_new(double fs, double f_hz, int D, double bw_hz) {
	chan_stream_t *c = calloc(1, sizeof *c);
	c->fs = fs;
	c->f_hz = f_hz;
	c->D = D;
	c->ntaps = 8 * D + 1;
	c->h = fir_lowpass(c->ntaps, bw_hz / 2, fs, 8.0);
	c->dph = -2 * M_PI * f_hz / fs;
	return c;
}

// hist holds mixed input samples with absolute indices [n_in - hlen, n_in).
size_t chan_stream_push(chan_stream_t *c, const cf_t *x, size_t n, cf_t *out, size_t max) {
	if (c->hlen + n > c->hcap) {
		c->hcap = (c->hlen + n) * 2 + 1024;
		c->hist = realloc(c->hist, sizeof(cf_t) * c->hcap);
	}
	// rotating phasor in double, renormalised every 1024 samples (exact phase recomputed)
	double complex rot = cexp(I * c->ph), step = cexp(I * c->dph);
	for (size_t i = 0; i < n; i++) {
		c->hist[c->hlen + i] = x[i] * (cf_t)rot;
		rot *= step;
		if (((c->n_in + (long long)i + 1) & 1023) == 0) {
			c->ph = fmod(c->dph * (double)(c->n_in + (long long)i + 1), 2 * M_PI);
			rot = cexp(I * c->ph);
		}
	}
	c->ph = carg(rot);
	c->hlen += n;
	c->n_in += (long long)n;
	int half = (c->ntaps - 1) / 2;
	long long first = c->n_in - (long long)c->hlen;   // absolute index of hist[0]
	size_t w = 0;
	while (w < max) {
		long long centre = c->k_out * c->D;           // output k aligned with input k*D
		if (centre + half >= c->n_in) break;          // needs samples up to centre+half
		float re = 0, im = 0;
		for (int j = 0; j < c->ntaps; j++) {
			long long idx = centre + half - j;
			if (idx < first) break;                   // before the stream start: zeros
			cf_t v = c->hist[idx - first];
			re += c->h[j] * crealf(v);
			im += c->h[j] * cimagf(v);
		}
		out[w++] = re + im * I;
		c->k_out++;
	}
	// keep only what future outputs need
	long long need_from = c->k_out * c->D - half;
	if (need_from > first) {
		size_t drop = (size_t)(need_from - first);
		if (drop > c->hlen) drop = c->hlen;
		memmove(c->hist, c->hist + drop, sizeof(cf_t) * (c->hlen - drop));
		c->hlen -= drop;
	}
	return w;
}

void chan_stream_free(chan_stream_t *c) {
	if (!c) return;
	free(c->h);
	free(c->hist);
	free(c);
}
