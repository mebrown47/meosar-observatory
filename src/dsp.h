// SPDX-License-Identifier: GPL-3.0-or-later
// Filters and channelisers. Complex samples are float complex.
#ifndef MEOSAR_DSP_H
#define MEOSAR_DSP_H

#include <complex.h>
#include <stddef.h>

typedef float complex cf_t;

// Low-pass FIR, Kaiser-windowed sinc, unit DC gain (as scipy.signal.firwin(ntaps, cutoff,
// fs=fs, window=("kaiser", beta))).
float *fir_lowpass(int ntaps, double cutoff_hz, double fs, double beta);

// Mix f_hz to 0 (phase referenced to sample index n0), low-pass to +-bw/2, decimate by D.
// Output sample k is aligned with input sample k*D (zero delay), like resample_poly.
// Returns the number of output samples written to out (ceil(n / D)).
size_t channelise(const cf_t *x, size_t n, double f_hz, double fs, int D, double bw_hz,
                  long long n0, cf_t *out);

// Streaming version for the live input: push input, pull decimated output with the same
// alignment (output k <-> input k*D), delayed by the filter half-length.
typedef struct {
	double fs, f_hz;
	int D, ntaps;
	float *h;
	cf_t *hist;          // last ntaps-1 + block mixed samples
	size_t hlen, hcap;
	long long n_in;      // input samples consumed so far
	long long k_out;     // next output index
	double ph, dph;      // NCO
} chan_stream_t;

chan_stream_t *chan_stream_new(double fs, double f_hz, int D, double bw_hz);
// Appends n input samples; writes up to max outputs, returns how many.
size_t chan_stream_push(chan_stream_t *c, const cf_t *x, size_t n, cf_t *out, size_t max);
void chan_stream_free(chan_stream_t *c);

#endif
