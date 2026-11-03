// SPDX-License-Identifier: GPL-3.0-or-later
// SGB receiver for weak downlink bursts.
//
// Acquisition (PMF-FFT): the 25-bit preamble (6400 known chips on I and Q, 166.7 ms) is split
// into 25 one-bit segments. The band stream is decimated to 62.5 kS/s; for every coarse
// frequency bin (150 Hz apart) each segment replica is correlated with the block at all lags by
// FFT. A cheap non-coherent sum of the 25 segment powers picks the best 32 lags per bin, where a
// 32-point DFT across the segments adds them coherently at fine frequency. Bins run in threads.
// Demodulation (125 kS/s band stream): chip integrals with fractional timing; timing and
// frequency refined on the preamble; per-bit correlations on I and Q; carrier frequency
// residual and drift fitted on the squared bit correlations of the whole burst; residual phase
// smoothed decision-directed; soft bits scaled by the quadrature noise.
#ifndef MEOSAR_SGB_RX_H
#define MEOSAR_SGB_RX_H

#include "dsp.h"
#include "sgb.h"

typedef struct sgb_acq sgb_acq_t;
typedef struct { double t, f, z; } sgb_peak_t;   // t: s from block start, f: Hz from centre

sgb_acq_t *sgb_acq_new(double fs_in, double span_hz, int threads);
// x: band stream block at fs_in. Searches burst starts in [0, t_max) s. Returns peaks sorted by
// z (coherent peak power over the mean non-coherent power; noise ~ max of 32 Exp(1)).
int sgb_acq_search(sgb_acq_t *a, const cf_t *x, size_t n, double t_max, sgb_peak_t *out, int max);
double sgb_acq_span(const sgb_acq_t *a);   // seconds of burst starts covered per call
void sgb_acq_free(sgb_acq_t *a);

typedef struct {
	int ok;
	double llr[SGB_MSG_BITS];
	double t0, f_hz, fdot, ebn0_db, pre_z;
} sgb_demod_t;

void sgb_demodulate(const cf_t *x, size_t n, double fs, double t0_guess, double f_guess, int mode,
                    sgb_demod_t *d);

// One decoded or attempted SGB burst (relays combined).
typedef struct {
	double t, f_hz;
	int nrelay;
	double rf[8], rdt_ms[8], reb[8], rpz[8];
	int ralone[8];
	int ok, hard, nerr, flips, accept;   // accept: 2 single, 1 weak
	char reasons[96];
	uint8_t bits[SGB_MSG_BITS];
	const char *band;
} sgb_rec_t;

// Searches a band-stream block for bursts starting in [t_lo, t_hi) (absolute s; y[0] at t_y0).
// Returns the number of records written to out.
int sgb_process_block(sgb_acq_t *a, const cf_t *y, size_t n, double f1, double mid, double t_y0,
                      double t_lo, double t_hi, double acq_thr, double prez_thr, sgb_rec_t *out, int max);

#endif
