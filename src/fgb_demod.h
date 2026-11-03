// SPDX-License-Identifier: GPL-3.0-or-later
// Coherent FGB demodulator for weak bursts, working on a whole burst held in memory.
//
// After channelising to ~20 kS/s around the burst carrier:
//   x(t) = A cos1.1 e^{j theta} + j A sin1.1 m(t) e^{j theta} + n(t)
// The residual carrier never goes away, so it is an absolute phase reference and the
// Manchester data m(t) sit in quadrature. Steps: carrier frequency and drift over the burst
// (dechirp grid + zero-padded FFT, then a local search); timing from the 24 known bits;
// decision-directed phase smoothed forwards and backwards; half-bit-difference soft bits,
// LLR = 2 mu r / sigma^2 with sigma from the noise before the burst.
#ifndef MEOSAR_FGB_DEMOD_H
#define MEOSAR_FGB_DEMOD_H

#include "dsp.h"
#include "fgb.h"

typedef struct {
	int ok;
	double llr[FGB_BITS];
	double f_hz, fdot, t_data, t_cw;
	int polarity, selftest_sync;
	double sync_score, cn0_dbhz, ebn0_db;
} fgb_demod_t;

void fgb_demodulate(const cf_t *x, size_t n, double fs, double t_cw_guess, fgb_demod_t *d);

// One candidate from a band stream (rate f1, centre mid, y[0] at t_y0).
typedef struct {
	int valid;              // demodulated
	double t_frame, f_det, det_snr;
	double t_cw, f_hz, fdot, ebn0, cn0, sync;
	int selftest_sync, pol;
	fgb_frame_t frame;
	double llr[FGB_BITS];
	const char *band;
	int accept;             // 0 none, 1 weak, 2 single
	int combined;           // >0: number of relays summed into this record
	int relay_idx[8];       // for combined records: indices into the per-block array
} fgb_rec_t;

int fgb_process_candidate(const cf_t *y, size_t n, double f1, double mid, double t_y0,
                          double t_frame, double f_abs, int L, fgb_rec_t *r);
int fgb_accept(const fgb_rec_t *r);
// Groups relays of one burst (within 40 ms, carriers > 20 Hz apart) and decodes their summed
// LLRs. Appends combined records after the n inputs; returns the new total.
int fgb_combine_relays(fgb_rec_t *recs, int n, int cap, int L);

#endif
