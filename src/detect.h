// SPDX-License-Identifier: GPL-3.0-or-later
// FGB burst detector: the 160 ms unmodulated carrier at the start of every FGB burst shows up as
// a narrow line in 131 ms FFTs (16384 points at 125 kS/s, Hann, every 4096 samples). A bin whose
// power exceeds its noise mean (median over the block / ln 2) by thr_db and is a local maximum
// is a candidate; candidates within +-40 Hz and 0.4 s of a stronger one are merged into it.
#ifndef MEOSAR_DETECT_H
#define MEOSAR_DETECT_H

#include "dsp.h"

#define DET_NFFT 16384
#define DET_HOP 4096

typedef struct {
	double t;       // start of the strongest detector frame, s (absolute)
	double f;       // Hz (absolute)
	double snr_db;
} cand_t;

// y: band stream at f1 Hz centred on mid; t_y0: time of y[0]. Returns the number of merged
// candidates (sorted by time) written to out.
int fgb_detect(const cf_t *y, size_t n, double f1, double mid, double lo, double hi,
               double t_y0, double thr_db, cand_t *out, int max);

#endif
