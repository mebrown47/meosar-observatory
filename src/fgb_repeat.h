// SPDX-License-Identifier: GPL-3.0-or-later
// Repeat combining: a beacon sends the same BCH-1 codeword (PDF-1 + parity, bits 25-106) in every
// burst, about every 50 s. Bursts too weak to decode alone are put to use in two ways (the C port
// of proto/combine.py, in streaming form):
//
//  1. Known beacons. A burst's soft bits are correlated with the BCH-1 codeword of each confirmed
//     beacon (z = sum(L s) / |L|, ~N(0,1) when unrelated); z >= 6, with odds of 1000:1 over every
//     other known beacon (log-likelihood ratio over the bits where the codewords differ),
//     attributes the burst to that beacon, however weak it is. (z alone cannot tell siblings
//     apart: it is capped at sqrt(82), and a perfect burst leads a sibling 7 bits away by 1.55.)
//  2. Unknown beacons. Unmatched bursts are kept for a window (20 min). A new one that correlates
//     with earlier ones (pairwise z >= 5, or >= 4 in a group of three bursts or more) is summed
//     with them and decoded. Every member must match the decoded codeword on its own (z >= 3) and
//     favour it over every other known beacon's (log-likelihood ratio >= 0); members that favour
//     another are left out and the rest decoded again. The result stands if the sum favours it
//     over every other known beacon by odds of 1000:1, BCH-1 needs at most 3 flips, and PDF-2,
//     where there is one, decodes. (Beacons of one family, consecutive serials, have codewords
//     only 7-12 bits apart: their weak bursts link and sum together.)
//
// A burst is a relay-combined record, or a relay that was not combined with others.
#ifndef MEOSAR_FGB_REPEAT_H
#define MEOSAR_FGB_REPEAT_H

#include "fgb_demod.h"

#define FGB_REP_MAXMEM 16

typedef struct {
	double t_cw, f_hz, ebn0, z;
	double lodds;             // log-likelihood ratio, this codeword against the likeliest other known one
	const char *band;
} fgb_rep_member_t;

typedef struct {
	int kind;                 // 1: burst matched to a known beacon; 2: new frame decoded from repeats
	fgb_frame_t frame;        // the beacon's frame (kind 1: the known beacon's last frame)
	double z, z2;             // kind 1: correlation with the beacon's codeword, and with the runner-up's
	double lodds;             // log-likelihood ratio against the likeliest other known codeword
	                          // (kind 1: of the burst; kind 2: of the summed members); capped at 999
	int nmem, bursts;         // members, and distinct transmissions among them
	fgb_rep_member_t mem[FGB_REP_MAXMEM];
} fgb_rep_t;

typedef struct fgb_repeat fgb_repeat_t;

fgb_repeat_t *fgb_repeat_new(double window_s, int L);
void fgb_repeat_free(fgb_repeat_t *R);
// One band's records for one block (after relay combining). Writes up to cap results; returns
// how many.
int fgb_repeat_block(fgb_repeat_t *R, const fgb_rec_t *recs, int n, fgb_rep_t *out, int cap);

#endif
