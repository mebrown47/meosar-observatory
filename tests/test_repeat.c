// SPDX-License-Identifier: GPL-3.0-or-later
// Repeat combining on synthetic soft bursts: weak repeats decode together, a weak burst of a
// confirmed beacon is matched to it, siblings are kept apart, and noise yields nothing.
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fgb_repeat.h"

static double gauss(void) {
	double u = (rand() + 1.0) / (RAND_MAX + 2.0), v = (rand() + 1.0) / (RAND_MAX + 2.0);
	return sqrt(-2 * log(u)) * cos(2 * M_PI * v);
}

// A burst's soft bits at a given Eb/N0 (dB), as the demodulator gives them: LLR = 2 r / sigma^2.
static void burst(const uint8_t *bits, double ebn0_db, double t, fgb_rec_t *r) {
	memset(r, 0, sizeof *r);
	double sig = sqrt(1 / (2 * pow(10, ebn0_db / 10.0)));
	for (int i = 0; i < FGB_BITS; i++) r->llr[i] = bits ? 2 * ((bits[i] ? 1 : -1) + sig * gauss()) / (sig * sig) : 2 * sig * gauss() / (sig * sig);
	r->valid = 1; r->sync = ebn0_db > 3 ? 2.5 : 1.0; r->t_cw = t; r->f_hz = 1544.1e6; r->ebn0 = ebn0_db; r->band = "galileo";
	fgb_decode_soft(r->llr, 10, &r->frame);
	r->accept = fgb_accept(r);
}

int main(void) {
	fgb_init();
	srand(11);
	int fails = 0;
	uint8_t a[FGB_BITS], b[FGB_BITS];
	fgb_build("CE013E96B18785040BA88E5A4EC14D", a);   // a long-frame beacon
	fgb_build("8E3301E2402B002BBA863609670908", b);   // another
	fgb_rep_t out[16];

	// 1. Repeats too weak alone (-3 dB) decode together, within a few bursts.
	int decoded = 0, alone = 0, nb = 0;
	for (int k = 0; k < 40; k++) {
		fgb_repeat_t *R = fgb_repeat_new(1200, 10);
		int done = 0;
		for (int i = 0; i < 6 && !done; i++) {
			fgb_rec_t r;
			burst(a, -3.0, 50.0 * i, &r);
			alone += r.frame.ok && !memcmp(r.frame.bits + 24, a + 24, 82);
			nb++;
			r.accept = 0;   // as in the receiver: far below its acceptance threshold
			int n = fgb_repeat_block(R, &r, 1, out, 16);
			for (int j = 0; j < n; j++) {
				if (out[j].kind != 2) continue;
				if (memcmp(out[j].frame.bits + 24, a + 24, 82)) { printf("repeat decode gave the WRONG frame\n"); fails++; }
				else done = 1;
			}
		}
		decoded += done;
		fgb_repeat_free(R);
	}
	printf("weak repeats at -3 dB: decoded together in %d/40 trials within 6 bursts (alone: %d of %d bursts)\n", decoded, alone, nb);
	if (decoded < 30) fails++;

	// 2. A confirmed beacon's weak bursts (-1 dB) are matched to it, and another beacon's never are
	//    (that one may well get decoded from its own repeats, and then matched to itself).
	{
		fgb_repeat_t *R = fgb_repeat_new(1200, 10);
		fgb_rec_t r;
		burst(a, 8.0, 0, &r);                      // strong: decoded, teaches the known table
		if (r.accept != 2) { printf("strong burst not accepted\n"); fails++; }
		fgb_repeat_block(R, &r, 1, out, 16);
		int match_a = 0, match_b = 0;
		for (int i = 1; i <= 50; i++) {
			burst(a, -1.0, 50.0 * i, &r); r.accept = 0;
			int n = fgb_repeat_block(R, &r, 1, out, 16);
			for (int j = 0; j < n; j++) match_a += out[j].kind == 1 && !memcmp(out[j].frame.bits + 24, a + 24, 82);
			burst(b, -1.0, 50.0 * i + 10, &r); r.accept = 0;
			n = fgb_repeat_block(R, &r, 1, out, 16);
			for (int j = 0; j < n; j++) match_b += !memcmp(out[j].frame.bits + 24, a + 24, 82);
		}
		printf("matching: beacon's own weak bursts %d/50 matched to it, another beacon's %d/50\n", match_a, match_b);
		if (match_a < 25 || match_b) fails++;
		fgb_repeat_free(R);
	}

	// 3. A family: two known beacons with consecutive serials, codewords 7 bits apart, both weak
	//    (-2 dB) and interleaved. No burst may be credited to the wrong sibling.
	{
		uint8_t x[FGB_BITS], y[FGB_BITS];
		fgb_build("8C9EF9C0011A2A0D8E7A11CC54D691", x);
		fgb_build("8C9EF9C0031A2A09C2EA154E4CD6D3", y);
		int wrong = 0, right = 0, by_match = 0;
		for (int trial = 0; trial < 20; trial++) {
			fgb_repeat_t *R = fgb_repeat_new(1200, 10);
			fgb_rec_t r;
			burst(x, 8.0, 0, &r);  fgb_repeat_block(R, &r, 1, out, 16);
			burst(y, 8.0, 5, &r);  fgb_repeat_block(R, &r, 1, out, 16);
			for (int i = 1; i <= 24; i++)
				for (int w = 0; w < 2; w++) {
					double t = 50.0 * i + 20 * w;   // x on the 50 s, y 20 s later
					burst(w ? y : x, -2.0, t, &r); r.accept = 0;
					int n = fgb_repeat_block(R, &r, 1, out, 16);
					for (int j = 0; j < n; j++) {
						int is_x = !memcmp(out[j].frame.bits + 24, x + 24, 82);
						for (int m = 0; m < out[j].nmem; m++) {
							int from_x = fmod(out[j].mem[m].t_cw, 50.0) < 10;
							if (is_x == from_x) right++; else wrong++;
							by_match += out[j].kind == 1;
						}
					}
				}
			fgb_repeat_free(R);
		}
		printf("family 7 bits apart at -2 dB: %d bursts credited to the right sibling (%d by matching), %d to the wrong one\n", right, by_match, wrong);
		if (wrong > right / 100) fails++;
	}

	// 4. Noise: no decodes, no matches.
	{
		fgb_repeat_t *R = fgb_repeat_new(1200, 10);
		fgb_rec_t r;
		burst(a, 8.0, 0, &r);
		fgb_repeat_block(R, &r, 1, out, 16);
		int bad = 0;
		for (int i = 0; i < 20000; i++) {
			burst(NULL, 0.0, 1.0 + 0.06 * i, &r);   // 20000 noise "bursts" over 20 min
			r.accept = 0;
			bad += fgb_repeat_block(R, &r, 1, out, 16);
		}
		printf("noise: %d results from 20000 noise bursts\n", bad);
		if (bad) fails++;
		fgb_repeat_free(R);
	}
	printf(fails ? "FAILED\n" : "all ok\n");
	return fails != 0;
}
