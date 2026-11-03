// SPDX-License-Identifier: GPL-3.0-or-later
// Codec checks: dec406's FGB test vector, Chase gain on noisy frames.
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fgb.h"

static double gauss(void) {
	double u = (rand() + 1.0) / (RAND_MAX + 2.0), v = (rand() + 1.0) / (RAND_MAX + 2.0);
	return sqrt(-2 * log(u)) * cos(2 * M_PI * v);
}

int main(void) {
	fgb_init();
	int fails = 0;
	const char *h = "8E3301E2402B002BBA863609670908";
	uint8_t bits[FGB_BITS];
	int n = fgb_build(h, bits);
	char out[40];
	fgb_hex(bits, n, out);
	printf("rebuild %s -> %s %s\n", h, out, strcmp(h, out) ? "FAIL" : "ok");
	fails += strcmp(h, out) != 0;
	srand(7);
	for (int eb = 1; eb <= 4; eb++) {
		double sig = sqrt(1 / (2 * pow(10, eb / 10.0)));
		int ok_h = 0, ok_c = 0, N = 300;
		for (int t = 0; t < N; t++) {
			double llr[FGB_BITS];
			for (int i = 0; i < FGB_BITS; i++) llr[i] = 2 * ((bits[i] ? 1 : -1) + sig * gauss()) / (sig * sig);
			fgb_frame_t f;
			fgb_decode_soft(llr, 0, &f);
			ok_h += f.nbits && !memcmp(f.bits + 24, bits + 24, 82);
			fgb_decode_soft(llr, 10, &f);
			ok_c += f.nbits && !memcmp(f.bits + 24, bits + 24, 82);
		}
		printf("Eb/N0 %d dB: BCH-1 hard %5.1f%%  Chase L=10 %5.1f%%\n", eb, 100.0 * ok_h / N, 100.0 * ok_c / N);
		if (eb == 3 && ok_c < N * 0.97) fails++;
	}
	printf("%s\n", fails ? "FAILED" : "all ok");
	return fails != 0;
}
