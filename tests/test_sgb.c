// SPDX-License-Identifier: GPL-3.0-or-later
// SGB codec checks against T.018 Table 2.2 and published vectors.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sgb.h"

static int hexbits(const char *h, uint8_t *b, int nb) {
	int n = 0;
	for (; *h && n < nb; h++) {
		int v = (*h >= 'A') ? *h - 'A' + 10 : *h - '0';
		for (int j = 3; j >= 0 && n < nb; j--) b[n++] = v >> j & 1;
	}
	return n;
}

int main(void) {
	sgb_init();
	int fails = 0;
	// T.018 Table 2.2: normal I channel, first 64 chips 8000 0108 4212 84A1 (logic 1 -> -1)
	const char *ref = "80000108421284A1";
	uint8_t rb[64];
	hexbits(ref, rb, 64);
	const float *ci = sgb_chips_i(SGB_NORMAL);
	int ok = 1;
	for (int i = 0; i < 64; i++) ok &= (ci[i] < 0) == rb[i];
	printf("PRN Table 2.2: %s\n", ok ? "ok" : "FAIL");
	fails += !ok;
	// sgb-codec BCH vector
	uint8_t info[202], par[48];
	hexbits("00E608F4C986196188A047C000000000000FFFC0100C1A00960", info, 202);
	sgb_parity(info, par);
	const char *ep = "010010010010101001001111110001010111101001001001";
	ok = 1;
	for (int i = 0; i < 48; i++) ok &= par[i] == ep[i] - '0';
	printf("BCH parity vector: %s\n", ok ? "ok" : "FAIL");
	fails += !ok;
	// real frame from dec406 logs (63 hex, two leading pad bits)
	uint8_t f252[252], cw[250];
	hexbits("BFFF67125B8BF83E07FFC1FE00000000001FFFF003FFFFFFE706FFEC0EF52B3", f252, 252);
	memcpy(cw, f252 + 2, 250);
	ok = sgb_syndromes_zero(cw);
	printf("real frame syndrome zero: %s\n", ok ? "ok" : "FAIL");
	fails += !ok;
	// 6 random errors corrected
	uint8_t full[250];
	memcpy(full, info, 202);
	memcpy(full + 202, par, 48);
	srand(3);
	int good = 0;
	for (int t = 0; t < 200; t++) {
		uint8_t e[250];
		memcpy(e, full, 250);
		for (int k = 0; k < 6; k++) e[rand() % 250] ^= 1;   // up to 6 distinct
		int n = sgb_bm_decode(e);
		good += n >= 0 && !memcmp(e, full, 250);
	}
	printf("BM corrects up to 6 errors: %d/200 %s\n", good, good == 200 ? "ok" : "FAIL");
	fails += good != 200;
	char hx[64];
	sgb_hex(cw, hx);
	// the log line had its two pad bits set (B = 1011); T.018 hex pads with zeros (3 = 0011)
	ok = !strcmp(hx, "3FFF67125B8BF83E07FFC1FE00000000001FFFF003FFFFFFE706FFEC0EF52B3");
	printf("hex round trip: %s\n", ok ? "ok" : "FAIL");
	fails += !ok;
	printf("%s\n", fails ? "FAILED" : "all ok");
	return fails != 0;
}
