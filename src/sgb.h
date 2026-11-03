// SPDX-License-Identifier: GPL-3.0-or-later
// Second-generation (SGB) 406 MHz beacon (C/S T.018): spreading codes, BCH(250,202), soft decoding.
//
// OQPSK DSSS, 38.4 kchip/s on each of I and Q, Q delayed by half a chip. Each channel: 150 bits
// of 256 chips (1.0 s), 25 preamble bits (0) then 125 message bits; message bits interleaved
// I0,Q0,I1,Q1,... Chips are the PRN, inverted for a data 1; PRN logic 0 -> +1. PRN: Fibonacci LFSR
// x^23 + x^18 + 1, output X0, feedback X0^X18 into X22, one continuous sequence per channel.
// BCH(250,202): shortened BCH(255,207), t = 6, GF(2^8) with primitive 0x11D; bit 0 = X^249.
#ifndef MEOSAR_SGB_H
#define MEOSAR_SGB_H

#include <stdint.h>

#define SGB_CHIP_RATE 38400.0
#define SGB_CHIPS_PER_BIT 256
#define SGB_BITS_PER_CH 150
#define SGB_PRE_BITS 25
#define SGB_MSG_BITS 250
#define SGB_NCHIPS (SGB_BITS_PER_CH * SGB_CHIPS_PER_BIT)

enum { SGB_NORMAL = 0, SGB_SELFTEST = 1 };

void sgb_init(void);
// +-1 chips for the whole burst (SGB_NCHIPS each), mode SGB_NORMAL or SGB_SELFTEST.
const float *sgb_chips_i(int mode);
const float *sgb_chips_q(int mode);

void sgb_parity(const uint8_t *info202, uint8_t *par48);
int sgb_syndromes_zero(const uint8_t *cw250);
// Berlekamp-Massey + Chien, in place; returns errors corrected (0..6) or -1.
int sgb_bm_decode(uint8_t *cw250);
// Chase-II over the algebraic decoder; llr > 0 means 1. Returns flips vs hard decisions or -1.
int sgb_chase(const double *llr, int L, uint8_t *out250, double *cost);
// T.018 hex: 63 characters, two zero bits padded on the left.
void sgb_hex(const uint8_t *bits250, char *out64);

#endif
