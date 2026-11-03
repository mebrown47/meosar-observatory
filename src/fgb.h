// SPDX-License-Identifier: GPL-3.0-or-later
// First-generation (FGB) 406 MHz beacon frame (C/S T.001): layout, BCH codes, soft decoding.
//
// Bits are 0-based as transmitted:
//   0..14 bit sync (15 ones)   15..23 frame sync (000101111 normal, 011010000 self-test)
//   24..105 BCH(82,61) t=3: PDF-1 (24 = format flag) + 21 parity
//   106..143 BCH(38,26) t=2: PDF-2 + 12 parity (long frames)
// A codeword is valid when its bits, first bit as the highest power, divide by the generator.
#ifndef FGB_H
#define FGB_H

#include <stdint.h>

#define FGB_BITS 144
#define FGB_BIT_RATE 400.0
#define FGB_CW_SECONDS 0.160
#define FGB_MOD_INDEX 1.1

extern const uint8_t FGB_SYNC_NORMAL[9], FGB_SYNC_SELFTEST[9];
extern const short FGB_MIDS[];
extern const int FGB_NMIDS;

void fgb_init(void);   // builds the syndrome tables; call once before decoding

typedef struct {
	int ok;
	uint8_t bits[FGB_BITS];
	int nbits;            // 112 or 144
	int long_frame, selftest;
	double bch1_cost, bch2_cost;
	int bch1_flips, bch2_flips;
	char reasons[160];    // why not ok, "; "-separated
} fgb_frame_t;

// Soft-decode 144 LLRs (positive = 1) with Chase-II of order L on BCH-1 (and BCH-2 for long
// frames; L capped at 8 there). L = 0 is plain hard-decision table decoding.
void fgb_decode_soft(const double *llr, int L, fgb_frame_t *out);

int fgb_country(const uint8_t *bits);
int fgb_is_orbitography(const uint8_t *bits);
int fgb_mid_valid(int mid);
// Hex of bits 24.. (30 characters for 144 bits, 22 for 112, the last digit zero-padded).
void fgb_hex(const uint8_t *bits, int nbits, char *out);
// Build a valid frame from a 30/15-hex payload (bits 24..), recomputing BCH; returns nbits.
int fgb_build(const char *payload_hex, uint8_t *bits);

#endif
