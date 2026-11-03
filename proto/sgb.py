"""Second-generation (SGB) 406 MHz beacon: spreading codes, frame, BCH(250,202), soft decoding.

Signal (C/S T.018): OQPSK DSSS, 38.4 kchip/s on each of I and Q, Q delayed by half a chip.
Each channel carries 150 bits of 256 chips (1.0 s): 25 preamble bits (all 0) then 125 message
bits; the 250 message bits are interleaved I0,Q0,I1,Q1,... Chips are the PRN, inverted for a
data 1; PRN logic 0 -> +1. PRN: Fibonacci LFSR x^23 + x^18 + 1, output X0, feedback X0^X18 into
X22, one continuous 38 400-chip sequence per channel per burst.
BCH(250,202): shortened BCH(255,207), t = 6, generator below (T.018 Appendix B); bit 0 is the
coefficient of X^249.
"""
from __future__ import annotations

import itertools

import numpy as np

CHIP_RATE = 38400.0
CHIPS_PER_BIT = 256
BITS_PER_CH = 150
PRE_BITS = 25
MSG_BITS = 250
SEEDS = {"normal": (0x000001, 0x1AC1FC), "selftest": (0x52C9F0, 0x3CE928)}
BCH_GEN = 0x1C7EB85DF3C97  # degree 48
NPAR = 48
NINFO = 202


def prn(seed: int, n: int = BITS_PER_CH * CHIPS_PER_BIT) -> np.ndarray:
    """Logic chips 0/1."""
    out = np.empty(n, dtype=np.uint8)
    s = seed
    for i in range(n):
        out[i] = s & 1
        fb = (s ^ (s >> 18)) & 1
        s = (s >> 1) | (fb << 22)
    return out


_prn_cache: dict = {}


def chips_pm(mode: str = "normal"):
    """(+-1 I chips, +-1 Q chips) for the whole burst, float32."""
    if mode not in _prn_cache:
        si, sq = SEEDS[mode]
        _prn_cache[mode] = tuple((1.0 - 2.0 * prn(s)).astype(np.float32) for s in (si, sq))
    return _prn_cache[mode]


def bch_parity(info) -> list[int]:
    rem = 0
    for b in info:
        rem = (rem << 1) | int(b)
        if rem >> NPAR & 1:
            rem ^= BCH_GEN
    for _ in range(NPAR):
        rem <<= 1
        if rem >> NPAR & 1:
            rem ^= BCH_GEN
    return [(rem >> (NPAR - 1 - i)) & 1 for i in range(NPAR)]


def encode(info) -> list[int]:
    info = list(map(int, info))
    assert len(info) == NINFO
    return info + bch_parity(info)


# ---- GF(2^8), primitive 0x11D, for the algebraic decoder ----
GF_EXP = [0] * 512
GF_LOG = [0] * 256
_x = 1
for _i in range(255):
    GF_EXP[_i] = _x
    GF_LOG[_x] = _i
    _x <<= 1
    if _x & 0x100:
        _x ^= 0x11D
for _i in range(255, 512):
    GF_EXP[_i] = GF_EXP[_i - 255]


def gmul(a, b):
    return 0 if a == 0 or b == 0 else GF_EXP[GF_LOG[a] + GF_LOG[b]]


def gdiv(a, b):
    return 0 if a == 0 else GF_EXP[(GF_LOG[a] - GF_LOG[b]) % 255]


# syndrome contribution of a 1 at position i (coefficient of X^(249-i)) for S1..S12
_SYN_POS = np.array([[GF_EXP[((249 - i) * j) % 255] for j in range(1, 13)] for i in range(250)],
                    dtype=np.int64)


def syndromes(cw) -> list[int]:
    idx = np.nonzero(np.asarray(cw, dtype=np.uint8))[0]
    if len(idx) == 0:
        return [0] * 12
    return [int(v) for v in np.bitwise_xor.reduce(_SYN_POS[idx], axis=0)]


def bm_decode(cw, syn=None):
    """Berlekamp-Massey + Chien. Returns (corrected list, n_errors) or (None, -1)."""
    s = syndromes(cw) if syn is None else syn
    if not any(s):
        return list(map(int, cw)), 0
    C, B = [1] + [0] * 12, [1] + [0] * 12
    L, m, b = 0, 1, 1
    for n in range(12):
        d = s[n]
        for i in range(1, L + 1):
            d ^= gmul(C[i], s[n - i])
        if d == 0:
            m += 1
        elif 2 * L <= n:
            T = C[:]
            coef = gdiv(d, b)
            for i in range(m, 13):
                C[i] ^= gmul(coef, B[i - m])
            L, B, b, m = n + 1 - L, T, d, 1
        else:
            coef = gdiv(d, b)
            for i in range(m, 13):
                C[i] ^= gmul(coef, B[i - m])
            m += 1
    if L > 6:
        return None, -1
    # Chien search over the 250 used positions: position i <-> X^(249-i), root at alpha^-(249-i)
    errs = []
    for i in range(250):
        e = 249 - i
        x = GF_EXP[(255 - e) % 255]
        v, xp = 0, 1
        for k in range(L + 1):
            v ^= gmul(C[k], xp)
            xp = gmul(xp, x)
        if v == 0:
            errs.append(i)
    if len(errs) != L:
        return None, -1
    out = list(map(int, cw))
    for i in errs:
        out[i] ^= 1
    return out, L


def chase(llr: np.ndarray, L: int = 8):
    """Chase-II over the algebraic decoder. llr > 0 means bit 1. Returns
    (bits, cost, n_flipped_vs_hard) or (None, inf, -1)."""
    hard = (np.asarray(llr) > 0).astype(np.uint8)
    rel = np.abs(llr)
    s0 = np.array(syndromes(hard), dtype=np.int64)
    weak = np.argsort(rel)[:L]
    best, best_cost = None, np.inf
    for mask in range(1 << L):
        flips = [int(weak[j]) for j in range(L) if mask >> j & 1]
        s = s0.copy()
        for p in flips:
            s ^= _SYN_POS[p]
        test = hard.copy()
        test[flips] ^= 1
        out, n = bm_decode(test, [int(v) for v in s])
        if out is None:
            continue
        diff = np.nonzero(np.asarray(out, dtype=np.uint8) ^ hard)[0]
        cost = float(rel[diff].sum())
        if cost < best_cost:
            best, best_cost = out, cost
    if best is None:
        return None, np.inf, -1
    return best, best_cost, int(np.sum(np.asarray(best, dtype=np.uint8) ^ hard))


def to_hex252(bits250) -> str:
    """T.018 convention: 63 hex characters, two zero bits padded on the left."""
    b = [0, 0] + list(map(int, bits250))
    return "".join("%X" % int("".join(map(str, b[i:i + 4])), 2) for i in range(0, 252, 4))


def hex23_id(bits250) -> str:
    """23-hex beacon ID: T.018 bits 1..202 main field -> ID is built from fixed bits; here we
    return bits 0..91 of the main field as 23 hex characters for grouping."""
    b = list(map(int, bits250[:92]))
    return "".join("%X" % int("".join(map(str, b[i:i + 4])), 2) for i in range(0, 92, 4))
