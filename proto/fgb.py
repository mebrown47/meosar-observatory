"""First-generation (FGB) 406 MHz beacon frame: layout, BCH codes, soft decoding.

Frame (C/S T.001), bit indices 0-based as transmitted:
  0..14    bit sync, 15 ones
  15..23   frame sync: 000101111 normal, 011010000 self-test
  24       format flag (0 short 112 bits, 1 long 144 bits)
  24..84   PDF-1 (61 bits)      85..105  BCH-1 parity (21 bits)   -> BCH(82,61), t=3
  106..131 PDF-2 (26 bits)      132..143 BCH-2 parity (12 bits)   -> BCH(38,26), t=2 (long frames)
A codeword is valid when its bits, read as a polynomial with the first bit as the highest power,
are divisible by the generator.
"""
from __future__ import annotations

import itertools
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

SYNC_NORMAL = (0, 0, 0, 1, 0, 1, 1, 1, 1)
SYNC_SELFTEST = (0, 1, 1, 0, 1, 0, 0, 0, 0)
PREAMBLE = (1,) * 15
BIT_RATE = 400.0
CW_SECONDS = 0.160
MOD_INDEX = 1.1  # rad

G1 = int("1001101101100111100011", 2)  # degree 21
G2 = int("1010100111001", 2)           # degree 12


def bits_to_int(bits) -> int:
    v = 0
    for b in bits:
        v = (v << 1) | int(b)
    return v


def int_to_bits(v: int, n: int) -> list[int]:
    return [(v >> (n - 1 - i)) & 1 for i in range(n)]


def polymod(v: int, nbits: int, g: int) -> int:
    deg = g.bit_length() - 1
    for i in range(nbits - 1, deg - 1, -1):
        if v >> i & 1:
            v ^= g << (i - deg)
    return v


class BchCode:
    """Shortened binary BCH code: syndrome-table hard decoder and Chase-II soft decoder."""

    def __init__(self, g: int, n: int, t: int):
        self.g, self.n, self.t = g, n, t
        self.r = g.bit_length() - 1
        self.k = n - self.r
        # remainder of a single error at position i (i = 0 is the first, highest-power bit)
        self.pos_rem = np.array([polymod(1 << (n - 1 - i), n, g) for i in range(n)], dtype=np.int64)
        table: dict[int, tuple[int, ...]] = {0: ()}
        for w in range(1, t + 1):
            for pos in itertools.combinations(range(n), w):
                s = 0
                for p in pos:
                    s ^= int(self.pos_rem[p])
                table.setdefault(s, pos)
        self.table = table

    def encode(self, msg_bits) -> list[int]:
        m = bits_to_int(msg_bits) << self.r
        return list(msg_bits) + int_to_bits(polymod(m, self.n, self.g), self.r)

    def syndrome(self, bits) -> int:
        s = 0
        for i, b in enumerate(bits):
            if b:
                s ^= int(self.pos_rem[i])
        return s

    def hard_decode(self, bits):
        """Returns (corrected bits, n corrected) or (None, -1)."""
        pos = self.table.get(self.syndrome(bits))
        if pos is None:
            return None, -1
        out = list(bits)
        for p in pos:
            out[p] ^= 1
        return out, len(pos)

    def chase(self, llr: np.ndarray, L: int = 10):
        """Chase-II: flip every subset of the L least reliable bits, hard-decode each test
        pattern, keep the codeword closest to the soft input. llr > 0 means bit 1.
        Returns (bits, cost, n_flipped_vs_hard, n_distinct) or (None, inf, -1, 0).
        cost = sum of |llr| over bits where the codeword disagrees with the hard decision."""
        hard = (llr > 0).astype(np.int64)
        rel = np.abs(llr)
        s0 = int(np.bitwise_xor.reduce(self.pos_rem[hard == 1])) if hard.any() else 0
        weak = np.argsort(rel)[:L]
        best, best_cost, seen = None, np.inf, set()
        for mask in range(1 << L):
            flips = [int(weak[j]) for j in range(L) if mask >> j & 1]
            s = s0
            for p in flips:
                s ^= int(self.pos_rem[p])
            pos = self.table.get(s)
            if pos is None:
                continue
            e = set(flips) ^ set(pos)  # total positions flipped vs hard decision
            key = frozenset(e)
            if key in seen:
                continue
            seen.add(key)
            cost = float(sum(rel[p] for p in e))
            if cost < best_cost:
                best_cost, best = cost, e
        if best is None:
            return None, np.inf, -1, 0
        out = hard.copy()
        for p in best:
            out[p] ^= 1
        return out.tolist(), best_cost, len(best), len(seen)


BCH1 = BchCode(G1, 82, 3)
BCH2 = BchCode(G2, 38, 2)


def _load_mids() -> set[int]:
    p = Path.home() / "dec406/include/country_codes.h"  # ITU MID table (facts, used as data)
    try:
        return {int(m) for m in re.findall(r"\{\s*(\d{3})\s*,", p.read_text())}
    except OSError:
        return set(range(201, 776))


MIDS = _load_mids()


def frame_hex(bits) -> str:
    """Hex of bits 25..144 (1-based), the usual 30-hex-char FGB message (15 for short)."""
    n = len(bits)
    payload = bits[24:n]
    payload = list(payload) + [0] * (-len(payload) % 4)
    return "".join("%X" % bits_to_int(payload[i:i + 4]) for i in range(0, len(payload), 4))


def hex_id15(bits) -> str:
    """15-hex beacon ID: bits 26..85 (1-based), with location bits defaulted by protocol is not
    applied here; this is the raw PDF-1 identification field for grouping."""
    b = list(bits[25:85])
    return "".join("%X" % bits_to_int(b[i:i + 4]) for i in range(0, 60, 4))


def country(bits) -> int:
    return bits_to_int(bits[26:36])


def is_orbitography(bits) -> bool:
    # same test as dec406: protocol flag 1 (user protocol) and protocol code 000
    return bool(bits[25]) and not bits[36] and not bits[37] and not bits[38]


@dataclass
class FrameDecode:
    ok: bool
    bits: list[int] | None = None
    long: bool = False
    selftest: bool = False
    bch1_cost: float = float("inf")
    bch1_flips: int = -1
    bch2_cost: float = 0.0
    bch2_flips: int = 0
    reasons: list[str] = field(default_factory=list)

    @property
    def hex(self) -> str:
        return frame_hex(self.bits) if self.bits else ""


def decode_soft(llr: np.ndarray, L: int = 10, max_cost: float | None = None) -> FrameDecode:
    """Soft-decode a 144-value LLR vector (positive = 1). Uses Chase-II on BCH-1 and, for long
    frames, BCH-2. Sync is checked by soft correlation; the caller sets max_cost from noise
    calibration."""
    llr = np.asarray(llr, dtype=float)
    hard = (llr > 0).astype(int)
    d = FrameDecode(ok=False)
    sync = hard[15:24]
    dn = sum(int(a != b) for a, b in zip(sync, SYNC_NORMAL))
    ds = sum(int(a != b) for a, b in zip(sync, SYNC_SELFTEST))
    d.selftest = ds < dn
    if min(dn, ds) > 1:
        d.reasons.append(f"sync errors {min(dn, ds)}")
    c1, cost1, f1, _ = BCH1.chase(llr[24:106], L=L)
    if c1 is None:
        d.reasons.append("bch1 no codeword")
        return d
    d.bch1_cost, d.bch1_flips = cost1, f1
    bits = list(hard[:24]) + list(c1)
    d.long = bool(c1[0])
    if d.long:
        if len(llr) < 144:
            d.reasons.append("long frame truncated")
            return d
        if is_orbitography(bits):
            bits += list(hard[106:144])  # bits 107-144 not BCH-2 protected for orbitography
        else:
            c2, cost2, f2, _ = BCH2.chase(llr[106:144], L=min(L, 8))
            if c2 is None:
                d.reasons.append("bch2 no codeword")
                return d
            d.bch2_cost, d.bch2_flips = cost2, f2
            bits += list(c2)
    else:
        bits += list(hard[106:112])
    d.bits = bits
    mid = country(bits)
    if mid not in MIDS:
        d.reasons.append(f"country {mid} not assigned")
    if max_cost is not None and d.bch1_cost + d.bch2_cost > max_cost:
        d.reasons.append(f"cost {d.bch1_cost + d.bch2_cost:.1f} > {max_cost:.1f}")
    d.ok = not d.reasons
    return d


def build_frame(payload_hex: str) -> list[int]:
    """Frame bits from a 15/30-hex payload (bits 25..), recomputing both BCH fields."""
    v = bin(int(payload_hex, 16))[2:].zfill(len(payload_hex) * 4)
    p = [int(c) for c in v]
    long = bool(p[0])
    pdf1 = p[:61]
    bits = list(PREAMBLE) + list(SYNC_NORMAL) + BCH1.encode(pdf1)
    if long:
        bits += BCH2.encode(p[82:108])
    else:
        bits += p[82:88] + [0] * max(0, 6 - len(p[82:88]))
    return bits
