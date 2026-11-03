"""Inject synthetic FGB bursts at a known C/N0 into real recorded noise (the quiet band).

usage: synth.py SRC.sigmf-data|cs16 OUT_BASE --cn0 30 [--n 30] [--spacing 1.6] [--src-start S]
Writes OUT_BASE.cs16 (ci16, 2.5 MS/s, centre as the source) and OUT_BASE.truth.json.
"""
from __future__ import annotations

import argparse
import json

import numpy as np

import dsp
import fgb
from run import QUIET, meta_of
from pathlib import Path


def random_frame(rng) -> list[int]:
    """Valid long frame: national-location-like PDF-1 with an assigned country, not orbitography."""
    mids = sorted(fgb.MIDS)
    while True:
        pdf1 = [1] + list(rng.integers(0, 2, 60))
        mid = int(rng.choice(mids))
        pdf1[2:12] = fgb.int_to_bits(mid, 10)
        bits = list(fgb.PREAMBLE) + list(fgb.SYNC_NORMAL) + fgb.BCH1.encode(pdf1)
        if fgb.is_orbitography(bits):
            continue
        bits += fgb.BCH2.encode(list(rng.integers(0, 2, 26)))
        return bits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("out")
    ap.add_argument("--cn0", type=float, required=True)
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--spacing", type=float, default=1.6)
    ap.add_argument("--src-start", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--fixed-f", action="store_true", help="all bursts within +-300 Hz of the quiet-band centre (dec406 needs one -f)")
    a = ap.parse_args()
    fs, fc, _ = meta_of(Path(a.src))
    rng = np.random.default_rng(a.seed + int(a.cn0 * 10))
    dur = a.n * a.spacing + 1.5
    # N0 from the quiet band on the first 2 s of the noise
    x2 = dsp.read_ci16(a.src, a.src_start, 2.0, fs)
    xb = dsp.channelise(x2, (QUIET[0] + QUIET[1]) / 2 - fc, fs, 20, 260e3)
    n0 = dsp.noise_psd(xb, fs / 20, -100e3, 100e3)   # unit-DC-gain FIR keeps the PSD per Hz
    del x2, xb
    amp = np.sqrt(10 ** (a.cn0 / 10) * n0)
    truth, waves = [], []
    for k in range(a.n):
        bits = random_frame(rng)
        t_cw = 0.6 + k * a.spacing + rng.uniform(-0.05, 0.05)
        f = ((QUIET[0] + QUIET[1]) / 2 + rng.uniform(-300, 300)) if a.fixed_f else rng.uniform(QUIET[0] + 20e3, QUIET[1] - 20e3)
        fdot = rng.uniform(-5, 5)
        waves.append((int(t_cw * fs), f, fdot, rng.uniform(0, 2 * np.pi), bits))
        truth.append({"t_cw": round(t_cw, 4), "f_hz": round(f, 1), "fdot": round(fdot, 2),
                      "hex": fgb.frame_hex(bits), "id": fgb.frame_hex(bits)[:15],
                      "country": fgb.country(bits)})
    blk = int(fs)
    nblk = int(np.ceil(dur * fs / blk))
    with open(a.out + ".cs16", "wb") as fo:
        for b in range(nblk):
            n0s = b * blk
            x = dsp.read_ci16(a.src, a.src_start + n0s / fs, blk / fs, fs)
            for (i0, f, fdot, ph, bits) in waves:
                s = dsp.fgb_baseband(bits, fs, f_off=f - fc, fdot=fdot, phase0=ph)
                lo, hi = max(i0, n0s), min(i0 + len(s), n0s + len(x))
                if lo >= hi:
                    continue
                x[lo - n0s:hi - n0s] += (amp * s[lo - i0:hi - i0]).astype(np.complex64)
            out = np.empty(2 * len(x), dtype=np.int16)
            out[0::2] = np.clip(np.round(x.real), -32768, 32767)
            out[1::2] = np.clip(np.round(x.imag), -32768, 32767)
            out.tofile(fo)
    json.dump({"cn0": a.cn0, "fs": fs, "fc": fc, "n0": n0, "amp": amp, "bursts": truth},
              open(a.out + ".truth.json", "w"), indent=0)
    print(f"{a.out}: {a.n} bursts at C/N0 {a.cn0} dB-Hz, N0 {n0:.3g}/Hz, amplitude {amp:.2f}")


if __name__ == "__main__":
    main()
