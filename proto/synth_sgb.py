"""Inject synthetic SGB bursts at a known C/N0 into real recorded noise (the quiet band).
usage: synth_sgb.py SRC OUT_BASE --cn0 28 [--n 6] [--spacing 3.0]"""
import argparse, json
from pathlib import Path
import numpy as np
import dsp, sgb, sgb_rx
from run import QUIET, meta_of
ap = argparse.ArgumentParser(); ap.add_argument("src"); ap.add_argument("out")
ap.add_argument("--cn0", type=float, required=True); ap.add_argument("--n", type=int, default=6)
ap.add_argument("--spacing", type=float, default=3.0); ap.add_argument("--src-start", type=float, default=0.0)
a = ap.parse_args()
fs, fc, _ = meta_of(Path(a.src)); rng = np.random.default_rng(int(a.cn0 * 10))
x2 = dsp.read_ci16(a.src, a.src_start, 2.0, fs)
xb = dsp.channelise(x2, (QUIET[0] + QUIET[1]) / 2 - fc, fs, 20, 110e3)   # 110 kHz at 125 kS/s: no aliasing
n0 = dsp.noise_psd(xb, fs / 20, -40e3, 40e3); amp = np.sqrt(10 ** (a.cn0 / 10) * n0)
mid = (QUIET[0] + QUIET[1]) / 2
dur = a.n * a.spacing + 2.0
x = dsp.read_ci16(a.src, a.src_start, dur, fs).astype(np.complex64)
truth = []
for k in range(a.n):
    mids = sorted(__import__("fgb").MIDS)
    info = list(rng.integers(0, 2, 202)); info[30:40] = [int(c) for c in format(int(rng.choice(mids)), "010b")]
    cw = sgb.encode(info)
    t0 = 0.5 + k * a.spacing + rng.uniform(0, 0.3); f = mid + rng.uniform(-7e3, 7e3)
    s = sgb_rx.baseband(cw, fs, f_off=f - fc, phase0=rng.uniform(0, 6.28), fdot=rng.uniform(-5, 5))
    i0 = int(t0 * fs); x[i0:i0 + len(s)] += (amp * s).astype(np.complex64)
    truth.append({"t": round(t0, 6), "f_hz": round(f, 1), "hex": sgb.to_hex252(cw)})
o = np.empty(2 * len(x), dtype=np.int16); o[0::2] = np.clip(np.round(x.real), -32768, 32767); o[1::2] = np.clip(np.round(x.imag), -32768, 32767)
o.tofile(a.out + ".cs16"); json.dump({"cn0": a.cn0, "bursts": truth}, open(a.out + ".truth.json", "w"))
print(f"{a.out}: {a.n} SGB bursts at {a.cn0} dB-Hz")
