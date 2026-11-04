"""Where does the FGB demodulator lose near threshold? Demodulate synthetic bursts with oracle
knowledge (true timing / frequency / drift) and compare bit error rates.
usage: oracle_fgb.py BASE   (BASE.cs16 + BASE.truth.json from synth.py)"""
import json, sys
import numpy as np
from scipy import ndimage
import dsp, fgb, demod
from math import erfc, sqrt
base = sys.argv[1]
T = json.load(open(base + ".truth.json")); fs, fc = T["fs"], T["fc"]
QM = 1543.75e6; f2 = 125000 / 6
def bits_of(h): return np.array([1] * 15 + list(fgb.SYNC_NORMAL) + [int(c) for c in bin(int(h, 16))[2:].zfill(120)])
def ber_from(z, n_data, pol, bits):
    spb = f2 / 400; r = np.zeros(144)
    for b in range(144):
        a = n_data + b * spb; h0, h1, h2 = int(round(a)), int(round(a + spb / 2)), int(round(a + spb))
        r[b] = (z.imag[h0:h1].sum() - z.imag[h1:h2].sum()) * pol
    return np.mean((r > 0).astype(int)[24:] != bits[24:])
res = {k: [] for k in ("default", "oracle f+fdot+t, carrier phase only", "oracle f+fdot+t, + DD smoothing", "default but true fdot")}
for b in T["bursts"]:
    bits = bits_of(b["hex"])
    t_win = b["t_cw"] - 0.25
    x = dsp.read_ci16(base + ".cs16", t_win, 1.05, fs)
    y = dsp.channelise(x, QM - fc, fs, 20, 126e3)
    xb = dsp.channelise(y, b["f_hz"] - QM, 125000, 6, 8000)        # carrier exactly at 0 Hz at t = t_cw
    d = demod.demodulate(xb, f2, 0.25)
    res["default"].append(np.mean((d.llr > 0).astype(int)[24:] != bits[24:]))
    # oracle: remove drift (fdot about the burst start) and carrier phase from the CW; true data start
    t = np.arange(len(xb)) / f2 - 0.25
    z = xb * np.exp(-1j * np.pi * b["fdot"] * t ** 2)
    i_cw, i_d = int(0.25 * f2), int((0.25 + 0.160) * f2)
    z = z * np.exp(-1j * np.angle(z[i_cw + 100:i_d - 100].sum()))
    res["oracle f+fdot+t, carrier phase only"].append(ber_from(z, i_d, 1, bits))
    # + decision-directed smoothing like the demodulator (30 ms window), true timing
    zz = z.copy()
    for _ in range(3):
        spb = f2 / 400; r = np.zeros(144)
        for k in range(144):
            a = i_d + k * spb; h0, h1, h2 = int(round(a)), int(round(a + spb / 2)), int(round(a + spb))
            r[k] = zz.imag[h0:h1].sum() - zz.imag[h1:h2].sum()
        hard = (r > 0).astype(int); hard[:24] = bits[:24]
        m = np.zeros(len(zz)); mm = demod.manchester(hard, spb, int(round(144 * spb)))
        m[i_d:i_d + len(mm)] = mm[:len(zz) - i_d]
        c = zz * np.exp(-1j * 1.1 * m); reg = np.zeros(len(zz), bool); reg[i_cw:i_d + len(mm)] = True
        cc = np.where(reg, c, 0); w = int(0.03 * f2) | 1
        sm = ndimage.uniform_filter1d(cc.real, w) + 1j * ndimage.uniform_filter1d(cc.imag, w)
        zz = zz * np.exp(-1j * np.angle(sm))
    res["oracle f+fdot+t, + DD smoothing"].append(ber_from(zz, i_d, 1, bits))
    # default demodulator but with the true drift removed first
    d2 = demod.demodulate(xb * np.exp(-1j * np.pi * b["fdot"] * t ** 2).astype(np.complex64), f2, 0.25)
    res["default but true fdot"].append(np.mean((d2.llr > 0).astype(int)[24:] != bits[24:]))
eb = T["cn0"] - 26.02 - 1.0
print(f"C/N0 {T['cn0']} dB-Hz, Eb/N0 {eb:.2f} dB, theory BER {0.5 * erfc(sqrt(10 ** (eb / 10))):.3f}")
for k, v in res.items():
    v = np.array(v); print(f"  {k:40s} BER median {np.median(v):.3f} mean {v.mean():.3f}  (>10%: {np.sum(v > 0.1)}/{len(v)})")
