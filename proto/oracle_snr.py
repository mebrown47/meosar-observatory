"""Matched-filter statistics of synthetic FGB bursts with oracle phase/timing: signal mean and
noise spread during the burst, vs noise before it, vs the generator's N0."""
import json, sys
import numpy as np
import dsp, fgb, demod
base = sys.argv[1]
T = json.load(open(base + ".truth.json")); fs, fc = T["fs"], T["fc"]
QM = 1543.75e6; f2 = 125000 / 6; spb = f2 / 400
mu_l, sd_l, sdn_l, A_l = [], [], [], []
for b in T["bursts"]:
    bits = np.array([1] * 15 + list(fgb.SYNC_NORMAL) + [int(c) for c in bin(int(b["hex"], 16))[2:].zfill(120)])
    x = dsp.read_ci16(base + ".cs16", b["t_cw"] - 0.25, 1.05, fs)
    y = dsp.channelise(x, QM - fc, fs, 20, 126e3)
    xb = dsp.channelise(y, b["f_hz"] - QM, 125000, 6, 8000)
    t = np.arange(len(xb)) / f2 - 0.25
    z = xb * np.exp(-1j * np.pi * b["fdot"] * t ** 2)
    i_cw, i_d = int(0.25 * f2), int(0.41 * f2)
    ph = np.angle(z[i_cw + 100:i_d - 100].sum()); z = z * np.exp(-1j * ph)
    A_l.append(np.mean(z[i_cw + 100:i_d - 100].real))
    def stat(zz, n0):
        r = []
        for k in range(144):
            a = n0 + k * spb; h0, h1, h2 = int(round(a)), int(round(a + spb / 2)), int(round(a + spb))
            r.append(zz.imag[h0:h1].sum() - zz.imag[h1:h2].sum())
        return np.array(r)
    r = stat(z, i_d) * (2 * bits - 1)
    rn = stat(z, 100)[:30]                        # same statistic on noise before the burst (~0.07 s)
    mu_l.append(r.mean()); sd_l.append(r.std()); sdn_l.append(np.sqrt(np.mean(rn ** 2)))
mu, sd, sdn, A = map(np.mean, (mu_l, sd_l, sdn_l, A_l))
amp = T["amp"]; n0 = T["n0"]
print(f"C/N0 {T['cn0']}: generator amplitude {amp:.3f}, N0 {n0:.4g}/Hz")
print(f"  measured CW amplitude {A:.3f} (x{A / amp:.3f} of generator, after both filters)")
print(f"  bit statistic: mean {mu:.2f}, spread in burst {sd:.2f}, spread on noise before burst {sdn:.2f}")
print(f"  expected mean A*sin(1.1)*spb = {A * np.sin(1.1) * spb:.2f}")
N0_meas = sdn ** 2 / (spb * f2 * 0.5) / f2   # noise var of the half-difference sum: spb samples * var_im; var_im = N0*fs/2
print(f"  N0 from the noise statistic {N0_meas:.4g}/Hz  (generator {n0:.4g}) -> {10*np.log10(N0_meas / n0):+.2f} dB")
print(f"  effective Eb/N0 (in-burst) {10*np.log10(mu**2 / (2*sd**2)):.2f} dB, from noise-only spread {10*np.log10(mu**2/(2*sdn**2)):.2f} dB, theory {T['cn0'] - 27.02:.2f} dB")
