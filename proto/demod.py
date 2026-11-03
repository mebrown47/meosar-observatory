"""Coherent FGB demodulator for weak bursts, working on a whole burst held in memory.

Model after channelising to fs (~20 kS/s) around the burst carrier:
    x(t) = A exp(j(theta(t) + 1.1 m(t))) + n(t),   m(t) = Manchester +-1, 0 during the 160 ms CW
         = A cos1.1 e^{j theta}  +  j A sin1.1 m(t) e^{j theta}
The carrier never disappears (residual carrier at -6.9 dB during data), so carrier phase is an
absolute reference and data sit in quadrature: after removing theta, Im(x) carries m(t).

Steps
 1. frequency + drift: dechirp over a grid of drift rates, zero-padded FFT, take the carrier peak.
 2. timing: correlate Im(x) with the 24 known bits (15 ones + frame sync), both sync words.
 3. phase: decision-directed, non-causal smoother (remove decided modulation, average the
    complex carrier over a centred window), iterate.
 4. soft bits: half-bit difference (matched filter for biphase-L), LLR = 2 mu r / sigma^2 with
    sigma^2 from the same statistic on the noise before the burst.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage

import fgb

KNOWN = np.array(fgb.PREAMBLE + fgb.SYNC_NORMAL)
KNOWN_ST = np.array(fgb.PREAMBLE + fgb.SYNC_SELFTEST)


def manchester(bits, spb: float, n: int | None = None) -> np.ndarray:
    nb = len(bits)
    n = n or int(round(nb * spb))
    t = np.arange(n) / spb
    k = np.minimum(t.astype(int), nb - 1)
    b = 2 * np.asarray(bits)[k] - 1
    return np.where((t - k) < 0.5, b, -b).astype(np.float32)


@dataclass
class Demod:
    ok: bool
    llr: np.ndarray | None = None
    f_hz: float = 0.0          # carrier offset in the input
    fdot: float = 0.0          # Hz/s
    t_data: float = 0.0        # data start (first preamble bit), s from window start
    t_cw: float = 0.0          # estimated CW start, s
    polarity: int = 1
    selftest_sync: bool = False
    sync_score: float = 0.0    # normalised correlation with the 24 known bits
    cn0_dbhz: float = 0.0
    ebn0_db: float = 0.0
    notes: list[str] = field(default_factory=list)


def estimate_carrier(x: np.ndarray, fs: float, fmax: float = 80.0,
                     fdots=np.arange(-12.0, 12.01, 1.5), pad: int = 4):
    """Max over drift grid of the zero-padded spectrum peak within +-fmax. Returns
    (f, fdot, complex amplitude at t=0 of the window)."""
    n = len(x)
    t = np.arange(n) / fs
    nfft = 1 << int(np.ceil(np.log2(n * pad)))
    freqs = np.fft.fftfreq(nfft, 1 / fs)
    sel = np.abs(freqs) <= fmax
    best = (-1.0, 0.0, 0.0)
    for fd in fdots:
        X = np.fft.fft(x * np.exp(-1j * np.pi * fd * t ** 2), nfft)
        p = np.abs(X[sel]) ** 2
        k = int(np.argmax(p))
        if p[k] > best[0]:
            best = (p[k], float(freqs[sel][k]), float(fd))
    _, f0, fd = best
    # fine: golden-ish local search on f and fdot maximising |sum|
    def score(f, d):
        return np.abs(np.sum(x * np.exp(-2j * np.pi * (f * t + 0.5 * d * t ** 2))))
    step_f, step_d = fs / nfft, 0.75
    for _ in range(25):
        moved = False
        for df, dd in ((step_f, 0), (-step_f, 0), (0, step_d), (0, -step_d)):
            if score(f0 + df, fd + dd) > score(f0, fd):
                f0, fd, moved = f0 + df, fd + dd, True
        if not moved:
            step_f, step_d = step_f / 2, step_d / 2
            if step_f < 0.02:
                break
    a = np.sum(x * np.exp(-2j * np.pi * (f0 * t + 0.5 * fd * t ** 2))) / n
    return f0, fd, a


def demodulate(x: np.ndarray, fs: float, t_cw_guess: float, noise: np.ndarray | None = None,
               search_s: float = 0.08, smooth_s: float = 0.030, iters: int = 3,
               nbits: int = 144) -> Demod:
    """x: complex baseband around the carrier (|offset| < ~80 Hz), long enough to hold
    [t_cw_guess - search_s, t_cw_guess + 0.16 + 0.36 + search_s]. noise: burst-free samples
    from the same channel for sigma; if None, the first 40 ms of x before the CW is used."""
    spb = fs / fgb.BIT_RATE
    d = Demod(ok=False)
    # 1. carrier: use the region that should hold the burst
    i0 = max(0, int((t_cw_guess - search_s) * fs))
    i1 = min(len(x), int((t_cw_guess + fgb.CW_SECONDS + nbits / fgb.BIT_RATE + search_s) * fs))
    seg = x[i0:i1]
    f0, fd, _ = estimate_carrier(seg, fs)
    t = (np.arange(len(x)) - i0) / fs
    z = x * np.exp(-2j * np.pi * (f0 * t + 0.5 * fd * t ** 2))
    d.f_hz, d.fdot = f0, fd

    # constant phase from the burst region (carrier is real-positive after this)
    z = z * np.exp(-1j * np.angle(np.sum(z[i0:i1])))

    # 2. timing via the 24 known bits on Im(z), both sync words and polarities
    ref_n = manchester(KNOWN, spb)
    ref_s = manchester(KNOWN_ST, spb)
    lo = max(0, int((t_cw_guess + fgb.CW_SECONDS - search_s) * fs))
    hi = min(len(z) - len(ref_n) - int(nbits * spb) // 2,
             int((t_cw_guess + fgb.CW_SECONDS + search_s) * fs))
    if hi <= lo:
        d.notes.append("window too short")
        return d
    y = z.imag.astype(np.float64)
    cn = np.correlate(y[lo:hi + len(ref_n)], ref_n, mode="valid")
    cs = np.correlate(y[lo:hi + len(ref_s)], ref_s, mode="valid")
    cand = [(np.max(np.abs(cn)), int(np.argmax(np.abs(cn))), cn, False),
            (np.max(np.abs(cs)), int(np.argmax(np.abs(cs))), cs, True)]
    _, k, c, st = max(cand, key=lambda q: q[0])
    pol = 1 if c[k] > 0 else -1
    n_data = lo + k
    d.polarity, d.selftest_sync = pol, st

    known = KNOWN_ST if st else KNOWN
    nd = int(round(nbits * spb))
    tmpl_idx = n_data + np.arange(nd)
    tmpl_idx = tmpl_idx[tmpl_idx < len(z)]

    def bit_stats(zz, nstart):
        yy = zz.imag.astype(np.float64) * pol
        r = np.zeros(nbits)
        for b in range(nbits):
            a = nstart + b * spb
            h0, h1, h2 = int(round(a)), int(round(a + spb / 2)), int(round(a + spb))
            if h2 > len(yy):
                break
            r[b] = yy[h0:h1].sum() - yy[h1:h2].sum()
        return r

    # 3. decision-directed phase smoothing
    ncw0 = max(0, n_data - int(fgb.CW_SECONDS * fs))
    w = max(3, int(smooth_s * fs)) | 1
    for _ in range(iters):
        r = bit_stats(z, n_data)
        hard = (r > 0).astype(int)
        hard[:24] = known
        m = np.zeros(len(z), dtype=np.float32)
        mm = manchester(hard, spb, nd)[: len(tmpl_idx)] * pol
        m[tmpl_idx] = mm
        c = z * np.exp(-1j * fgb.MOD_INDEX * m)
        # weight: CW samples carry full carrier power, data samples carry cos(1.1) after removal
        # (after removing decided modulation, data samples carry the full carrier as well)
        region = np.zeros(len(z), dtype=bool)
        region[ncw0:tmpl_idx[-1] + 1] = True
        cc = np.where(region, c, 0)
        sm = ndimage.uniform_filter1d(cc.real, w) + 1j * ndimage.uniform_filter1d(cc.imag, w)
        ph = np.angle(sm)
        # hold the edges
        ph[: ncw0] = ph[ncw0]
        ph[tmpl_idx[-1] + 1:] = ph[tmpl_idx[-1]]
        z = z * np.exp(-1j * ph)

    r = bit_stats(z, n_data)
    # amplitude from the CW (real part), noise from the burst-free region
    cw = z[ncw0 + int(0.01 * fs): n_data - int(0.01 * fs)]
    A = float(np.mean(cw.real)) if len(cw) > 10 else float(np.abs(np.mean(z[tmpl_idx])) / np.cos(1.1))
    if noise is None:
        nz = x[: max(0, ncw0 - int(0.02 * fs))]
    else:
        nz = noise
    if len(nz) < 4 * spb:
        d.notes.append("no noise reference")
        return d
    yn = nz.imag.astype(np.float64)
    nsb = int((len(yn) - spb) // spb)
    rn = np.array([yn[int(i * spb): int(i * spb + spb / 2)].sum() - yn[int(i * spb + spb / 2): int((i + 1) * spb)].sum()
                   for i in range(nsb)])
    # also half-bit-shifted samples for more estimates
    sig2 = float(np.mean(rn ** 2))
    mu = A * np.sin(fgb.MOD_INDEX) * spb
    d.llr = 2 * mu * r / sig2
    n0 = float(np.mean(np.abs(nz) ** 2)) / fs * 1.0  # power per Hz if noise were white over fs
    # better N0 from the bit statistic: sigma^2 = N0 * spb / fs... use the measured Eb/N0
    d.ebn0_db = 10 * np.log10(mu ** 2 / (2 * sig2))
    d.cn0_dbhz = d.ebn0_db + 10 * np.log10(fgb.BIT_RATE) - 10 * np.log10(np.sin(fgb.MOD_INDEX) ** 2)
    d.t_data = n_data / fs
    d.t_cw = d.t_data - fgb.CW_SECONDS
    kn = r[:24] * (2 * known - 1)
    d.sync_score = float(np.mean(kn) / np.sqrt(sig2))
    d.ok = True
    return d
