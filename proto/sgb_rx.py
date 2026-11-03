"""SGB receiver for weak downlink bursts: synthetic modulator, preamble acquisition, coherent
despreading, soft decoding.

Acquisition (PMF-FFT): the 25-bit preamble (6400 chips, 166.7 ms of known chips on I and Q) is
split into 25 one-bit segments. For each coarse frequency bin (150 Hz apart, one segment's
bandwidth) the received block is correlated with every segment replica by FFT at all lags; an
FFT across the 25 segment outputs then adds them coherently with fine frequency resolution. The
statistic is normalised by its median over lags.

Demodulation: chips integrated with fractional timing (cumulative sum + interpolation), timing
and frequency refined on the preamble, then per-bit complex correlations on I and Q. Phase is
tracked decision-directed with a non-causal smoother; soft bits are scaled by the noise
measured in the quadrature component.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.fft as sfft
from scipy import ndimage

import sgb

TC = 1.0 / sgb.CHIP_RATE
NSEG = sgb.PRE_BITS


def baseband(bits250, fs: float, f_off: float = 0.0, phase0: float = 0.0, mode: str = "normal",
             fdot: float = 0.0) -> np.ndarray:
    """Unit-amplitude complex OQPSK burst (rectangular chips), 1.0 s, at sample rate fs.
    Total power 1: I and Q each carry +-1/sqrt(2)."""
    cI, cQ = sgb.chips_pm(mode)
    bI = np.zeros(sgb.BITS_PER_CH, dtype=np.float32)
    bQ = np.zeros(sgb.BITS_PER_CH, dtype=np.float32)
    b = np.asarray(bits250, dtype=np.float32)
    bI[sgb.PRE_BITS:] = b[0::2]
    bQ[sgb.PRE_BITS:] = b[1::2]
    sI = cI * (1 - 2 * np.repeat(bI, sgb.CHIPS_PER_BIT))
    sQ = cQ * (1 - 2 * np.repeat(bQ, sgb.CHIPS_PER_BIT))
    n = int(round((sgb.BITS_PER_CH * sgb.CHIPS_PER_BIT + 0.5) * TC * fs))
    t = np.arange(n) / fs
    ki = np.minimum((t / TC).astype(int), len(sI) - 1)
    kq = np.floor(t / TC - 0.5).astype(int)
    q = np.where((kq >= 0) & (kq < len(sQ)), sQ[np.clip(kq, 0, len(sQ) - 1)], 0.0)
    i = np.where(t / TC < len(sI), sI[ki], 0.0)
    ph = phase0 + 2 * np.pi * (f_off * t + 0.5 * fdot * t * t)
    return ((i + 1j * q) / np.sqrt(2) * np.exp(1j * ph)).astype(np.complex64)


# ---------------- acquisition ----------------

class Acquirer:
    def __init__(self, fs: float, nfft: int, span_hz: float = 10e3, mode: str = "normal"):
        self.fs, self.nfft = fs, nfft
        cI, cQ = sgb.chips_pm(mode)
        nseg_samp = sgb.CHIPS_PER_BIT * TC * fs
        self.seg_T = sgb.CHIPS_PER_BIT * TC
        npre = int(np.ceil(NSEG * nseg_samp)) + 2
        t = np.arange(npre) / fs
        ki = (t / TC).astype(int)
        kq = np.floor(t / TC - 0.5).astype(int)
        G = np.zeros((NSEG, nfft), dtype=np.complex64)
        for j in range(NSEG):
            lo, hi = j * sgb.CHIPS_PER_BIT, (j + 1) * sgb.CHIPS_PER_BIT
            m_i = (ki >= lo) & (ki < hi)
            m_q = (kq >= lo) & (kq < hi)
            g = np.zeros(npre, dtype=np.complex64)
            g[m_i] += cI[ki[m_i]]
            g[m_q] += 1j * cQ[kq[m_q]]
            G[j] = np.conj(np.fft.fft(g, nfft))
        self.G = G
        self.npre = npre
        df = 1.0 / self.seg_T
        self.bins = np.arange(-span_hz, span_hz + df / 2, df)
        self.M = 32

    def search(self, x: np.ndarray, f_center: float = 0.0, top: int = 3):
        """x: complex block at fs (len <= nfft - npre usable lags). Returns list of
        (z, lag_samples, f_hz) best peaks, z = peak / median power over lags (per bin best)."""
        n = len(x)
        X = sfft.fft(x, self.nfft, workers=-1)
        fres = self.fs / self.nfft
        nlag = n - self.npre
        if nlag <= 0:
            return []
        best = []
        jj = np.arange(NSEG)[:, None]
        for fb in self.bins:
            sh = int(round((fb + f_center) / fres))
            Xs = np.roll(X, -sh)
            Z = sfft.ifft(Xs[None, :] * self.G, axis=1, workers=-1)[:, :nlag]   # (25, lags)
            Y = sfft.fft(Z, self.M, axis=0, workers=-1)                         # fine frequency
            P = np.abs(Y) ** 2
            pm = P.max(axis=0)
            k = int(np.argmax(pm))
            med = float(np.median(pm))
            m = int(np.argmax(P[:, k]))
            mf = m if m < self.M // 2 else m - self.M
            f_fine = fb + f_center + mf / (self.M * self.seg_T)
            best.append((pm[k] / med, k, f_fine))
        best.sort(key=lambda b: -b[0])
        return best[:top]


# ---------------- demodulation ----------------

def chip_integrals(x: np.ndarray, fs: float, t0: float, nchips: int, offset_chips: float = 0.0,
                   rate_scale: float = 1.0) -> np.ndarray:
    """Integrate x over chip intervals starting at t0 (s) + offset_chips*Tc, chip period
    Tc*rate_scale, using the cumulative sum with linear interpolation."""
    cs = np.concatenate([[0], np.cumsum(x.astype(np.complex128))])
    edges = (t0 + (np.arange(nchips + 1) + offset_chips) * TC * rate_scale) * fs
    edges = np.clip(edges, 0, len(x))
    n = np.arange(len(cs))
    re = np.interp(edges, n, cs.real)
    im = np.interp(edges, n, cs.imag)
    v = re + 1j * im
    return (v[1:] - v[:-1]).astype(np.complex128)


@dataclass
class SgbDemod:
    ok: bool
    llr: np.ndarray | None = None
    t0: float = 0.0
    f_hz: float = 0.0
    mode: str = "normal"
    ebn0_db: float = -99.0
    pre_z: float = 0.0
    fdot: float = 0.0
    notes: list[str] = field(default_factory=list)


def demodulate(x: np.ndarray, fs: float, t0_guess: float, f_guess: float, mode: str = "normal",
               smooth_bits: int = 61) -> SgbDemod:
    """x: complex baseband (carrier near f_guess), t0_guess: start of the first chip (s)."""
    d = SgbDemod(ok=False, mode=mode)
    cI, cQ = sgb.chips_pm(mode)
    npc = sgb.PRE_BITS * sgb.CHIPS_PER_BIT
    tt = np.arange(len(x)) / fs

    def pre_corr(t0, f):
        y = x * np.exp(-2j * np.pi * f * tt)
        a = chip_integrals(y, fs, t0, npc)
        b = chip_integrals(y, fs, t0, npc, 0.5)
        # per-bit correlation, both channels, then coherent sum
        zi = (a * cI[:npc]).reshape(sgb.PRE_BITS, -1).sum(1)
        zq = (b * cQ[:npc] * -1j).reshape(sgb.PRE_BITS, -1).sum(1)
        return zi + zq

    # refine timing (+-0.6 chip) and frequency (+-20 Hz) on the preamble
    best = (-1, t0_guess, f_guess)
    for dtc in np.arange(-1.5, 1.51, 0.25):
        for dfh in np.arange(-20, 20.1, 4):
            z = pre_corr(t0_guess + dtc * TC, f_guess + dfh)
            # linear phase trend across bits = residual frequency; use FFT magnitude over bits
            p = np.abs(np.fft.fft(z, 64)).max()
            if p > best[0]:
                best = (p, t0_guess + dtc * TC, f_guess + dfh)
    _, t0, f = best
    z = pre_corr(t0, f)
    # residual frequency from the phase slope across the preamble bits
    ph = np.unwrap(np.angle(z))
    slope = np.polyfit(np.arange(len(z)), ph, 1)[0]
    f += slope / (2 * np.pi * sgb.CHIPS_PER_BIT * TC)
    for dtc in np.arange(-0.15, 0.151, 0.05):
        zz = pre_corr(t0 + dtc * TC, f)
        if np.abs(zz.sum()) > np.abs(z.sum()):
            z, t0 = zz, t0 + dtc * TC
    d.t0, d.f_hz = t0, f

    # whole burst
    y = x * np.exp(-2j * np.pi * f * tt)
    nch = sgb.BITS_PER_CH * sgb.CHIPS_PER_BIT
    a = chip_integrals(y, fs, t0, nch)
    b = chip_integrals(y, fs, t0, nch, 0.5)
    zi = (a * cI).reshape(sgb.BITS_PER_CH, -1).sum(1)              # data on Re after derotation
    zq = (b * cQ * -1j).reshape(sgb.BITS_PER_CH, -1).sum(1)        # data on Re after derotation
    known = np.zeros(sgb.BITS_PER_CH, dtype=bool)
    known[:sgb.PRE_BITS] = True
    # carrier phase model over the whole burst: squaring removes the data (BPSK on each
    # channel), so frequency residual and drift are fitted on all 300 bit correlations; the
    # 180-degree ambiguity of squaring is resolved on the known preamble.
    tb = (np.arange(sgb.BITS_PER_CH) + 0.5) * sgb.CHIPS_PER_BIT * TC
    tb = tb - tb.mean()
    w = zi ** 2 + zq ** 2
    dfs = np.arange(-3.0, 3.01, 0.1)
    fds = np.arange(-10.0, 10.01, 0.25)
    DF, FD = np.meshgrid(dfs, fds, indexing="ij")
    th = 2 * np.pi * (DF[..., None] * tb + 0.5 * FD[..., None] * tb ** 2)   # (ndf, nfd, 150)
    score = np.abs(np.sum(w * np.exp(-2j * th), axis=-1))
    a, b2 = np.unravel_index(int(np.argmax(score)), score.shape)
    theta = th[a, b2]
    phi0 = np.angle(np.sum((zi[:sgb.PRE_BITS] + zq[:sgb.PRE_BITS]) * np.exp(-1j * theta[:sgb.PRE_BITS])))
    d.f_hz += dfs[a]
    d.fdot = float(fds[b2])
    phase = phi0 + theta
    smooth_bits = min(smooth_bits, 31)
    for _ in range(3):
        rot = np.exp(-1j * phase)
        si = np.where(known, 1.0, np.sign((zi * rot).real))
        sq = np.where(known, 1.0, np.sign((zq * rot).real))
        pil = zi * si + zq * sq
        # smooth only the residual around the fitted model
        res = pil * np.exp(-1j * (phi0 + theta))
        sm = ndimage.uniform_filter1d(res.real, smooth_bits, mode="nearest") + \
            1j * ndimage.uniform_filter1d(res.imag, smooth_bits, mode="nearest")
        phase = phi0 + theta + np.angle(sm)
    rot = np.exp(-1j * phase)
    ri, rq = (zi * rot), (zq * rot)
    # noise from the quadrature component (both channels, message part)
    sig2 = float(np.mean(np.concatenate([ri.imag[sgb.PRE_BITS:], rq.imag[sgb.PRE_BITS:]]) ** 2))
    mu = float(np.mean(np.concatenate([np.abs(ri.real[sgb.PRE_BITS:]), np.abs(rq.real[sgb.PRE_BITS:])])))
    # bit 0 -> chips as PRN -> positive correlation; LLR > 0 means bit 1
    llr = np.empty(sgb.MSG_BITS)
    llr[0::2] = -2 * mu * ri.real[sgb.PRE_BITS:] / sig2
    llr[1::2] = -2 * mu * rq.real[sgb.PRE_BITS:] / sig2
    d.llr = llr
    d.ebn0_db = 10 * np.log10(max(mu ** 2 / (2 * sig2), 1e-6))
    d.pre_z = float(np.abs(z.sum()) / np.sqrt(sig2 * sgb.PRE_BITS * 2))
    d.ok = True
    return d


def decode(d: SgbDemod, L: int = 8):
    bits, cost, flips = sgb.chase(d.llr, L=L)
    return bits, cost, flips
