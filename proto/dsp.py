"""Signal helpers: ci16 reading, channelising, FGB modulation for synthetic tests."""
from __future__ import annotations

import numpy as np
from scipy import signal

import fgb


def read_ci16(path, start_s: float, dur_s: float, fs: float) -> np.ndarray:
    n0, n = int(round(start_s * fs)), int(round(dur_s * fs))
    raw = np.fromfile(path, dtype=np.int16, count=2 * n, offset=4 * n0)
    return (raw[0::2].astype(np.float32) + 1j * raw[1::2].astype(np.float32)).astype(np.complex64)


def mix(x: np.ndarray, f_hz: float, fs: float, t0_samples: int = 0) -> np.ndarray:
    """Shift by -f_hz (brings f_hz to 0). Phase continuous across calls through t0_samples."""
    n = np.arange(len(x), dtype=np.float64) + t0_samples
    return (x * np.exp(-2j * np.pi * f_hz / fs * n).astype(np.complex64)).astype(np.complex64)


_fir_cache: dict = {}


def channelise(x: np.ndarray, f_hz: float, fs: float, decim: int, bw_hz: float,
               t0_samples: int = 0) -> np.ndarray:
    """Mix f_hz to 0, low-pass to +-bw_hz/2, decimate. Linear phase FIR, group delay removed."""
    key = (fs, decim, bw_hz)
    if key not in _fir_cache:
        ntaps = 8 * decim + 1
        _fir_cache[key] = signal.firwin(ntaps, bw_hz / 2, fs=fs, window=("kaiser", 8.0)).astype(np.float32)
    h = _fir_cache[key]
    y = signal.resample_poly(mix(x, f_hz, fs, t0_samples), 1, decim, window=h)
    return y.astype(np.complex64)


def fgb_baseband(bits, fs: float, f_off: float = 0.0, fdot: float = 0.0, phase0: float = 0.0,
                 cw_s: float = fgb.CW_SECONDS, polarity: int = 1) -> np.ndarray:
    """Complex FGB burst at unit amplitude: CW then biphase-L +-1.1 rad, with frequency
    offset f_off and linear drift fdot (Hz/s). polarity=+1: bit 1 -> +phase then -phase."""
    nb = len(bits)
    ncw = int(round(cw_s * fs))
    nbit = fs / fgb.BIT_RATE
    ndat = int(round(nb * nbit))
    m = np.zeros(ncw + ndat, dtype=np.float32)
    t = np.arange(ndat) / nbit  # time in bits
    k = np.minimum(t.astype(int), nb - 1)
    first_half = (t - k) < 0.5
    b = np.asarray(bits)[k] * 2 - 1
    m[ncw:] = np.where(first_half, b, -b) * polarity
    tt = np.arange(len(m)) / fs
    ph = phase0 + 2 * np.pi * (f_off * tt + 0.5 * fdot * tt ** 2) + fgb.MOD_INDEX * m
    return np.exp(1j * ph).astype(np.complex64)


def noise_psd(x: np.ndarray, fs: float, f_lo: float, f_hi: float, nfft: int = 8192) -> float:
    """Median-based noise PSD (power per Hz) inside [f_lo, f_hi] of baseband x."""
    f, p = signal.welch(x, fs=fs, nperseg=nfft, return_onesided=False, scaling="density",
                        average="median")
    m = (f >= f_lo) & (f < f_hi)
    return float(np.median(p[m]))
