"""Cut a window from a ci16 recording and resample it to 2.4576 MS/s (dec406's validated rate).
usage: win2457.py SRC t_onset f_off_hz OUT.cs16  -> writes OUT (ci16, x16 scale), prints new t, f"""
import sys
import numpy as np
from scipy import signal
import dsp
src, t, f, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
fs = 2.5e6
pre = 0.4
x = dsp.read_ci16(src, t - pre, 1.6, fs)
# shift the burst near 0 Hz first so the resampler's band holds it, then resample
y = signal.resample_poly(dsp.mix(x, f, fs), 3072, 3125).astype(np.complex64) * 16
if len(sys.argv) > 5:   # optional narrow low-pass (Hz, one side) to keep other relays out
    h = signal.firwin(4001, float(sys.argv[5]), fs=2457600.0)
    y = signal.fftconvolve(y, h, mode="same").astype(np.complex64)
o = np.empty(2 * len(y), dtype=np.int16)
o[0::2] = np.clip(np.round(y.real), -32768, 32767); o[1::2] = np.clip(np.round(y.imag), -32768, 32767)
o.tofile(out)
print(f"{pre:.4f} 0.0")
