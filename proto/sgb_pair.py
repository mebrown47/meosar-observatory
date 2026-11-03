"""Demodulate given SGB relays (t, f) of one burst, combine their soft bits, decode.
usage: sgb_pair.py FILE t1 f1 [t2 f2 ...]"""
import sys
from pathlib import Path
import numpy as np
import dsp, sgb, sgb_rx
from run import BANDS, meta_of, D1
path = Path(sys.argv[1]); fs, fc, _ = meta_of(path)
lo, hi = BANDS["galileo"]; mid = (lo + hi) / 2; f1 = fs / D1
pairs = [(float(sys.argv[i]), float(sys.argv[i + 1])) for i in range(2, len(sys.argv), 2)]
import os
t0 = float(os.environ.get("T0", min(t for t, _ in pairs) - 0.3))
x = dsp.read_ci16(path, t0, float(os.environ.get("WIN", 2.45)), fs)
y = dsp.channelise(x, mid - fc, fs, D1, (hi - lo) * 1.05)
llrs = []
for t, f in pairs:
    d = sgb_rx.demodulate(y, f1, t - t0, f - mid)
    b, c, fl = sgb_rx.decode(d, L=8)
    hb, n = sgb.bm_decode(list((d.llr > 0).astype(int)))
    print(f"relay t={t} f={f:.1f}: Eb/N0 {d.ebn0_db:.2f} pre_z {d.pre_z:.1f} | hard BCH {'ok, %d errors' % n if hb else 'fail'} | chase flips {fl}")
    llrs.append(d.llr)
L = np.sum(llrs, axis=0)
hb, n = sgb.bm_decode(list((L > 0).astype(int)))
b, c, fl = sgb.chase(L, L=8)
print(f"combined: hard BCH {'ok, %d errors' % n if hb else 'fail'}  chase flips {fl} cost {c:.1f}")
for name, bits in (("hard", hb), ("chase", b)):
    if bits:
        print(f"  {name}: {sgb.to_hex252(bits)}")
