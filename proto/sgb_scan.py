"""Search a recording for SGB bursts in one SAR band: PMF-FFT acquisition per 1 s step, then
demodulate and decode every peak above --zthr.
usage: sgb_scan.py FILE.sigmf-data [--start S] [--dur S] [--band galileo] > out.jsonl"""
import argparse, datetime as dt, json, sys, time
from pathlib import Path
import numpy as np
import dsp, sgb, sgb_rx
from run import BANDS, meta_of, D1

ap = argparse.ArgumentParser()
ap.add_argument("file"); ap.add_argument("--start", type=float, default=0.0)
ap.add_argument("--dur", type=float, default=60.0); ap.add_argument("--band", default="galileo")
ap.add_argument("--zthr", type=float, default=5.0); ap.add_argument("--prez", type=float, default=7.0); ap.add_argument("--span", type=float, default=10e3)
a = ap.parse_args()
path = Path(a.file); fs, fc, t0 = meta_of(path)
lo, hi = BANDS[a.band]; mid = (lo + hi) / 2
f1 = fs / D1
acq = sgb_rx.Acquirer(f1, 1 << 18, span_hz=a.span)
step, extra = 1.0, 0.25 + 1.2   # search lags within 1 s; keep the whole burst for demod
t = a.start; tic = time.time(); n_dec = 0
while t < a.start + a.dur:
    x = dsp.read_ci16(path, t, step + extra, fs)
    y = dsp.channelise(x, mid - fc, fs, D1, (hi - lo) * 1.05)
    nsearch = int((step + sgb_rx.NSEG * 256 / 38400 + 0.01) * f1)
    peaks = acq.search(y[:nsearch], top=6)
    relays = []
    for z, lag, f in peaks:
        if z < a.zthr or lag / f1 >= step:
            continue
        d = sgb_rx.demodulate(y, f1, lag / f1, f)
        if d.pre_z < a.prez:
            continue
        hb, nerr = sgb.bm_decode(list((d.llr > 0).astype(int)))
        relays.append((d, float(z), hb, nerr))
    # relays of one burst: within 40 ms, different carriers
    used = set()
    for i, (d, z, hb, nerr) in enumerate(relays):
        if i in used:
            continue
        grp = [i] + [j for j in range(i + 1, len(relays)) if j not in used
                     and abs(relays[j][0].t0 - d.t0) < 0.040 and abs(relays[j][0].f_hz - d.f_hz) > 20]
        used.update(grp)
        L = np.sum([relays[k][0].llr for k in grp], axis=0)
        hb2, n2 = sgb.bm_decode(list((L > 0).astype(int)))
        if hb2 is None:
            b, cost, fl = sgb.chase(L, L=8)
            ok = b is not None and fl <= 8
            bits, how = (b, f"chase {fl} flips") if ok else (None, "fail")
        else:
            bits, how = hb2, f"hard {n2} errors"
        rec = {"t": round(float(t + d.t0), 6), "band": a.band, "relays": [
                   {"f_hz": round(float(mid + relays[k][0].f_hz), 1), "dt_ms": round(1e3 * (relays[k][0].t0 - d.t0), 2),
                    "ebn0": round(float(relays[k][0].ebn0_db), 2), "pre_z": round(float(relays[k][0].pre_z), 1),
                    "alone": relays[k][2] is not None} for k in grp],
               "ok": bits is not None, "how": how, "hex": sgb.to_hex252(bits) if bits else ""}
        if t0 is not None:
            rec["utc"] = (t0 + dt.timedelta(seconds=rec["t"])).isoformat(timespec="milliseconds")
        n_dec += rec["ok"]
        print(json.dumps(rec), flush=True)
    t += step
    print(f"[{t:.0f}/{a.start + a.dur:.0f}] decodes {n_dec}  {(t - a.start) / (time.time() - tic):.2f}x real time", file=sys.stderr, flush=True)
