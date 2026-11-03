"""FGB pipeline over a ci16 SigMF recording: coherent CW detector -> demod -> soft decode.

usage: run.py FILE.sigmf-data [--start S] [--dur S] [--thr DB] > out.jsonl
Detector: each SAR band is channelised to 125 kS/s; 131 ms Hann FFTs every 33 ms; a bin whose
power exceeds the per-bin noise mean (median over the block / ln2) by --thr dB is a CW candidate.
Candidates are merged within +-40 Hz / 0.4 s and each goes to the demodulator.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import time
from pathlib import Path

import numpy as np

import demod
import dsp
import fgb

BANDS = {"galileo": (1544.04e6, 1544.16e6), "glonass": (1544.84e6, 1544.96e6)}
D1 = 20          # 2.5 MS/s -> 125 kS/s for detection
D2 = 125         # 2.5 MS/s -> 20 kS/s for demodulation
NFFT, HOP = 16384, 4096


def meta_of(path: Path):
    m = path.with_suffix(".sigmf-meta")
    try:
        j = json.loads(m.read_text())
        cap = j["captures"][0]
        return (float(j["global"]["core:sample_rate"]), float(cap["core:frequency"]),
                dt.datetime.fromisoformat(cap["core:datetime"].replace("Z", "+00:00")))
    except (OSError, KeyError, ValueError):
        return 2.5e6, 1544.5e6, None


def band_stream(x, fs, fc, lo, hi):
    """Channelise one SAR band to fs/D1. Returns (y, f1, mid)."""
    mid = (lo + hi) / 2
    return dsp.channelise(x, mid - fc, fs, D1, (hi - lo) * 1.05), fs / D1, mid


def detect(y, f1, mid, lo, hi, t_y0, thr_db, band):
    """CW candidates in a band stream: list of (t_frame_start_abs, f_abs, snr_db, band)."""
    out = []
    nfr = (len(y) - NFFT) // HOP + 1
    if nfr < 8:
        return out
    win = np.hanning(NFFT).astype(np.float32)
    idx = np.arange(NFFT)[None, :] + HOP * np.arange(nfr)[:, None]
    P = np.abs(np.fft.fftshift(np.fft.fft(y[idx] * win, axis=1), axes=1)) ** 2
    fr = np.fft.fftshift(np.fft.fftfreq(NFFT, 1 / f1)) + mid
    keep = (fr >= lo) & (fr < hi)
    P, fr = P[:, keep], fr[keep]
    R = P / (np.median(P, axis=0) / np.log(2))[None, :]
    for i, j in np.argwhere(R > 10 ** (thr_db / 10)):
        if R[i, j] < R[max(0, i - 1): i + 2, max(0, j - 2): j + 3].max():
            continue
        out.append((t_y0 + i * HOP / f1, float(fr[j]), 10 * np.log10(R[i, j]), band))
    return out


def merge(cands, df=40.0, dtm=0.4):
    cands.sort(key=lambda c: -c[2])
    kept = []
    for c in cands:
        if any(abs(c[1] - k[1]) < df and abs(c[0] - k[0]) < dtm for k in kept):
            continue
        kept.append(c)
    return sorted(kept)


D3 = 6  # band stream 125 kS/s -> 20.8 kS/s for demodulation


def process_candidate(y, f1, mid, t_y0, t_frame, f_abs, L):
    """Demodulate one candidate from a band stream. t_frame = start of its detector frame."""
    t_cw_guess = t_frame + 0.01
    t_win = t_cw_guess - 0.25
    i0 = int(round((t_win - t_y0) * f1))
    i1 = i0 + int(1.05 * f1)
    if i0 < 0 or i1 > len(y):
        return None
    seg = y[i0:i1]
    f2 = f1 / D3
    xb = dsp.channelise(seg, f_abs - mid, f1, D3, 8000)
    s = xb[int(0.25 * f2): int(0.39 * f2)]
    S = np.abs(np.fft.fft(s, 1 << 15))
    fr = np.fft.fftfreq(1 << 15, 1 / f2)
    m = np.abs(fr) < 40
    f_fine = float(fr[m][np.argmax(S[m])])
    xb = dsp.channelise(seg, f_abs - mid + f_fine, f1, D3, 8000)
    d = demod.demodulate(xb, f2, 0.25)
    if not d.ok:
        return {"demod": "fail", "notes": d.notes}
    r = fgb.decode_soft(d.llr, L=L)
    return {
        "t_cw": round(t_win + d.t_cw, 4), "f_hz": round(f_abs + f_fine + d.f_hz, 1),
        "fdot": round(d.fdot, 2), "ebn0": round(d.ebn0_db, 2), "cn0": round(d.cn0_dbhz, 2),
        "sync": round(d.sync_score, 2), "selftest_sync": d.selftest_sync, "pol": d.polarity,
        "ok": r.ok, "hex": r.hex, "long": r.long, "cost": round(r.bch1_cost + r.bch2_cost, 2),
        "flips": [r.bch1_flips, r.bch2_flips], "reasons": r.reasons,
        "hard_hex": fgb.frame_hex(list((d.llr > 0).astype(int))),
        "_llr": d.llr,
    }


def combine_relays(recs, L, dt_max=0.040, df_same=20.0):
    """Group per-relay records of one burst and decode their summed LLRs (maximum-ratio
    combining: each LLR is already scaled by its own SNR). Returns new combined records."""
    pool = [r for r in recs if "_llr" in r and r["sync"] >= 0.8]
    pool.sort(key=lambda r: r["t_cw"])
    out, used = [], set()
    for i, r in enumerate(pool):
        if i in used:
            continue
        grp = [i]
        for j in range(i + 1, len(pool)):
            q = pool[j]
            if q["t_cw"] - r["t_cw"] > dt_max:
                break
            if j in used:
                continue
            # same relay detected twice (sideband lines): same refined carrier -> skip
            if any(abs(q["f_hz"] - pool[k]["f_hz"]) < df_same for k in grp):
                used.add(j)
                continue
            grp.append(j)
        if len(grp) < 2:
            continue
        used.update(grp)
        rs = [pool[k] for k in grp]
        llr = np.sum([q["_llr"] for q in rs], axis=0)
        dec = fgb.decode_soft(llr, L=L)
        eb = 10 * np.log10(sum(10 ** (q["ebn0"] / 10) for q in rs))
        sy = float(np.sqrt(sum(q["sync"] ** 2 for q in rs)))
        out.append({
            "t": rs[0]["t"], "band": rs[0]["band"], "combined": len(rs),
            "relays": [{"f_hz": q["f_hz"], "ebn0": q["ebn0"], "dt_ms": round(1e3 * (q["t_cw"] - rs[0]["t_cw"]), 2),
                        "hex": q["hex"] if q.get("accept") else ""} for q in rs],
            "t_cw": rs[0]["t_cw"], "f_hz": rs[0]["f_hz"], "ebn0": round(eb, 2), "sync": round(sy, 2),
            "ok": dec.ok, "hex": dec.hex, "long": dec.long, "cost": round(dec.bch1_cost + dec.bch2_cost, 2),
            "flips": [dec.bch1_flips, dec.bch2_flips], "reasons": dec.reasons,
        })
    return out


def accept(res) -> str:
    """'single' when one decode is trustworthy alone, 'weak' when it needs confirmation by a
    second decode of the same ID, '' otherwise."""
    if not res.get("ok"):
        return ""
    f1, f2 = res["flips"]
    if res["ebn0"] >= 1.5 and res["sync"] >= 1.3 and f1 <= 3 and f2 <= 2:
        return "single"
    if res["ebn0"] >= 0.0 and res["sync"] >= 1.0:
        return "weak"
    return ""


QUIET = (1543.62e6, 1543.88e6)   # no signals seen here: noise-only calibration band
PRE, POST = 0.3, 0.9


def noise_test(path, fs, fc, n, L, rng):
    """Demodulate n random candidates in a band with no signal: every acceptance is false."""
    tot = path.stat().st_size / 4 / fs
    counts = {"single": 0, "weak": 0, "": 0}
    done = 0
    while done < n:
        t = rng.uniform(1.0, tot - 10.0)
        x = dsp.read_ci16(path, t - PRE, 8.0 + PRE + POST, fs)
        y, f1, mid = band_stream(x, fs, fc, *QUIET)
        for _ in range(min(50, n - done)):
            tf = t + rng.uniform(0.0, 8.0)
            f = rng.uniform(QUIET[0] + 2e3, QUIET[1] - 2e3)
            res = process_candidate(y, f1, mid, t - PRE, tf, f, L)
            if res is None or "ebn0" not in res:
                continue
            res.pop("_llr", None)
            done += 1
            a = accept(res)
            counts[a] += 1
            if a:
                print(json.dumps({"noise": True, "accept": a, **res}), flush=True)
    print(json.dumps({"noise_test": n, **{k or "rejected": v for k, v in counts.items()}}), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--dur", type=float, default=None)
    ap.add_argument("--thr", type=float, default=13.0)
    ap.add_argument("--L", type=int, default=10)
    ap.add_argument("--block", type=float, default=8.0)
    ap.add_argument("--noise-test", type=int, default=0)
    ap.add_argument("--no-combine", action="store_true")
    ap.add_argument("--bands", default="galileo,glonass", help="galileo,glonass,quiet")
    a = ap.parse_args()
    path = Path(a.file)
    fs, fc, t0 = meta_of(path)
    if a.noise_test:
        return noise_test(path, fs, fc, a.noise_test, a.L, np.random.default_rng(3))
    total = path.stat().st_size / 4 / fs
    end = total if a.dur is None else min(total, a.start + a.dur)
    t = a.start
    tic = time.time()
    n_c = n_ok = 0
    while t < end - 1.0:
        blk = min(a.block, end - t)
        t_x0 = max(0.0, t - PRE)
        x = dsp.read_ci16(path, t_x0, (t - t_x0) + blk + POST, fs)
        for band in a.bands.split(","):
            lo, hi = QUIET if band == "quiet" else BANDS[band]
            if lo < fc - fs / 2 or hi > fc + fs / 2:
                continue
            y, f1, mid = band_stream(x, fs, fc, lo, hi)
            recs = []
            for (tf, f_abs, snr, _) in merge(detect(y, f1, mid, lo, hi, t_x0, a.thr, band)):
                if not (t <= tf < t + blk):       # owned by another block
                    continue
                res = process_candidate(y, f1, mid, t_x0, tf, f_abs, a.L)
                if res is None:
                    continue
                n_c += 1
                recs.append({"t": round(tf, 3), "band": band, "det_snr": round(snr, 1), **res, "accept": accept(res)})
            if not a.no_combine:
                recs += [{**c, "accept": accept(c)} for c in combine_relays(recs, a.L)]
            for rec in recs:
                n_ok += bool(rec["accept"])
                if t0 is not None and "t_cw" in rec:
                    rec["utc"] = (t0 + dt.timedelta(seconds=rec["t_cw"])).isoformat(timespec="milliseconds")
                if "_llr" in rec:
                    rec["llr"] = [round(float(v), 1) for v in rec.pop("_llr")]
                print(json.dumps(rec), flush=True)
        t += blk
        el = time.time() - tic
        print(f"[{t:7.1f}/{end:.0f} s] candidates {n_c} accepted {n_ok}  speed {(t - a.start) / el:.1f}x",
              file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
