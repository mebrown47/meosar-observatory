"""Where is the receiving station? Fit its position from relayed SGB bursts.

Each copy of a burst came through a known satellite (satid.py --dump, run with any rough station
position: frequency alone attributes copies well enough), so for a trial station position every
copy implies the burst's transmit time and frequency; at the true position they agree. A coarse
grid then a fine one around the best cell; scored by the median absolute deviation of implied
transmit time and frequency within bursts. On 18 Nov 2026 data this fits to ~4 us and ~0.6 Hz,
against 0.43 ms and 79 Hz for a position ~260 km off.

usage: stationfit.py dump.jsonl frame_pos.json tle/galileo.tle tle/glo-ops.tle
       [--box LAT0 LAT1 LON0 LON1] (default: 15-60 N, 135-55 W)
"""
from __future__ import annotations

import argparse
import collections
import json
from datetime import datetime

import numpy as np
from skyfield.api import load, wgs84
from skyfield.framelib import itrs

from satid import C, F_NOM, T_XLATE, load_sats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("frame_pos")
    ap.add_argument("tles", nargs="+")
    ap.add_argument("--box", type=float, nargs=4, default=[15, 60, -135, -55])
    a = ap.parse_args()
    ts = load.timescale()
    sats = {s["label"]: s["sat"] for v in load_sats(None, a.tles, ts).values() for s in v}
    pos = json.load(open(a.frame_pos))
    bursts = collections.defaultdict(list)
    for r in map(json.loads, open(a.dump)):
        if r["gen"] == 2 and r["margin"] > 9:
            bursts[(r["utc"], r["hex"])].append(r)
    bursts = [(k, v) for k, v in bursts.items() if len(v) >= 2]
    P, V, BP, F, DT, T, idx = [], [], [], [], [], [], []
    for bi, ((utc, hx), v) in enumerate(bursts):
        t = ts.from_datetime(datetime.fromisoformat(utc))
        bpos = wgs84.latlon(pos[hx][0], pos[hx][1]).itrs_xyz.km
        for x in v:
            p, vel = sats[x["sat"]].at(t).frame_xyz_and_velocity(itrs)
            P.append(p.km); V.append(vel.km_per_s); BP.append(bpos)
            F.append(x["f_rx"]); DT.append(x["dt_ms"]); T.append(T_XLATE[x["band"]]); idx.append(bi)
    P, V, BP, F, DT, T, idx = map(np.array, (P, V, BP, F, DT, T, idx))
    nb = len(bursts)
    du = P - BP
    r_up = np.linalg.norm(du, axis=1)
    rr_up = np.sum(du * V, axis=1) / r_up

    def dem(x):
        return x - (np.bincount(idx, x, nb) / np.bincount(idx, None, nb))[idx]

    def spread(lat, lon):
        st = wgs84.latlon(lat, lon, elevation_m=100).itrs_xyz.km
        dd = P - st
        r_dn = np.linalg.norm(dd, axis=1)
        rr_dn = np.sum(dd * V, axis=1) / r_dn
        f0 = F - T + (T + F_NOM[2]) * rr_dn / C + F_NOM[2] * rr_up / C
        tx = DT - (r_up + r_dn) / C * 1e3
        mt, mf = np.median(np.abs(dem(tx))), np.median(np.abs(dem(f0)))
        return (mt / 0.15) ** 2 + (mf / 15) ** 2, mt, mf

    print(f"{nb} bursts, {len(F)} copies")
    lat0, lat1, lon0, lon1 = a.box
    coarse = min((spread(la, lo) + (la, lo) for la in np.arange(lat0, lat1, 1.0) for lo in np.arange(lon0, lon1, 1.0)))
    la, lo = coarse[3], coarse[4]
    fine = sorted(spread(x, y) + (x, y) for x in np.arange(la - 1.5, la + 1.5, 0.02) for y in np.arange(lo - 1.5, lo + 1.5, 0.02))
    for sc, mt, mf, x, y in fine[:5]:
        print(f"  {x:7.3f} {y:8.3f}: timing MAD {mt * 1e3:5.1f} us, frequency MAD {mf:5.1f} Hz")


if __name__ == "__main__":
    main()
