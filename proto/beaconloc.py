"""Beacon geolocation from the relayed copies of its bursts (prototype).

A copy relayed by satellite s arrives at
    f_rx = f0 + D_up(s, P) + T + D_dn(s)        t_rx = t_tx + (r_up(s, P) + r_dn(s)) / c
with P the beacon's position (see observatory/satattr.py). With P unknown, the copies of a burst
still have to agree on one f0 and one t_tx under the right assignment to satellites, and only at
the right P. So: over a grid of trial positions, score each burst by its best assignment's chi2
(sigma 15 Hz, 0.15 ms; capped, so a bad burst cannot dominate), sum over the beacon's bursts, and
take the minimum; then refine on finer grids around the best coarse minima. A fixed beacon's
bursts over hours see the satellites move, which pins the position down.

Needs bursts with two copies or more (the observatory database keeps copies), the station's
position, and TLEs (python3-sgp4 + numpy).

usage: beaconloc.py DB TLE_DIR --station LAT,LON BID [BID ...] [--max-bursts 60] [--truth]
  --truth: for a beacon whose frames carry a position, also report the error (the position is
           never used in the search).
"""
from __future__ import annotations

import argparse
import itertools
import math
import os
import sqlite3
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "observatory"))
import satattr as SA   # noqa: E402

C = SA.C
CAP = 60.0          # per-burst chi2 cap
MIN_EL_BC = 0.0


def bursts_of(con, bid, max_bursts):
    gen = 2 if "-" in bid else 1
    table = "sgb_bursts" if gen == 2 else "fgb_bursts"
    rows = con.execute(
        f"SELECT b.id, b.t, b.band, f.lat, f.lon FROM {table} b JOIN frames f ON f.hex=b.hex WHERE b.bid=? "
        "AND (SELECT COUNT(*) FROM copies c WHERE c.gen=? AND c.burst=b.id) >= 2 ORDER BY b.t", (bid, gen)).fetchall()
    if len(rows) > max_bursts:      # spread over the whole span: geometry diversity
        idx = np.linspace(0, len(rows) - 1, max_bursts).round().astype(int)
        rows = [rows[i] for i in sorted(set(idx))]
    out = []
    for burst, t, band, lat, lon in rows:
        cp = con.execute("SELECT f_hz, tc FROM copies WHERE gen=? AND burst=? ORDER BY k", (gen, burst)).fetchall()
        t0 = min(c[1] for c in cp)
        out.append({"t": t, "band": band, "gen": gen, "f": np.array([c[0] for c in cp]),
                    "dt": np.array([(c[1] - t0) * 1e3 for c in cp]), "lat": lat, "lon": lon})
    return gen, out


def ecef_grid(lat, lon):
    """ECEF (G,3) and up vectors (G,3) of geodetic points (deg), on the ellipsoid."""
    a, f = 6378.137, 1 / 298.257223563
    e2 = f * (2 - f)
    la, lo = np.radians(lat), np.radians(lon)
    n = a / np.sqrt(1 - e2 * np.sin(la) ** 2)
    up = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)
    xyz = np.stack([n * up[:, 0], n * up[:, 1], n * (1 - e2) * up[:, 2]], -1)
    return xyz, up


def prepare(bursts, tles, station):
    """Per burst: satellites in view of the station (positions, velocities, downlink terms)."""
    out = []
    for b in bursts:
        sats = tles.at(b["t"], b["band"])
        R, V, dn = [], [], []
        f_dn = SA.T_XLATE[b["band"]] + SA.F_NOM[b["gen"]]
        for label, name, sat in sats:
            st = SA.sat_state(sat, b["t"])
            if st is None:
                continue
            r, v = st
            r_dn, rr_dn, el = SA.look(r, v, station)
            if el < SA.MIN_EL:
                continue
            R.append(r); V.append(v); dn.append((label, -f_dn * rr_dn / C, r_dn / C * 1e3))
        n = len(b["f"])
        if len(R) < n:
            continue
        perms = np.array(list(itertools.permutations(range(len(R)), n)), dtype=np.int16)
        out.append({**b, "R": np.array(R), "V": np.array(V), "labels": [d[0] for d in dn],
                    "dop_dn": np.array([d[1] for d in dn]), "del_dn": np.array([d[2] for d in dn]), "perms": perms})
    return out


def burst_chi(b, xyz, up, want_perm=False):
    """Best-assignment chi2 of burst b at each trial position (G,)."""
    d = b["R"][None, :, :] - xyz[:, None, :]                         # G,S,3
    r_up = np.linalg.norm(d, axis=-1)                                # G,S
    rr_up = np.einsum("gsk,sk->gs", d, b["V"]) / r_up
    el_ok = np.einsum("gsk,gk->gs", d, up) / r_up >= math.sin(math.radians(MIN_EL_BC))
    F = SA.F_NOM[b["gen"]]
    dop = b["dop_dn"][None, :] - F * rr_up / C                        # G,S
    dly = b["del_dn"][None, :] + r_up / C * 1e3
    T = SA.T_XLATE[b["band"]]
    n = len(b["f"])
    f0 = b["f"][None, None, :] - T - dop[:, :, None]                   # G,S,n (copy k via sat s)
    tx = b["dt"][None, None, :] - dly[:, :, None]
    bad = ~el_ok
    P = b["perms"]
    ks = np.arange(n)
    best = np.full(len(xyz), np.inf)
    bestp = np.zeros(len(xyz), dtype=np.int32)
    for i in range(0, len(P), 2000):
        p = P[i:i + 2000]                                            # Q,n
        fs = f0[:, p, ks]                                            # G,Q,n
        ts = tx[:, p, ks]
        chi = (fs.var(-1) * n / SA.SIG_F ** 2 + ts.var(-1) * n / SA.SIG_T ** 2)
        chi = np.where(bad[:, p].any(-1), np.inf, chi)
        j = chi.argmin(1)
        c = chi[np.arange(len(xyz)), j]
        upd = c < best
        best[upd] = c[upd]
        bestp[upd] = j[upd] + i
    return (best, bestp) if want_perm else best


def score(prep, lat, lon):
    xyz, up = ecef_grid(lat, lon)
    tot = np.zeros(len(lat))
    fit = np.zeros(len(lat), dtype=int)
    for b in prep:
        c = burst_chi(b, xyz, up)
        n = len(b["f"])
        tot += np.minimum(c, CAP)
        fit += c < 3 * n
    return tot, fit


def gc_km(lat1, lon1, lat2, lon2):
    r = math.radians
    h = math.sin(r(lat2 - lat1) / 2) ** 2 + math.cos(r(lat1)) * math.cos(r(lat2)) * math.sin(r(lon2 - lon1) / 2) ** 2
    return 12742 * math.asin(math.sqrt(min(1, h)))


def locate(prep, station_ll, coarse=2.0):
    # Coarse: everything within 80 deg of the station (a beacon further away shares no MEO
    # satellite with it above 5 deg).
    la = np.arange(-80, 80.01, coarse)
    lats, lons = [], []
    for x in la:
        w = coarse / max(0.2, math.cos(math.radians(x)))
        for y in np.arange(-180, 180, w):
            if gc_km(station_ll[0], station_ll[1], x, y) < 80 * 111.2:
                lats.append(x); lons.append(y)
    lats, lons = np.array(lats), np.array(lons)
    tot, fit = score(prep, lats, lons)
    order = np.argsort(tot)
    seeds = []
    for i in order:
        if all(gc_km(lats[i], lons[i], s[0], s[1]) > 1500 for s in seeds):
            seeds.append((lats[i], lons[i], tot[i]))
        if len(seeds) == 3:
            break
    results = []
    for la0, lo0, _ in seeds:
        c = (la0, lo0)
        r95 = None
        for span, step in ((3.0, 0.25), (0.4, 0.02), (0.04, 0.002)):
            g = np.arange(-span, span + step / 2, step)
            LA, LO = np.meshgrid(c[0] + g, c[1] + g / max(0.2, math.cos(math.radians(c[0]))), indexing="ij")
            t2, f2 = score(prep, LA.ravel(), LO.ravel())
            i = t2.argmin()
            c = (LA.ravel()[i], LO.ravel()[i])
            if step < 0.1:
                # 95% region (delta chi2 6 for two parameters): its farthest point from the best
                inside = t2 <= t2[i] + 6.0
                d = max(gc_km(c[0], c[1], x, y) for x, y in zip(LA.ravel()[inside], LO.ravel()[inside]))
                # the finest grid's estimate, unless the region runs off that grid
                if r95 is None or d < span * 111.2 * 0.8:
                    r95 = max(d, step * 111.2 / 2)
        results.append((t2[i], f2[i], c, r95))
    results.sort(key=lambda r: r[0])
    return results, (lats, lons, tot)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("db")
    ap.add_argument("tle_dir")
    ap.add_argument("bids", nargs="+")
    ap.add_argument("--station", required=True, help="LAT,LON")
    ap.add_argument("--max-bursts", type=int, default=60)
    ap.add_argument("--max-copies", type=int, default=5, help="drop bursts with more copies (assignment cost)")
    ap.add_argument("--truth", action="store_true")
    a = ap.parse_args()
    st_ll = tuple(float(x) for x in a.station.split(","))
    station = SA.ecef(st_ll[0], st_ll[1], 0.1)
    tles = SA.Tles(a.tle_dir)
    con = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    for bid in a.bids:
        t_start = time.time()
        gen, bursts = bursts_of(con, bid, a.max_bursts * 3)
        bursts = [b for b in bursts if len(b["f"]) <= a.max_copies]
        if len(bursts) > a.max_bursts:
            idx = np.linspace(0, len(bursts) - 1, a.max_bursts).round().astype(int)
            bursts = [bursts[i] for i in sorted(set(idx))]
        prep = prepare(bursts, tles, station)
        if len(prep) < 3:
            print(f"{bid}: only {len(prep)} usable bursts")
            continue
        span_h = (prep[-1]["t"] - prep[0]["t"]) / 3600
        res, _ = locate(prep, st_ll)
        (tot, fit, (la, lo), r95) = res[0]
        alt = "; ".join(f"{r[2][0]:.2f},{r[2][1]:.2f} score {r[0]:.0f} fit {r[1]}" for r in res[1:])
        line = (f"{bid}: {la:.3f}, {lo:.3f} +-{r95:.1f} km  score {tot:.0f}, {fit}/{len(prep)} bursts fit  "
                f"({len(prep)} bursts over {span_h:.1f} h, {time.time() - t_start:.0f} s)  next: {alt}")
        if a.truth:
            tr = [(b["lat"], b["lon"]) for b in prep if b["lat"] is not None]
            if tr:
                tl = tr[len(tr) // 2]
                line += f"  | truth {tl[0]:.3f},{tl[1]:.3f}  error {gc_km(la, lo, *tl):.1f} km"
        print(line, flush=True)


if __name__ == "__main__":
    main()
