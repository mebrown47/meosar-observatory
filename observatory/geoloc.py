#!/usr/bin/env python3
# Copyright (c) 2026 Mike Brown
# SPDX-License-Identifier: GPL-3.0-or-later
"""
MEOSAR Observatory: beacon positions estimated from the relayed copies of their bursts

A copy relayed by satellite s arrives at
    f_rx = f0 + D_up(s, P) + T + D_dn(s)        t_rx = t_tx + (r_up(s, P) + r_dn(s)) / c
with P the beacon's position (see satattr.py). With P unknown, the copies of a burst still have to
agree on one f0 and one t_tx under the right assignment to satellites, and only at the right P.
So: over a grid of trial positions, score each burst by its best assignment's chi2 (capped, so a
bad burst cannot dominate), sum over the beacon's bursts, and take the minimum; then refine on
finer grids around the best coarse minima. A fixed beacon's bursts over hours see the satellites
move, which pins the position down. A moving beacon fits nowhere and gets no estimate.

Checked against beacons whose frames carry a precise position: second generation within ~2 km,
first generation within 3-13 km. The 95% radius (from the chi2 surface) is never shown smaller
than 3 and 15 km. (Many first-generation frames carry only a coarse position, to 0.5 degrees, or
come from a beacon on a boat that moves: the estimate is then often the better one.)

Needs bursts with two copies or more (the database keeps copies), the station's position, TLEs
(rollup.py fetches them) and python3-sgp4 + python3-numpy.

usage:
  geoloc.py DB TLE_DIR --station LAT,LON --update [--budget S]
      the service: (re)estimate the beacons that are due, within S seconds, into the geoloc table
  geoloc.py DB TLE_DIR --station LAT,LON BID [BID ...] [--truth]
      print estimates for these beacons (nothing stored); --truth: for a beacon whose frames carry
      a position, also the error (the position is never used in the search)
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import satattr as SA   # noqa: E402

C = SA.C
CAP = 60.0          # per-burst chi2 cap
# Noise models: (sigma_f Hz, sigma_t ms by generation). LOOSE for the global search, where a sharp
# minimum could fall between grid points; TIGHT for the refinement, from the scatter of clearly
# attributed copies (median-based: ~1.3 Hz; ~9 us SGB, ~45 us FGB), with room for the heavy tails
# (the per-burst cap discards a burst with a bad copy).
LOOSE = (15.0, {1: 0.15, 2: 0.15})
TIGHT = (2.5, {1: 0.060, 2: 0.015})
MIN_EL_BC = 0.0


def bursts_of(con, bid, max_bursts, since=0.0):
    gen = 2 if "-" in bid else 1
    table = "sgb_bursts" if gen == 2 else "fgb_bursts"
    rows = con.execute(
        f"SELECT b.id, b.t, b.band, f.lat, f.lon FROM {table} b JOIN frames f ON f.hex=b.hex WHERE b.bid=? AND b.t>=? "
        "AND (SELECT COUNT(*) FROM copies c WHERE c.gen=? AND c.burst=b.id) >= 2 ORDER BY b.t", (bid, since, gen)).fetchall()
    if len(rows) > max_bursts:      # spread over the whole span: geometry diversity
        idx = np.linspace(0, len(rows) - 1, max_bursts).round().astype(int)
        rows = [rows[i] for i in sorted(set(idx))]
    out = []
    for burst, t, band, lat, lon in rows:
        cp = con.execute("SELECT f_hz, tc, ebn0, band FROM copies WHERE gen=? AND burst=? ORDER BY k", (gen, burst)).fetchall()
        t0 = min(c[1] for c in cp)
        # a first-generation burst's copies can be in both bands (Galileo and GLONASS relays)
        out.append({"t": t, "band": band, "gen": gen, "f": np.array([c[0] for c in cp]),
                    "bands": [SA.band_of(c[3], c[0]) for c in cp],
                    "eb": np.array([c[2] if c[2] is not None else -99.0 for c in cp]),
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


def make_perms(copy_bands, sat_bands):
    """Every assignment of copies to distinct satellites of their own band: (P, n) satellite
    indices, built band by band (far fewer than all permutations of all satellites)."""
    n = len(copy_bands)
    parts = []
    for band in sorted(set(copy_bands)):
        ks = [k for k in range(n) if copy_bands[k] == band]
        ss = [j for j, b in enumerate(sat_bands) if b == band]
        if len(ss) < len(ks):
            return None
        parts.append((ks, list(itertools.permutations(ss, len(ks)))))
    out = []
    for combo in itertools.product(*[p[1] for p in parts]):
        row = [0] * n
        for (ks, _), sel in zip(parts, combo):
            for k, j in zip(ks, sel):
                row[k] = j
        out.append(row)
    return np.array(out, dtype=np.int16)


def prepare(bursts, tles, station):
    """Per burst: satellites of its copies' bands in view of the station (positions, velocities,
    downlink terms, translation)."""
    out = []
    for b in bursts:
        R, V, dn = [], [], []
        for band in sorted(set(b["bands"])):
            f_dn = SA.T_XLATE[band] + SA.F_NOM[b["gen"]]
            for label, name, sat in tles.at(b["t"], band):
                st = SA.sat_state(sat, b["t"])
                if st is None:
                    continue
                r, v = st
                r_dn, rr_dn, el = SA.look(r, v, station)
                if el < SA.MIN_EL:
                    continue
                R.append(r); V.append(v); dn.append((label, -f_dn * rr_dn / C, r_dn / C * 1e3, band))
        perms = make_perms(b["bands"], [d[3] for d in dn])
        if perms is None or not len(perms):
            continue
        out.append({**b, "R": np.array(R), "V": np.array(V), "labels": [d[0] for d in dn],
                    "dop_dn": np.array([d[1] for d in dn]), "del_dn": np.array([d[2] for d in dn]),
                    "T": np.array([SA.T_XLATE[d[3]] for d in dn]), "sat_bands": [d[3] for d in dn], "perms": perms})
    return out


def burst_chi(b, xyz, up, want_perm=False, model=LOOSE):
    """Best-assignment chi2 of burst b at each trial position (G,)."""
    d = b["R"][None, :, :] - xyz[:, None, :]                         # G,S,3
    r_up = np.linalg.norm(d, axis=-1)                                # G,S
    rr_up = np.einsum("gsk,sk->gs", d, b["V"]) / r_up
    el_ok = np.einsum("gsk,gk->gs", d, up) / r_up >= math.sin(math.radians(MIN_EL_BC))
    F = SA.F_NOM[b["gen"]]
    dop = b["dop_dn"][None, :] - F * rr_up / C                        # G,S
    dly = b["del_dn"][None, :] + r_up / C * 1e3
    n = len(b["f"])
    f0 = b["f"][None, None, :] - b["T"][None, :, None] - dop[:, :, None]   # G,S,n (copy k via sat s)
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
        chi = (fs.var(-1) * n / model[0] ** 2 + ts.var(-1) * n / model[1][b["gen"]] ** 2)
        chi = np.where(bad[:, p].any(-1), np.inf, chi)
        j = chi.argmin(1)
        c = chi[np.arange(len(xyz)), j]
        upd = c < best
        best[upd] = c[upd]
        bestp[upd] = j[upd] + i
    return (best, bestp) if want_perm else best


def score(prep, lat, lon, model=LOOSE):
    xyz, up = ecef_grid(lat, lon)
    tot = np.zeros(len(lat))
    fit = np.zeros(len(lat), dtype=int)
    for b in prep:
        c = burst_chi(b, xyz, up, model=model)
        n = len(b["f"])
        tot += np.minimum(c, CAP)
        fit += c < 3 * n
    return tot, fit


def gc_km(lat1, lon1, lat2, lon2):
    r = math.radians
    h = math.sin(r(lat2 - lat1) / 2) ** 2 + math.cos(r(lat1)) * math.cos(r(lat2)) * math.sin(r(lon2 - lon1) / 2) ** 2
    return 12742 * math.asin(math.sqrt(min(1, h)))


def strongest(prep, m=3):
    """The bursts cut down to their m strongest copies: any subset of a burst's copies must agree
    too, and m copies cost far fewer assignments (the coarse global search)."""
    out = []
    for b in prep:
        n = len(b["f"])
        if n <= m:
            out.append(b)
            continue
        keep = np.sort(np.argsort(-b["eb"])[:m])
        bands = [b["bands"][k] for k in keep]
        out.append({**b, "f": b["f"][keep], "dt": b["dt"][keep], "eb": b["eb"][keep], "bands": bands,
                    "perms": make_perms(bands, b["sat_bands"])})
    return out


def locate(prep, station_ll, coarse=2.0):
    # Coarse: everything within 150 deg of the station. A MEO satellite (23,200 km up) is above
    # 5 deg at the station within ~73 deg of the sub-satellite point, and above the horizon at the
    # beacon within ~77 deg: the two can share one up to ~150 deg apart.
    la = np.arange(-88, 88.01, coarse)
    lats, lons = [], []
    for x in la:
        w = coarse / max(0.2, math.cos(math.radians(x)))
        for y in np.arange(-180, 180, w):
            if gc_km(station_ll[0], station_ll[1], x, y) < 150 * 111.2:
                lats.append(x); lons.append(y)
    lats, lons = np.array(lats), np.array(lons)
    tot, fit = score(strongest(prep), lats, lons)
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
        grids = []
        for span, step in ((3.0, 0.25), (0.4, 0.02), (0.04, 0.002)):
            g = np.arange(-span, span + step / 2, step)
            LA, LO = np.meshgrid(c[0] + g, c[1] + g / max(0.2, math.cos(math.radians(c[0]))), indexing="ij")
            t2, f2 = score(prep, LA.ravel(), LO.ravel(), LOOSE if step > 0.1 else TIGHT)
            i = t2.argmin()
            c = (LA.ravel()[i], LO.ravel()[i])
            grids.append((span, step, LA.ravel(), LO.ravel(), t2))
        best, fit = t2[i], f2[i]
        # The noise model is a guess: scale chi2 by the bursts that fit at the best
        # position (reduced chi2; each burst of n copies has 2(n-1) degrees of freedom).
        xyz, up = ecef_grid(np.array([c[0]]), np.array([c[1]]))
        chis = [(burst_chi(b, xyz, up, model=TIGHT)[0], len(b["f"])) for b in prep]
        ok = [(x, n) for x, n in chis if x < 3 * n]
        scale = sum(x for x, _ in ok) / max(1, sum(2 * (n - 1) for _, n in ok)) if ok else 1.0
        scale = max(scale, 0.05)
        # 95% region: delta chi2 6 (two parameters) after scaling; its farthest point from the best,
        # from the finest grid it fits inside (or the coarse one's extent, as a lower bound).
        r95, edge = None, False
        for span, step, la_, lo_, t_ in reversed(grids[1:]):
            inside = t_ <= t_.min() + 6.0 * scale
            d = max(gc_km(c[0], c[1], x, y) for x, y in zip(la_[inside], lo_[inside]))
            r95, edge = max(d, step * 111.2 / 2), d >= span * 111.2 * 0.8
            if not edge:
                break
        results.append((best, fit, c, r95, edge, scale))
    results.sort(key=lambda r: r[0])
    return results, (lats, lons, tot)


# ---- the service: estimates kept in the database ------------------------------------------------

WINDOW_S = 72 * 3600     # a beacon's bursts from its last 72 h: plenty of geometry, and recent
MIN_BURSTS = 10          # bursts with two copies or more
MIN_SPAN_S = 3600        # the satellites have to move
REDO_S = 6 * 3600        # re-estimate a beacon still heard after this long
FLOOR_KM = {1: 15.0, 2: 3.0}   # validation errors (see the top): the radius is never smaller


def init_table(con):
    con.execute(
        """CREATE TABLE IF NOT EXISTS geoloc (
            bid TEXT PRIMARY KEY, gen INTEGER, quality TEXT, lat REAL, lon REAL, r95 REAL, r95_fit REAL,
            bursts INTEGER, fit INTEGER, span_h REAL, t_first REAL, t_last REAL, score REAL, alt_score REAL,
            alt_lat REAL, alt_lon REAL, scale REAL, truth_km REAL, computed REAL, secs REAL)""")


def quality(n, fit, best, alt, edge):
    """good: most bursts fit and nowhere else comes close; rough: a clear minimum, but many bursts
    disagree (a beacon that moved, or a noisy one); none: no position."""
    sep = (alt - best) / max(1, n)          # per burst; a burst that fits nowhere scores CAP
    if fit >= 0.8 * n and sep >= 20 and not edge:
        return "good"
    if fit >= 0.5 * n and sep >= 10:
        return "rough"
    return "none"


def locate_bid(con, tles, station, st_ll, bid, max_bursts=30, max_copies=5, since=0.0):
    """One beacon's estimate as a dict, or None (too few usable bursts)."""
    t_start = time.time()
    gen, bursts = bursts_of(con, bid, max_bursts * 3, since)
    bursts = [b for b in bursts if len(b["f"]) <= max_copies]
    if len(bursts) > max_bursts:
        idx = np.linspace(0, len(bursts) - 1, max_bursts).round().astype(int)
        bursts = [bursts[i] for i in sorted(set(idx))]
    prep = prepare(bursts, tles, station)
    if len(prep) < 3:
        return None
    res, _ = locate(prep, st_ll)
    best, fit, (la, lo), r95, edge, scale = res[0]
    n = len(prep)
    alt = res[1] if len(res) > 1 else (n * CAP, 0, (None, None))
    tr = sorted((b["lat"], b["lon"]) for b in prep if b["lat"] is not None)
    truth = None
    if tr:
        tl = tr[len(tr) // 2]
        truth = gc_km(la, lo, tl[0], tl[1])
    return {"bid": bid, "gen": gen, "quality": quality(n, fit, best, alt[0], edge), "lat": float(la), "lon": float(lo),
            "r95": max(float(r95), FLOOR_KM[gen]), "r95_fit": float(r95), "edge": edge, "bursts": n, "fit": int(fit),
            "span_h": (prep[-1]["t"] - prep[0]["t"]) / 3600, "t_first": prep[0]["t"], "t_last": prep[-1]["t"],
            "score": float(best), "alt_score": float(alt[0]), "alt_lat": alt[2][0], "alt_lon": alt[2][1],
            "scale": float(scale), "truth_km": truth, "computed": time.time(), "secs": time.time() - t_start,
            "truth": tr[len(tr) // 2] if tr else None}


def due(con, now):
    """Beacons to (re)estimate, most wanted first: never estimated, then without a position of
    their own, then the longest since."""
    newest = con.execute("SELECT MAX(t) FROM (SELECT MAX(t) t FROM fgb_bursts UNION ALL SELECT MAX(t) FROM sgb_bursts)").fetchone()[0]
    if newest is None:
        return [], 0.0
    since = newest - WINDOW_S
    have = {r[0]: r[1:] for r in con.execute("SELECT bid, t_last, computed FROM geoloc")}
    out = []
    for gen, table in ((1, "fgb_bursts"), (2, "sgb_bursts")):
        for bid, n, t0, t1, haspos in con.execute(
                f"SELECT b.bid, COUNT(*), MIN(b.t), MAX(b.t), MAX(f.lat IS NOT NULL) FROM {table} b JOIN frames f ON f.hex=b.hex "
                "WHERE b.t>=? AND (SELECT COUNT(*) FROM copies c WHERE c.gen=? AND c.burst=b.id) >= 2 GROUP BY b.bid", (since, gen)):
            if n < MIN_BURSTS or t1 - t0 < MIN_SPAN_S:
                continue
            h = have.get(bid)
            if h and not (t1 > h[0] + 60 and now - h[1] > REDO_S):
                continue
            out.append(((h is not None, bool(haspos), h[1] if h else 0), bid))
    out.sort()
    return [bid for _, bid in out], since


def update(db, tle_dir, st_ll, budget, verbose=False):
    con = sqlite3.connect(db, timeout=60)
    init_table(con)
    con.commit()
    tles = SA.Tles(tle_dir)
    station = SA.ecef(st_ll[0], st_ll[1], 0.1)
    t_end = time.time() + budget
    bids, since = due(con, time.time())
    done = {}
    for bid in bids:
        if time.time() > t_end:
            break
        r = locate_bid(con, tles, station, st_ll, bid, since=since)
        if r is None:
            # too few usable bursts (no satellites, too many copies): look again after REDO_S
            r = {"bid": bid, "gen": 2 if "-" in bid else 1, "quality": "none", "computed": time.time(),
                 "t_last": con.execute("SELECT MAX(t) FROM (SELECT MAX(t) t FROM fgb_bursts WHERE bid=? "
                                       "UNION ALL SELECT MAX(t) FROM sgb_bursts WHERE bid=?)", (bid, bid)).fetchone()[0]}
        cols = ["bid", "gen", "quality", "lat", "lon", "r95", "r95_fit", "bursts", "fit", "span_h", "t_first", "t_last",
                "score", "alt_score", "alt_lat", "alt_lon", "scale", "truth_km", "computed", "secs"]
        con.execute(f"INSERT OR REPLACE INTO geoloc({','.join(cols)}) VALUES({','.join('?' * len(cols))})",
                    [r.get(c) for c in cols])
        con.commit()
        done[r["quality"]] = done.get(r["quality"], 0) + 1
        if verbose:
            print(f"geoloc: {bid}: {r['quality']}" + (f" {r['lat']:.3f},{r['lon']:.3f} +-{r['r95']:.1f} km, "
                  f"{r['fit']}/{r['bursts']} fit, {r['secs']:.0f} s" if "lat" in r else ""), flush=True)
    con.close()
    left = len(bids) - sum(done.values())
    print(f"geoloc: {sum(done.values())} beacon(s) estimated" + "".join(f", {v} {k}" for k, v in sorted(done.items())) +
          (f"; {left} more due next time" if left else ""))


def main():
    ap = argparse.ArgumentParser(description="MEOSAR beacon positions from relayed copies")
    ap.add_argument("db")
    ap.add_argument("tle_dir")
    ap.add_argument("bids", nargs="*")
    ap.add_argument("--station", required=True, help="the station's position, LAT,LON")
    ap.add_argument("--update", action="store_true", help="estimate the beacons that are due into the database")
    ap.add_argument("--budget", type=float, default=900, help="--update: stop starting beacons after this many seconds")
    ap.add_argument("--max-bursts", type=int, default=30)
    ap.add_argument("--max-copies", type=int, default=5, help="drop bursts with more copies (assignment cost)")
    ap.add_argument("--truth", action="store_true")
    ap.add_argument("--verbose", "-v", action="store_true")
    a = ap.parse_args()
    try:
        st_ll = tuple(float(x) for x in a.station.split(","))
        assert len(st_ll) == 2
    except (ValueError, AssertionError):
        print(f"geoloc: --station {a.station!r} is not LAT,LON: nothing to do", file=sys.stderr)
        return 0
    if a.update:
        update(a.db, a.tle_dir, st_ll, a.budget, a.verbose)
        return 0
    station = SA.ecef(st_ll[0], st_ll[1], 0.1)
    tles = SA.Tles(a.tle_dir)
    con = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    for bid in a.bids:
        r = locate_bid(con, tles, station, st_ll, bid, a.max_bursts, a.max_copies)
        if r is None:
            print(f"{bid}: too few usable bursts")
            continue
        line = (f"{bid}: {r['quality']} {r['lat']:.3f}, {r['lon']:.3f} +-{'>' if r['edge'] else ''}{r['r95_fit']:.1f} km "
                f"(shown +-{r['r95']:.0f}; chi2 scale {r['scale']:.2f})  score {r['score']:.0f}, {r['fit']}/{r['bursts']} bursts fit  "
                f"({r['bursts']} bursts over {r['span_h']:.1f} h, {r['secs']:.0f} s)  next: "
                + (f"{r['alt_lat']:.2f},{r['alt_lon']:.2f} score {r['alt_score']:.0f}" if r["alt_lat"] is not None else "none"))
        if a.truth and r["truth"]:
            line += f"  | truth {r['truth'][0]:.3f},{r['truth'][1]:.3f}  error {r['truth_km']:.1f} km"
        print(line, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
