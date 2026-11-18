"""Satellite attribution: which MEOSAR satellite relayed each copy of a burst.

A copy relayed by satellite s arrives at
    f_rx = f0 + D_up(s) + T + D_dn(s)        t_rx = t_tx + (r_up(s) + r_dn(s)) / c + tau
f0 the beacon's frequency, D_up the Doppler beacon -> satellite at 406 MHz, T the transponder's
translation (Galileo 1138.05 MHz, GLONASS-K 1138.85 MHz), D_dn the Doppler satellite -> station
at L band, r the ranges, tau the transponder delay. Satellite orbits come from TLEs, so for a
beacon with a known position every term but f0 and t_tx is known: under the right assignment of
copies to satellites, every copy implies the same f0 and the same t_tx.

This script checks that on beacons whose frames carry a position (FGB location protocols, SGB
GNSS): for each burst it scores every assignment of its copies to distinct visible satellites by
the spread of implied f0 (sigma_f) and of implied t_tx (sigma_t), and reports the best, its
margin over the runner-up, and the implied f0 offsets (which also measure T).

Per-satellite translation and delay offsets are estimated from clearly attributed bursts and
removed (three passes). With the station's true position they come out at +-3 Hz and +-10 us, i.e.
nil: offsets of +-100-230 Hz and +-0.7 ms, varying with geometry, mean the station position is
wrong (stationfit.py finds it from the bursts). The common implied-f0 offset (-67 to -93 Hz from
406.05 MHz, drifting between days) is the receiver's own reference error.

TLEs:  curl -s "https://celestrak.org/NORAD/elements/gp.php?GROUP=galileo&FORMAT=tle" > tle/galileo.tle
       curl -s "https://celestrak.org/NORAD/elements/gp.php?GROUP=glo-ops&FORMAT=tle" > tle/glo-ops.tle
Needs skyfield (pip install skyfield).

usage: satid.py records.jsonl frame_pos.json tle/galileo.tle tle/glo-ops.tle --lat LAT --lon LON
       [--alt M] [--sar meosar_visible.py] [--sigma-f 15] [--sigma-t 0.15] [--max 40] [--dump FILE]
"""
from __future__ import annotations

import argparse
import collections
import itertools
import json
import math
import re
from datetime import datetime

import numpy as np
from skyfield.api import EarthSatellite, load, wgs84

C = 299792.458                       # km/s
T_XLATE = {"galileo": 1138.05e6, "glonass": 1138.85e6}
CONSTELLATION = {"galileo": "Galileo", "glonass": "GLONASS"}
F_NOM = {1: 406.04e6, 2: 406.05e6}   # for the uplink Doppler's scale (and SGB's nominal f0)


def load_sats(sar_file, tle_files, ts):
    """Every satellite in the TLE files, by constellation (galileo.tle -> Galileo, glo-*.tle ->
    GLONASS). dec406's SAR list only supplies Cospas-Sarsat IDs: it is incomplete (GSAT0222, for
    one, relays bursts but is not on it), so it does not restrict the candidates."""
    src = open(sar_file).read() if sar_file else ""
    cs_of = {int(b): int(a) for a, b, c, d in re.findall(r'SarSatellite\((\d+), (\d+), "(\w+)", "([^"]+)"', src)}
    out = collections.defaultdict(list)
    for f in tle_files:
        con = "GLONASS" if "glo" in f.lower() else "Galileo"
        L = [l.rstrip() for l in open(f) if l.strip()]
        for i in range(0, len(L) - 2, 3):
            n, l1, l2 = L[i].strip(), L[i + 1], L[i + 2]
            norad = int(l1[2:7])
            out[con].append({"cs": cs_of.get(norad), "label": n.split(" (")[0], "norad": norad,
                             "sat": EarthSatellite(l1, l2, n, ts)})
    return out


def geometry(sat, place, t):
    """(range km, range rate km/s, elevation deg) of sat seen from place at time t."""
    top = (sat - place).at(t)
    alt, _, rng, _, _, rr = top.frame_latlon_and_rates(place)
    return rng.km, rr.km_per_s, alt.degrees


def bursts_from(records, pos):
    """Bursts of beacons with a known position: [(utc, band, gen, hex, lat, lon, [(f_hz, dt_s)])]."""
    out = []
    fgb = collections.defaultdict(list)
    for r in records:
        k = r.get("kind")
        hx = r.get("hex")
        if not hx or hx not in pos or not r.get("utc"):
            continue
        lat, lon = pos[hx][0], pos[hx][1]
        if abs(lat - 1.0) < 1e-6 and abs(lon) < 1e-6:
            continue   # test transmitters' placeholder position (1 N 0 E), not where they are
        if k == "sgb" and r.get("accept") and r.get("relays"):
            rel = [(x["f_hz"], x["dt_ms"] / 1000) for x in r["relays"]]
            out.append((r["utc"], r["band"], 2, hx, lat, lon, rel))
        elif k is None and r.get("accept") and "combined" not in r:
            fgb[hx].append(r)
    # FGB: copies of one burst are separate records within 40 ms
    for hx, rs in fgb.items():
        rs.sort(key=lambda r: r["t_cw"])
        grp = []
        for r in rs + [None]:
            if r is not None and grp and r["t_cw"] - grp[0]["t_cw"] < 0.040 and r["band"] == grp[0]["band"]:
                if all(abs(r["f_hz"] - q["f_hz"]) > 20 for q in grp):
                    grp.append(r)
                continue
            if grp:
                g0 = grp[0]
                out.append((g0["utc"], g0["band"], 1, hx, pos[hx][0], pos[hx][1],
                            [(q["f_hz"], q["t_cw"] - g0["t_cw"]) for q in grp]))
            grp = [r] if r is not None else []
    return sorted(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("records")
    ap.add_argument("frame_pos")
    ap.add_argument("tles", nargs="+")
    ap.add_argument("--sar", default="/home/mike/dec406/tools/meosar_visible.py")
    ap.add_argument("--lat", type=float, required=True, help="station latitude (deg)")
    ap.add_argument("--lon", type=float, required=True, help="station longitude (deg)")
    ap.add_argument("--alt", type=float, default=100, help="station altitude (m)")
    ap.add_argument("--sigma-f", type=float, default=15.0, help="Hz")
    ap.add_argument("--sigma-t", type=float, default=0.15, help="ms")
    ap.add_argument("--min-relays", type=int, default=2)
    ap.add_argument("--max", type=int, default=40, help="bursts to print")
    ap.add_argument("--dump", help="write one JSON line per copy of each scored burst (residuals, geometry)")
    a = ap.parse_args()

    ts = load.timescale()
    sats = load_sats(a.sar, a.tles, ts)
    station = wgs84.latlon(a.lat, a.lon, elevation_m=a.alt)
    pos = json.load(open(a.frame_pos))
    recs = [json.loads(l) for l in open(a.records)]
    bursts = [b for b in bursts_from(recs, pos) if len(b[6]) >= a.min_relays]
    print(f"{len(bursts)} bursts with >= {a.min_relays} copies from beacons with known positions; "
          f"satellites: {', '.join(f'{k} {len(v)}' for k, v in sats.items())}")

    # Geometry once per burst: the candidate satellites with their Doppler and delay.
    prepared = []
    for utc, band, gen, hx, lat, lon, rel in bursts:
        t = ts.from_datetime(datetime.fromisoformat(utc))
        beacon = wgs84.latlon(lat, lon)
        cand = []
        for s in sats[CONSTELLATION[band]]:
            r_dn, rr_dn, el_st = geometry(s["sat"], station, t)
            r_up, rr_up, el_bc = geometry(s["sat"], beacon, t)
            if el_st < 5 or el_bc < 0:
                continue
            f_dn = T_XLATE[band] + F_NOM[gen]
            # implied f0 for a copy at f_rx: f_rx - T - D_dn - D_up
            dop = -f_dn * rr_dn / C - F_NOM[gen] * rr_up / C
            cand.append({"cs": s["cs"], "label": s["label"], "dop": dop, "delay": (r_up + r_dn) / C, "el": el_st,
                         "el_bc": el_bc, "r_up": r_up, "r_dn": r_dn, "rr_up": rr_up, "rr_dn": rr_dn})
        if len(cand) >= len(rel):
            prepared.append((utc, band, gen, hx, rel, cand))

    def score_all(fbias, tbias):
        """Best and runner-up assignment per burst, with per-satellite offsets removed.
        Copy 0's satellite anchors the beacon frequency; every other copy then only tries
        satellites whose implied f0 lies within WIN_HZ of it (all others would score a chi2 in the
        hundreds), so bursts with many copies and many visible satellites stay cheap."""
        WIN_HZ = 500.0
        res = []
        for utc, band, gen, hx, rel, cand in prepared:
            fb = [fbias.get((band, c["label"]), 0.0) for c in cand]
            tb = [tbias.get((band, c["label"]), 0.0) for c in cand]
            f0m = np.array([[rel[k][0] - T_XLATE[band] - c["dop"] - fb[j] for j, c in enumerate(cand)] for k in range(len(rel))])
            txm = np.array([[rel[k][1] * 1e3 - c["delay"] * 1e3 - tb[j] for j, c in enumerate(cand)] for k in range(len(rel))])
            top = []   # best two: (chi, chi_f, chi_t, perm, f0, tx)
            for j0 in range(len(cand)):
                ref = f0m[0, j0]
                opts = [[j0]] + [[j for j in range(len(cand)) if abs(f0m[k, j] - ref) < WIN_HZ] for k in range(1, len(rel))]
                for perm in itertools.product(*opts):
                    if len(set(perm)) < len(perm):
                        continue
                    f0 = f0m[range(len(rel)), perm]
                    tx = txm[range(len(rel)), perm]
                    chi_f = float(np.sum((f0 - f0.mean()) ** 2)) / a.sigma_f ** 2
                    chi_t = float(np.sum((tx - tx.mean()) ** 2)) / a.sigma_t ** 2
                    top.append((chi_f + chi_t, chi_f, chi_t, perm, f0, tx))
                    top.sort(key=lambda x: x[0])
                    del top[2:]
            if top:
                res.append((utc, band, gen, hx, rel, cand, top[0], top[1] if len(top) > 1 else None))
        return res

    def biases(res):
        """Per-satellite offsets (median residual from the burst mean) from clearly attributed bursts."""
        fr, tr = collections.defaultdict(list), collections.defaultdict(list)
        for utc, band, gen, hx, rel, cand, best, second in res:
            if second is not None and second[0] - best[0] < 25:
                continue
            for k in range(len(rel)):
                name = cand[best[3][k]]["label"]
                fr[(band, name)].append(best[4][k] - best[4].mean())
                tr[(band, name)].append(best[5][k] - best[5].mean())
        return ({k: float(np.median(v)) for k, v in fr.items() if len(v) >= 5},
                {k: float(np.median(v)) for k, v in tr.items() if len(v) >= 5}, fr, tr)

    def clear(res):
        return sum(1 for r in res if r[6][0] < 3 * len(r[4]) and (r[7] is None or r[7][0] - r[6][0] > 9))

    results = score_all({}, {})
    print(f"before calibration: {clear(results)} of {len(results)} bursts fit cleanly")
    fbias, tbias = {}, {}
    for it in range(3):   # offsets relative to the burst mean: iterate so they settle
        fb, tb, fr, tr = biases(score_all(fbias, tbias))
        for k, v in fb.items():
            fbias[k] = fbias.get(k, 0.0) + v
        for k, v in tb.items():
            tbias[k] = tbias.get(k, 0.0) + v
        results = score_all(fbias, tbias)
        print(f"calibration pass {it + 1}: {clear(results)} of {len(results)} bursts fit cleanly")
    print("per-satellite offsets (frequency Hz, delay ms; relative to the mean of the satellites seen):")
    for k in sorted(fbias):
        print(f"  {k[0]:7s} {k[1]:12s} {fbias[k]:+7.0f} Hz  {tbias.get(k, 0.0):+6.2f} ms  (n={len(fr.get(k, []))})")
    resid_f = collections.defaultdict(list)
    for utc, band, gen, hx, rel, cand, best, second in results:
        for k in range(len(rel)):
            resid_f[(band, gen)].append(best[4][k] - F_NOM[gen])

    if a.dump:
        with open(a.dump, "w") as out:
            for utc, band, gen, hx, rel, cand, best, second in results:
                margin = (second[0] - best[0]) if second else 1e9
                for k in range(len(rel)):
                    c = cand[best[3][k]]
                    key = (band, c["label"])
                    out.write(json.dumps({
                        "utc": utc, "band": band, "gen": gen, "hex": hx, "n": len(rel), "k": k, "sat": c["label"],
                        "f_rx": rel[k][0], "dt_ms": rel[k][1] * 1e3, "el": round(c["el"], 1), "el_bc": round(c["el_bc"], 1),
                        "r_up": round(c["r_up"], 1), "r_dn": round(c["r_dn"], 1), "rr_up": round(c["rr_up"], 4),
                        "rr_dn": round(c["rr_dn"], 4), "delay_ms": c["delay"] * 1e3,
                        "f_res": float(best[4][k] - best[4].mean()), "t_res": float(best[5][k] - best[5].mean()),
                        "fbias": fbias.get(key, 0.0), "tbias": tbias.get(key, 0.0),
                        "chi_f": float(best[1]), "chi_t": float(best[2]), "margin": float(margin)}) + "\n")

    for utc, band, gen, hx, rel, cand, best, second in results[:a.max]:
        names = [cand[i]["label"] for i in best[3]]
        margin = (second[0] - best[0]) if second else float("inf")
        print(f"{utc[11:19]} {band:7s} G{gen} {hx[:15]} {len(rel)} copies -> {', '.join(names)}"
              f"  chi2 f {best[1]:6.1f} t {best[2]:6.1f}  margin {margin:7.1f}"
              f"  f0-nom {np.mean(best[4]) - F_NOM[gen]:+8.0f} Hz")
    ok = [r for r in results if r[6][0] < 3 * len(r[4]) and (r[7] is None or r[7][0] - r[6][0] > 9)]
    print(f"\n{len(results)} bursts scored; {len(ok)} with a clear best assignment "
          f"(chi2 < 3 per copy, runner-up worse by > 9)")
    for key, v in sorted(resid_f.items()):
        v = np.array(v)
        print(f"implied f0 - nominal, {key[0]} G{key[1]}: median {np.median(v):+.0f} Hz, "
              f"spread (IQR) {np.subtract(*np.percentile(v, [75, 25])):.0f} Hz, n={len(v)}")


if __name__ == "__main__":
    main()
