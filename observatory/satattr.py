# Copyright (c) 2026 Mike Brown
# SPDX-License-Identifier: GPL-3.0-or-later
"""
MEOSAR Observatory: satellite attribution (which satellite relayed each copy of a burst)

The rollup keeps every relayed copy of a burst (its carrier frequency and its arrival time
relative to the burst's other copies). A copy relayed by satellite s arrives at

    f_rx = f0 + D_up(s) + T + D_dn(s)        t_rx = t_tx + (r_up(s) + r_dn(s)) / c

f0 the beacon's frequency, D_up the Doppler beacon -> satellite at 406 MHz, T the transponder's
translation (Galileo 1138.05 MHz, GLONASS-K 1138.85 MHz), D_dn the Doppler satellite -> station at
L band, r the ranges. For a beacon whose frame carries its position, orbits (TLEs, propagated with
SGP4) give every term but f0 and t_tx; under the right assignment of copies to satellites, every
copy implies the same f0 and the same t_tx. Each burst's assignments are scored by the spread of
implied f0 and t_tx (sigmas from the measured scatter: SIG); one is "clear" when it fits (chi2 < 3 per copy)
and the runner-up is worse by more than 9. The station's position must be right (to ~1 km):
proto/satid.py and proto/stationfit.py are the offline versions, and found that no per-satellite
offsets are needed. The implied f0's common offset from nominal is the receiver's reference error.

A first-generation burst often reaches the station through both constellations: its copies then
lie in both bands (Galileo 1544.1 MHz, GLONASS 1544.9 MHz) and are scored together, each copy only
against satellites of its own band and with that band's translation, all sharing one f0 and t_tx.

A burst with a single copy has no spread to score; it is attributed by frequency alone when the
same beacon's f0 is known from a clearly attributed burst within an hour (beacon oscillators drift
by a few Hz in that time; candidate satellites' Doppler usually differs by hundreds).

TLEs (CelesTrak's galileo and glo-ops groups) are fetched once a day into TLE_DIR/YYYY-MM-DD/; each
burst uses the snapshot nearest its time, and none older than 30 days. Needs python3-sgp4; without
it (or without the station's position) the rollup goes on without attribution.
"""

import bisect
import glob
import math
import os
import statistics
import sys
import time
import urllib.request
from datetime import datetime, timezone

try:
    from sgp4.api import Satrec, jday
except ImportError:                     # attribution is optional
    Satrec = None

C = 299792.458                          # km/s
OMEGA = 7.292115e-5                     # Earth's rotation, rad/s
# Translation (uplink 406 MHz -> downlink 1544 MHz). GLONASS-K's is 4 Hz above its nominal 1138.85 MHz
# relative to Galileo's: copies of one burst relayed through both imply a beacon frequency 3.98 Hz
# higher through GLONASS (326 bursts, 3 days; +2.0 to +4.8 Hz by satellite).
T_XLATE = {"galileo": 1138.05e6, "glonass": 1138.85e6 + 4.0}
F_NOM = {1: 406.04e6, 2: 406.05e6}      # for the uplink Doppler's scale
GROUPS = {"galileo": "galileo", "glonass": "glo-ops"}   # band -> CelesTrak group (= file name)
# Noise model by generation: (sigma of a copy's implied f0 in Hz, of its implied t_tx in ms), from
# the measured robust scatter of clearly attributed copies (FGB 1.3 Hz, 36 us; SGB 0.7 Hz, 4 us;
# hardly dependent on Eb/N0), widened for the heavy tails: FGB frequency most (p99 9 Hz), so a
# two-copy FGB burst whose copies differ by 10 Hz still fits. Against 15 Hz / 0.15 ms for both
# (the first model), on 3 days of bursts: ambiguous 471 -> 93, clear +237, and no burst clear
# under both changed satellites.
SIG = {1: (3.5, 0.070), 2: (1.5, 0.008)}
# A single copy against the beacon's f0 from bursts within REF_S: its own scatter plus the
# beacon's frequency wander within the hour (FGB 1.0 Hz, SGB 0.4 Hz; p99 9 and 8 Hz). Candidate
# satellites usually imply f0s hundreds of Hz apart, so this can afford the tails.
SIG_REF = {1: 5.0, 2: 3.0}
MODEL = "2"     # change with the model: every burst is attributed again under the new one
WIN_HZ = 500.0                          # copies' implied f0 must lie this close to the anchor's
MIN_EL = 5.0                            # station elevation, deg
SETTLE_S = 60                           # attribute a burst once every copy of it is in
REF_S = 3600                            # single copies: f0 from the beacon's bursts this close
TLE_MAX_AGE_D = 30
FETCH_EVERY_S = 86400
RETRY_S = 3600
CELESTRAK = "https://celestrak.org/NORAD/elements/gp.php?GROUP={}&FORMAT=tle"


def init_tables(con):
    con.executescript(
        """
        -- Each relayed copy of a burst. gen 1: burst = fgb_bursts.id, tc = the copy's t_cw (s);
        -- gen 2: burst = sgb_bursts.id, tc = its offset from the burst's first copy (s).
        -- sat ... salt: the attribution (satellite, its elevation at the station, residuals,
        -- sub-satellite point and altitude in km), when the burst is attributed.
        CREATE TABLE IF NOT EXISTS copies (
            gen INTEGER, burst INTEGER, k INTEGER, tc REAL, f_hz REAL, ebn0 REAL,
            sat TEXT, el REAL, f_res REAL, t_res REAL, slat REAL, slon REAL, salt REAL,
            band TEXT, PRIMARY KEY (gen, burst, k));
        -- One row per burst looked at. state: clear, ambiguous (the runner-up fits nearly as well),
        -- nofit, noref (a single copy, and no f0 for the beacon), nosat (too few satellites in view).
        -- f0: the beacon frequency the copies imply (Hz); sats: the satellites, comma-separated.
        CREATE TABLE IF NOT EXISTS attr (
            gen INTEGER, burst INTEGER, bid TEXT, t REAL, band TEXT, n INTEGER, state TEXT,
            chi_f REAL, chi_t REAL, margin REAL, f0 REAL, sats TEXT,
            PRIMARY KEY (gen, burst));
        CREATE INDEX IF NOT EXISTS ix_attr_t ON attr(t);
        CREATE INDEX IF NOT EXISTS ix_attr_bid ON attr(bid, t);
        -- Where the satellites that have relayed bursts are now (refreshed every run).
        CREATE TABLE IF NOT EXISTS sat_now (
            sat TEXT PRIMARY KEY, name TEXT, constellation TEXT, t REAL, lat REAL, lon REAL, alt REAL,
            el REAL, az REAL);
        """
    )
    # Added after schema 2: each copy's band (a first-generation burst's copies can be in both).
    if "band" not in [c[1] for c in con.execute("PRAGMA table_info(copies)")]:
        con.execute("ALTER TABLE copies ADD COLUMN band TEXT")


def band_of(band, f_hz):
    """A copy's band: as recorded, or from its carrier (Galileo 1544.1 MHz, GLONASS 1544.9 MHz)."""
    return band or ("galileo" if f_hz < 1544.5e6 else "glonass")


# ---------- geometry ----------

def ecef(lat, lon, h_km=0.0):
    a, f = 6378.137, 1 / 298.257223563
    e2 = f * (2 - f)
    la, lo = math.radians(lat), math.radians(lon)
    n = a / math.sqrt(1 - e2 * math.sin(la) ** 2)
    up = (math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la))
    return ((n + h_km) * up[0], (n + h_km) * up[1], (n * (1 - e2) + h_km) * up[2]), up, (la, lo)


def gmst(jd, fr):
    """Greenwich mean sidereal time (IAU 1982, as SGP4's TEME frame uses), radians. UTC stands in
    for UT1: under a second, ~1 km at the satellites, a few microseconds of delay."""
    t = (jd - 2451545.0 + fr) / 36525.0
    s = 67310.54841 + (876600.0 * 3600 + 8640184.812866) * t + 0.093104 * t * t - 6.2e-6 * t ** 3
    return math.radians((s % 86400.0) / 240.0)


def sat_state(sat, t):
    """Satellite position (km) and velocity (km/s) in the Earth-fixed frame at epoch time t."""
    d = datetime.fromtimestamp(t, timezone.utc)
    jd, fr = jday(d.year, d.month, d.day, d.hour, d.minute, d.second + d.microsecond * 1e-6)
    e, r, v = sat.sgp4(jd, fr)
    if e:
        return None
    g = gmst(jd, fr)
    c, s = math.cos(g), math.sin(g)
    x, y = c * r[0] + s * r[1], -s * r[0] + c * r[1]
    vx, vy = c * v[0] + s * v[1] + OMEGA * y, -s * v[0] + c * v[1] - OMEGA * x
    return (x, y, r[2]), (vx, vy, v[2])


def look(r, v, obs):
    """(range km, range rate km/s, elevation deg) of a satellite at r, v from obs = ecef(...)."""
    o, up = obs[0], obs[1]
    d = (r[0] - o[0], r[1] - o[1], r[2] - o[2])
    rng = math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
    rr = (d[0] * v[0] + d[1] * v[1] + d[2] * v[2]) / rng
    el = math.degrees(math.asin((d[0] * up[0] + d[1] * up[1] + d[2] * up[2]) / rng))
    return rng, rr, el


def azimuth(r, obs):
    o, _, (la, lo) = obs
    d = (r[0] - o[0], r[1] - o[1], r[2] - o[2])
    e = -math.sin(lo) * d[0] + math.cos(lo) * d[1]
    n = -math.sin(la) * math.cos(lo) * d[0] - math.sin(la) * math.sin(lo) * d[1] + math.cos(la) * d[2]
    return math.degrees(math.atan2(e, n)) % 360


def subpoint(r):
    """Geodetic latitude, longitude and altitude (km) under a satellite at r (Earth-fixed)."""
    a, f = 6378.137, 1 / 298.257223563
    e2 = f * (2 - f)
    p = math.hypot(r[0], r[1])
    lat = math.atan2(r[2], p * (1 - e2))
    for _ in range(5):
        n = a / math.sqrt(1 - e2 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - n
        lat = math.atan2(r[2], p * (1 - e2 * n / (n + h)))
    return math.degrees(lat), math.degrees(math.atan2(r[1], r[0])), h


# ---------- TLEs ----------

class Tles:
    """TLE snapshots in tle_dir/YYYY-MM-DD/<group>.tle: per band, the satellites of the snapshot
    nearest a given time."""

    def __init__(self, tle_dir):
        self.snaps = []      # [(epoch of the date, {band: [(label, name, satrec)]})]
        for d in sorted(glob.glob(os.path.join(tle_dir, "????-??-??"))):
            try:
                t = datetime.strptime(os.path.basename(d), "%Y-%m-%d").replace(tzinfo=timezone.utc).timestamp()
            except ValueError:
                continue
            sats = {band: self.read(os.path.join(d, grp + ".tle")) for band, grp in GROUPS.items()}
            if any(sats.values()):
                self.snaps.append((t, sats))
        self.times = [s[0] for s in self.snaps]

    @staticmethod
    def read(path):
        try:
            lines = [l.rstrip() for l in open(path) if l.strip()]
        except OSError:
            return []
        out = []
        for i in range(0, len(lines) - 2, 3):
            name, l1, l2 = lines[i].strip(), lines[i + 1], lines[i + 2]
            if not (l1.startswith("1 ") and l2.startswith("2 ")):
                continue
            try:
                out.append((name.split(" (")[0], name, Satrec.twoline2rv(l1, l2)))
            except Exception:  # noqa: one bad entry
                continue
        return out

    def at(self, t, band):
        if not self.snaps:
            return []
        i = bisect.bisect_left(self.times, t)
        j = min((k for k in (i - 1, i) if 0 <= k < len(self.snaps)), key=lambda k: abs(self.times[k] - t))
        if abs(self.times[j] - t) > TLE_MAX_AGE_D * 86400:
            return []
        return self.snaps[j][1].get(band, [])

    def latest(self):
        return self.snaps[-1][1] if self.snaps else {}


def fetch_tles(con, tle_dir, verbose):
    """Today's snapshot, once a day; a failure is retried an hour later. Returns True if new."""
    now = time.time()
    snaps = sorted(glob.glob(os.path.join(tle_dir, "????-??-??")))
    newest = os.path.getmtime(snaps[-1]) if snaps else 0
    last_try = float((con.execute("SELECT value FROM meta WHERE key='tle_try'").fetchone() or ["0"])[0])
    if now - newest < FETCH_EVERY_S or now - last_try < RETRY_S:
        return False
    con.execute("INSERT INTO meta(key,value) VALUES('tle_try',?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
                (str(now),))
    con.commit()
    day = os.path.join(tle_dir, datetime.now(timezone.utc).strftime("%Y-%m-%d"))
    got = {}
    for grp in GROUPS.values():
        try:
            req = urllib.request.Request(CELESTRAK.format(grp), headers={"User-Agent": "meosar-observatory"})
            with urllib.request.urlopen(req, timeout=20) as r:
                body = r.read().decode("ascii", "replace")
        except Exception as e:  # noqa: offline is fine, the old snapshot stays in use
            print(f"satattr: could not fetch TLEs ({grp}): {e}", file=sys.stderr)
            return False
        if body.count("\n1 ") < 3:
            print(f"satattr: CelesTrak's {grp} answer has no TLEs: {body[:80]!r}", file=sys.stderr)
            return False
        got[grp] = body
    os.makedirs(day, exist_ok=True)
    for grp, body in got.items():
        tmp = os.path.join(day, grp + ".tle.tmp")
        with open(tmp, "w") as f:
            f.write(body)
        os.replace(tmp, os.path.join(day, grp + ".tle"))
    os.utime(day)
    if verbose:
        print(f"  TLEs fetched into {day}")
    return True


# ---------- attribution ----------

def candidates(sats, band, gen, t, beacon, station):
    """Satellites in view of both the station and the beacon at t, with each one's Doppler (Hz,
    uplink + downlink), delay (ms) and geometry."""
    out = []
    f_dn = T_XLATE[band] + F_NOM[gen]
    for label, name, sat in sats:
        st = sat_state(sat, t)
        if st is None:
            continue
        r, v = st
        r_dn, rr_dn, el_st = look(r, v, station)
        if el_st < MIN_EL:
            continue
        r_up, rr_up, el_bc = look(r, v, beacon)
        if el_bc < 0:
            continue
        out.append({"label": label, "dop": -f_dn * rr_dn / C - F_NOM[gen] * rr_up / C,
                    "delay": (r_up + r_dn) / C * 1e3, "el": el_st, "r": r, "band": band, "T": T_XLATE[band]})
    return out


def score(copies, cand, bands=None, sig=SIG[1]):
    """Best and runner-up assignment of copies [(f_hz, dt_ms)] to distinct candidates:
    (chi, chi_f, chi_t, perm, implied f0s, implied t_tx's). bands: each copy's band (a copy only
    goes to a satellite of its band); None: all one band. sig: the noise model, (Hz, ms). Copy 0's
    satellite anchors f0; the others only try satellites within WIN_HZ of it."""
    sig_f, sig_t = sig
    n, m = len(copies), len(cand)
    ok = [[bands is None or bands[k] == c["band"] for c in cand] for k in range(n)]
    f0 = [[copies[k][0] - c["T"] - c["dop"] for c in cand] for k in range(n)]
    tx = [[copies[k][1] - c["delay"] for c in cand] for k in range(n)]
    top = []

    def walk(k, perm, ref):
        if k == n:
            fs = [f0[i][j] for i, j in enumerate(perm)]
            ts = [tx[i][j] for i, j in enumerate(perm)]
            mf, mt = sum(fs) / n, sum(ts) / n
            cf = sum((x - mf) ** 2 for x in fs) / sig_f ** 2
            ct = sum((x - mt) ** 2 for x in ts) / sig_t ** 2
            top.append((cf + ct, cf, ct, tuple(perm), fs, ts))
            top.sort(key=lambda x: x[0])
            del top[2:]
            return
        for j in range(m):
            if ok[k][j] and j not in perm and abs(f0[k][j] - ref) < WIN_HZ:
                perm.append(j)
                walk(k + 1, perm, ref)
                perm.pop()

    for j0 in range(m):
        if ok[0][j0]:
            walk(1, [j0], f0[0][j0])
    return top


def store(con, gen, burst, bid, t, band, copies, cand, state, best=None, second=None, f0=None):
    n = len(copies)
    sats = ""
    if best is not None:
        chi, cf, ct, perm, fs, ts = best
        mf, mt = sum(fs) / n, sum(ts) / n
        if n == 1 and f0 is not None:
            mf = f0         # a single copy: its residual from the beacon's reference frequency
        names = []
        for k, (key, _, _) in enumerate(copies):
            c = cand[perm[k]]
            names.append(c["label"])
            la, lo, h = subpoint(c["r"])
            con.execute("UPDATE copies SET sat=?, el=?, f_res=?, t_res=?, slat=?, slon=?, salt=? "
                        "WHERE gen=? AND burst=? AND k=?",
                        (c["label"], round(c["el"], 1), round(fs[k] - mf, 1), round(ts[k] - mt, 3),
                         round(la, 3), round(lo, 3), round(h), gen, burst, key))
        sats = ",".join(names)
        margin = (second[0] - chi) if second else None
        con.execute("INSERT OR REPLACE INTO attr(gen, burst, bid, t, band, n, state, chi_f, chi_t, margin, f0, sats) "
                    "VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
                    (gen, burst, bid, t, band, n, state, round(cf, 2), round(ct, 2),
                     None if margin is None else round(min(margin, 1e6), 1), round(f0 if f0 is not None else mf, 1), sats))
    else:
        con.execute("INSERT OR REPLACE INTO attr(gen, burst, bid, t, band, n, state, chi_f, chi_t, margin, f0, sats) "
                    "VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
                    (gen, burst, bid, t, band, n, state, None, None, None, None, ""))


def pending(con, gen, settle_t):
    """Bursts not yet attributed whose beacon reports a position: [(burst, bid, t, band, lat, lon,
    [(k, f_hz, tc, band)])]."""
    table = "fgb_bursts" if gen == 1 else "sgb_bursts"
    rows = con.execute(
        f"SELECT b.id, b.bid, b.t, b.band, f.lat, f.lon FROM {table} b JOIN frames f ON f.hex=b.hex "
        "WHERE f.lat IS NOT NULL AND b.t<=? AND NOT EXISTS (SELECT 1 FROM attr a WHERE a.gen=? AND a.burst=b.id) "
        "AND EXISTS (SELECT 1 FROM copies c WHERE c.gen=? AND c.burst=b.id) ORDER BY b.t",
        (settle_t, gen, gen)).fetchall()
    out = []
    for bid_, bid, t, band, lat, lon in rows:
        if abs(lat - 1.0) < 1e-6 and abs(lon) < 1e-6:
            continue        # test transmitters' placeholder position (1 N 0 E), not where they are
        cp = [(k, f, tc, band_of(b, f)) for k, f, tc, b in con.execute(
            "SELECT k, f_hz, tc, band FROM copies WHERE gen=? AND burst=? ORDER BY k", (gen, bid_))]
        out.append((bid_, bid, t, band, lat, lon, cp))
    return out


def run(con, where, tle_dir, verbose=False):
    """Fetch TLEs when due, attribute the bursts that are ready, refresh sat_now. Returns a short
    summary for the rollup's log line, or None when attribution is off."""
    if Satrec is None or not where:
        return None
    init_tables(con)
    if (con.execute("SELECT value FROM meta WHERE key='attr_model'").fetchone() or [None])[0] != MODEL:
        con.execute("DELETE FROM attr")
        con.execute("UPDATE copies SET sat=NULL, el=NULL, f_res=NULL, t_res=NULL, slat=NULL, slon=NULL, salt=NULL")
        con.execute("INSERT INTO meta(key,value) VALUES('attr_model',?) ON CONFLICT(key) DO UPDATE SET value=excluded.value", (MODEL,))
        if verbose:
            print(f"satattr: noise model {MODEL}: attributing every burst again")
    try:
        os.makedirs(tle_dir, exist_ok=True)
        fetch_tles(con, tle_dir, verbose)
    except OSError as e:
        print(f"satattr: {tle_dir}: {e}", file=sys.stderr)
    tles = Tles(tle_dir)
    if not tles.snaps:
        return "no TLEs yet"
    station = ecef(where[0], where[1], 0.1)
    newest = con.execute("SELECT MAX(t) FROM (SELECT MAX(t) t FROM fgb_bursts UNION ALL "
                         "SELECT MAX(t) FROM sgb_bursts)").fetchone()[0] or 0
    counts = {}
    singles = []
    for gen in (1, 2):
        for burst, bid, t, band, lat, lon, cp in pending(con, gen, newest - SETTLE_S):
            bands = sorted({c[3] for c in cp})
            if any(b not in T_XLATE or not tles.at(t, b) for b in bands):
                continue        # no TLEs near enough in time: left for later
            if len(cp) == 1:
                singles.append((gen, burst, bid, t, cp[0][3], lat, lon, cp))
                continue
            t0 = min(c[2] for c in cp)
            copies = [(k, f, (tc - t0) * 1e3) for k, f, tc, _ in cp]
            cand = [c for b in bands for c in candidates(tles.at(t, b), b, gen, t, ecef(lat, lon), station)]
            if any(sum(c["band"] == b for c in cand) < sum(x[3] == b for x in cp) for b in bands):
                store(con, gen, burst, bid, t, band, copies, cand, "nosat")
                counts["nosat"] = counts.get("nosat", 0) + 1
                continue
            top = score([(f, dt) for _, f, dt in copies], cand, [x[3] for x in cp], SIG[gen])
            if not top:
                state, best, second = "nofit", None, None
            else:
                best, second = top[0], (top[1] if len(top) > 1 else None)
                fits = best[0] < 3 * len(copies)
                state = "nofit" if not fits else "clear" if (second is None or second[0] - best[0] > 9) else "ambiguous"
                if state == "nofit":
                    best = None
            store(con, gen, burst, bid, t, band, copies, cand, state, best, second)
            counts[state] = counts.get(state, 0) + 1
    # Single copies, by the beacon's f0 from its clearly attributed bursts within REF_S.
    for gen, burst, bid, t, band, lat, lon, cp in singles:
        ref = [r[0] for r in con.execute("SELECT f0 FROM attr WHERE gen=? AND bid=? AND state='clear' AND n>1 "
                                         "AND t BETWEEN ? AND ?", (gen, bid, t - REF_S, t + REF_S))]
        k, f, tc, _ = cp[0]
        copies = [(k, f, 0.0)]
        if not ref:
            if t < newest - REF_S:
                store(con, gen, burst, bid, t, band, copies, [], "noref")
                counts["noref"] = counts.get("noref", 0) + 1
            continue        # else wait: a later burst may give the reference
        f0 = statistics.median(ref)
        cand = candidates(tles.at(t, band), band, gen, t, ecef(lat, lon), station)
        fit = sorted(((f - T_XLATE[band] - c["dop"] - f0) ** 2 / SIG_REF[gen] ** 2, j) for j, c in enumerate(cand))
        if not fit or fit[0][0] > 9:
            state = "nofit" if fit else "nosat"
            store(con, gen, burst, bid, t, band, copies, cand, state)
            counts[state] = counts.get(state, 0) + 1
            continue
        chi, j = fit[0]
        second = fit[1][0] if len(fit) > 1 else None
        state = "clear" if second is None or second - chi > 9 else "ambiguous"
        best = (chi, chi, 0.0, (j,), [f0 + (f - T_XLATE[band] - cand[j]["dop"] - f0)], [0.0])
        store(con, gen, burst, bid, t, band, copies, cand, state, best,
              None if second is None else (second,), f0)
        counts[state] = counts.get(state, 0) + 1
    refresh_now(con, tles, station)
    con.commit()
    return ", ".join(f"{v} {k}" for k, v in sorted(counts.items())) or "nothing new"


def refresh_now(con, tles, station):
    """Current position of every satellite that has relayed a clearly attributed burst."""
    seen = {r[0] for r in con.execute("SELECT DISTINCT sat FROM copies WHERE sat IS NOT NULL")}
    now = time.time()
    con.execute("DELETE FROM sat_now")
    for band, sats in tles.latest().items():
        for label, name, sat in sats:
            if label not in seen:
                continue
            st = sat_state(sat, now)
            if st is None:
                continue
            r, v = st
            _, _, el = look(r, v, station)
            la, lo, h = subpoint(r)
            con.execute("INSERT OR REPLACE INTO sat_now VALUES(?,?,?,?,?,?,?,?,?)",
                        (label, name, "Galileo" if band == "galileo" else "GLONASS", now, round(la, 3), round(lo, 3),
                         round(h), round(el, 1), round(azimuth(r, station))))
