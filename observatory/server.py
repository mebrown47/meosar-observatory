#!/usr/bin/env python3
# Copyright (c) 2026 Mike Brown
# SPDX-License-Identifier: GPL-3.0-or-later
"""
MEOSAR Observatory: the dashboard server

A small read-only HTTP server over the database rollup.py builds from meosar-rx's archive.
Stdlib only, no external assets. Follows inmarsat-rx's observatory server.

  GET /                 -> the dashboard (tables)
  GET /globe            -> the globe: beacons on the left, their positions and tracks on a globe
  GET /api/overview     -> the dashboard's data (?hours=N, 0 = everything)
  GET /api/map          -> the globe's data (?hours=N)
  GET /api/beacon?id=B  -> one beacon: decoded fields and its bursts, with the satellites that relayed each
  GET /static/...       -> the globe's script and its vendored assets (three.js, Earth texture, coastlines)
  GET /healthz          -> ok

Usage:
    server.py [--db PATH] [--port 8896] [--host 0.0.0.0] [--station NAME] [--where LAT,LON]
"""

import argparse
import collections
import json
import os
import sqlite3
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

HERE = os.path.dirname(os.path.abspath(__file__))
DB_PATH = "observatory.db"
PORT = 8896


def connect_ro(path):
    return sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=5)


def newest(con):
    """The newest time in the database: windows count back from it, not the clock, so imported
    recordings show too."""
    r = con.execute("SELECT MAX(t) FROM (SELECT MAX(t) t FROM fgb_bursts UNION ALL SELECT MAX(t) FROM sgb_bursts "
                    "UNION ALL SELECT MAX(t) FROM rx_events)").fetchone()
    return r[0]


def frame_info(con, hexs):
    r = con.execute("SELECT protocol, btype, country, country_name, ident, lat, lon, fields FROM frames WHERE hex=?",
                    (hexs,)).fetchone()
    if not r:
        return {}
    return {"protocol": r[0], "btype": r[1], "country": r[2], "country_name": r[3], "ident": r[4],
            "lat": r[5], "lon": r[6], "fields": json.loads(r[7] or "[]")}


def has_attr(con):
    return con.execute("SELECT 1 FROM sqlite_master WHERE name='attr'").fetchone() is not None


def build_sats(con, t0):
    """Satellite attribution in the window: per satellite, the bursts it relayed (clearly
    attributed ones), and where it is now; and how the bursts with a known beacon position fared."""
    if not has_attr(con):
        return None
    states = dict(con.execute("SELECT state, COUNT(*) FROM attr WHERE t>=? GROUP BY state", (t0,)).fetchall())
    per = {}
    for sat, band, bid, burst, gen, t, el in con.execute(
            "SELECT c.sat, a.band, a.bid, a.burst, a.gen, a.t, c.el FROM attr a JOIN copies c ON c.gen=a.gen AND c.burst=a.burst "
            "WHERE a.t>=? AND a.state='clear' AND c.sat IS NOT NULL", (t0,)):
        p = per.setdefault(sat, {"sat": sat, "band": band, "bursts": 0, "beacons": set(), "els": [], "last": 0})
        p["bursts"] += 1
        p["beacons"].add(bid)
        p["els"].append(el)
        p["last"] = max(p["last"], t)
    now = {r[0]: r for r in con.execute("SELECT sat, name, constellation, t, lat, lon, alt, el, az FROM sat_now")}
    out = []
    for sat in set(per) | set(now):
        p = per.get(sat, {"sat": sat, "band": None, "bursts": 0, "beacons": set(), "els": [], "last": None})
        n = now.get(sat)
        els = sorted(p["els"])
        out.append({"sat": sat, "name": n[1] if n else sat,
                    "constellation": n[2] if n else {"galileo": "Galileo", "glonass": "GLONASS"}.get(p["band"], ""),
                    "bursts": p["bursts"], "beacons": len(p["beacons"]), "el_med": els[len(els) // 2] if els else None,
                    "last": p["last"], "now": {"t": n[3], "lat": n[4], "lon": n[5], "alt": n[6], "el": n[7], "az": n[8]} if n else None})
    out.sort(key=lambda x: (-x["bursts"], x["sat"]))
    # The beacon frequency the copies imply, from second-generation bursts (406.05 MHz nominal,
    # and tightly held): its offset is the receiver's own frequency reference error.
    offs = sorted(r[0] - 406.05e6 for r in con.execute(
        "SELECT f0 FROM attr WHERE gen=2 AND state='clear' AND n>1 AND t>=?", (t0,)))
    return {"states": states, "sats": out, "ref_hz": offs[len(offs) // 2] if offs else None, "ref_n": len(offs)}


def burst_sats(con, gen, ids):
    """{burst id: (state, "SAT1,SAT2")} for bursts that were looked at."""
    if not ids or not has_attr(con):
        return {}
    out = {}
    ids = list(ids)
    for i in range(0, len(ids), 500):
        chunk = ids[i:i + 500]
        out.update({r[0]: (r[1], r[2]) for r in con.execute(
            f"SELECT burst, state, sats FROM attr WHERE gen=? AND burst IN ({','.join('?' * len(chunk))})", (gen, *chunk))})
    return out


def build_overview(db, hours):
    con = connect_ro(db)
    try:
        t1 = newest(con)
        if t1 is None:
            return {"empty": True, "generated": time.time()}
        hours = int(hours or 24)
        t0 = t1 - hours * 3600 if hours > 0 else 0
        span_h = hours if hours > 0 else max(1, (t1 - (con.execute("SELECT MIN(t) FROM fgb_bursts").fetchone()[0] or t1)) / 3600)
        # Bin width for the activity chart: about 100 bars.
        bin_s = next(b for b in (300, 900, 1800, 3600, 3 * 3600, 6 * 3600, 86400, 7 * 86400) if span_h * 3600 / b <= 120)
        tstart = (t0 if t0 else t1 - span_h * 3600) // bin_s * bin_s

        # First-generation beacons: confirmed by a strict single decode or two separate bursts.
        rows = con.execute(
            "SELECT bid, COUNT(*), MAX(single), MAX(ebn0), MAX(relays), MIN(t), MAX(t), SUM(selftest), "
            "SUM(combined_only), (SELECT hex FROM fgb_bursts f2 WHERE f2.bid=f.bid AND t>=? ORDER BY t DESC LIMIT 1) "
            "FROM fgb_bursts f WHERE t>=? GROUP BY bid", (t0, t0)).fetchall()
        beacons = []
        for bid, n, single, eb, rel, first, last, st, co, hexs in rows:
            if not single and n < 2:
                continue
            fi = frame_info(con, hexs)
            beacons.append({"id": bid, "hex": hexs, "bursts": n, "ebn0": eb, "relays": rel, "first": first,
                            "last": last, "selftest": st > 0.8 * n, "combined_only": co, **fi})
        beacons.sort(key=lambda b: -b["last"])
        ok_ids = {b["id"] for b in beacons}

        sg = con.execute(
            "SELECT bid, COUNT(*), MAX(ebn0), MAX(relays), MIN(t), MAX(t), MAX(test), SUM(alone=0), tac, serial, "
            "country, (SELECT hex FROM sgb_bursts s2 WHERE s2.bid=s.bid AND t>=? ORDER BY t DESC LIMIT 1), MAX(single) "
            "FROM sgb_bursts s WHERE t>=? GROUP BY bid", (t0, t0)).fetchall()
        sgb = []
        for bid, n, eb, rel, first, last, test, co, tac, serial, country, hexs, single in sg:
            if not single and n < 2:
                continue
            fi = frame_info(con, hexs)
            sgb.append({"id": bid, "hex": hexs, "bursts": n, "ebn0": eb, "relays": rel, "first": first, "last": last,
                        "test": bool(test), "combined_only": co, "tac": tac, "serial": serial, **fi,
                        "country": country})
        sgb.sort(key=lambda b: -b["last"])

        # Activity: confirmed-beacon bursts and SGB bursts per bin.
        act = collections.defaultdict(lambda: [0, 0])
        for t, bid in con.execute("SELECT t, bid FROM fgb_bursts WHERE t>=?", (tstart,)):
            if bid in ok_ids:
                act[int(t // bin_s * bin_s)][0] += 1
        sgb_ids = {s["id"] for s in sgb}
        for t, bid in con.execute("SELECT t, bid FROM sgb_bursts WHERE t>=?", (tstart,)):
            if bid in sgb_ids:
                act[int(t // bin_s * bin_s)][1] += 1
        nb = int((t1 - tstart) // bin_s) + 1
        activity = [[int(tstart + i * bin_s), *act.get(int(tstart + i * bin_s), [0, 0])] for i in range(nb)]

        relays = collections.Counter()
        n_bursts = only_comb = by_repeat = 0
        for rel, co, bid, via in con.execute("SELECT relays, combined_only, bid, COALESCE(via, '') FROM fgb_bursts WHERE t>=?", (t0,)):
            if bid in ok_ids:
                relays[min(rel, 8)] += 1
                n_bursts += 1
                if via:
                    by_repeat += 1
                else:
                    only_comb += co
        cand = con.execute("SELECT band, gen, SUM(cand), SUM(acc) FROM cand_hourly WHERE hour>=? GROUP BY band, gen",
                           (int(t0) // 3600 * 3600,)).fetchall()
        bands = collections.defaultdict(lambda: {"fgb_cand": 0, "fgb_acc": 0, "sgb_cand": 0, "sgb_acc": 0})
        for band, gen, c, a in cand:
            k = "fgb" if gen == 1 else "sgb"
            bands[band][k + "_cand"] += c or 0
            bands[band][k + "_acc"] += a or 0
        bt = collections.Counter(b.get("btype") or "other" for b in beacons)
        countries = {b.get("country") for b in beacons} | {b.get("country") for b in sgb}

        # Receiver: its last start and status record.
        rx = {}
        st = con.execute("SELECT t, cand, acc, sgb, dropped_s, cpu FROM rx_events WHERE kind='stat' ORDER BY t DESC LIMIT 1").fetchone()
        sr = con.execute("SELECT t, src FROM rx_events WHERE kind='start' ORDER BY t DESC LIMIT 1").fetchone()
        if st:
            rx = {"last": st[0], "cand": st[1], "acc": st[2], "sgb": st[3], "dropped_s": st[4], "cpu": st[5]}
        if sr:
            rx.update(started=sr[0], src=sr[1])
        rx["starts"] = con.execute("SELECT COUNT(*) FROM rx_events WHERE kind='start' AND t>=?", (t0,)).fetchone()[0]
        meta = dict(con.execute("SELECT key, value FROM meta").fetchall())
        return {
            "generated": time.time(), "newest": t1, "hours": hours, "bin_s": bin_s,
            "totals": {"beacons": len(beacons), "bursts": n_bursts, "only_combined": only_comb, "by_repeat": by_repeat,
                       "sgb_beacons": len(sgb), "sgb_bursts": sum(s["bursts"] for s in sgb),
                       "countries": len(countries - {None}), "types": dict(bt.most_common())},
            "beacons": beacons, "sgb": sgb, "activity": activity, "relays": dict(sorted(relays.items())),
            "bands": dict(bands), "rx": rx, "last_rollup": float(meta.get("last_run", 0)),
            "dec406": meta.get("dec406") == "yes", "satellites": build_sats(con, t0),
        }
    finally:
        con.close()


def build_beacon(db, bid):
    con = connect_ro(db)
    try:
        if "-" in (bid or ""):
            rows = con.execute("SELECT b.t, b.band, b.hex, b.ebn0, b.relays, b.alone=0, b.f_hz, b.ebn0s, b.how, "
                               "f.lat, f.lon, b.id FROM sgb_bursts b LEFT JOIN frames f ON f.hex=b.hex "
                               "WHERE b.bid=? ORDER BY b.t DESC LIMIT 1000", (bid,)).fetchall()
            sats = burst_sats(con, 2, [r[11] for r in rows])
            bursts = [{"t": r[0], "band": r[1], "hex": r[2], "ebn0": r[3], "relays": r[4], "combined_only": r[5],
                       "f_hz": r[6], "detail": f"{r[8]} · relays {r[7]} dB", "lat": r[9], "lon": r[10],
                       "attr": sats.get(r[11])} for r in rows]
        else:
            rows = con.execute("SELECT b.t, b.band, b.hex, b.ebn0, b.relays, b.combined_only, b.f_hz, b.cn0, b.single, "
                               "f.lat, f.lon, COALESCE(b.via, ''), b.id FROM fgb_bursts b LEFT JOIN frames f ON f.hex=b.hex "
                               "WHERE b.bid=? ORDER BY b.t DESC LIMIT 1000", (bid,)).fetchall()
            sats = burst_sats(con, 1, [r[12] for r in rows])
            bursts = [{"t": r[0], "band": r[1], "hex": r[2], "ebn0": r[3], "relays": r[4], "combined_only": r[5],
                       "f_hz": r[6], "detail": {"match": "matched to this beacon's frame", "repeat": "decoded with its repeats"}.get(
                           r[11], ("strict" if r[8] else "weak") + (f" · C/N0 {r[7]:.1f}" if r[7] else "")),
                       "lat": r[9], "lon": r[10], "attr": sats.get(r[12])} for r in rows]
        if not bursts:
            return {"id": bid, "bursts": []}
        # The latest frame's fields; a beacon that updates its position (or an orbitography
        # beacon's time) changes its frame, and the bursts carry each frame's position.
        nframes = len({b["hex"] for b in bursts})
        return {"id": bid, "bursts": bursts, "latest": {"hex": bursts[0]["hex"], **frame_info(con, bursts[0]["hex"])},
                "frames": nframes}
    finally:
        con.close()


def build_map(db, hours, where):
    """Confirmed beacons of both generations with their last known position and track."""
    ov = build_overview(db, hours)
    if ov.get("empty"):
        return {**ov, "station": where}
    t0 = ov["newest"] - ov["hours"] * 3600 if ov["hours"] > 0 else 0
    con = connect_ro(db)
    try:
        tracks = {}
        for table in ("fgb_bursts", "sgb_bursts"):
            for bid, t, lat, lon in con.execute(
                    f"SELECT b.bid, b.t, f.lat, f.lon FROM {table} b JOIN frames f ON f.hex=b.hex "
                    "WHERE b.t>=? AND f.lat IS NOT NULL ORDER BY b.t", (t0,)):
                tr = tracks.setdefault(bid, [])
                p = [round(t), round(lat, 4), round(lon, 4)]
                if tr and tr[-1][1:] == p[1:]:
                    tr[-1][0] = p[0]          # same place again: keep the latest time
                else:
                    tr.append(p)
        # Each beacon's latest clearly attributed burst: the satellites (sub-satellite points at the
        # time) that relayed it, for the globe's beacon -> satellite -> station links.
        relay = {}
        if has_attr(con):
            latest = {}
            for gen, bid, burst, t in con.execute("SELECT gen, bid, burst, t FROM attr WHERE state='clear' AND t>=?", (t0,)):
                if bid not in latest or t > latest[bid][2]:
                    latest[bid] = (gen, burst, t)
            for bid, (gen, burst, t) in latest.items():
                relay[bid] = {"t": t, "sats": [{"sat": r[0], "lat": r[1], "lon": r[2], "alt": r[3], "el": r[4]} for r in con.execute(
                    "SELECT sat, slat, slon, salt, el FROM copies WHERE gen=? AND burst=? AND sat IS NOT NULL ORDER BY k", (gen, burst))]}
        sats = build_sats(con, t0)
    finally:
        con.close()
    items = []
    for gen, lst in ((1, ov["beacons"]), (2, ov["sgb"])):
        for b in lst:
            tr = tracks.get(b["id"], [])
            items.append({"id": b["id"], "gen": gen, "btype": "test" if gen == 2 and b.get("test") else b.get("btype"),
                          "ident": b.get("ident") or b.get("protocol"), "country": b.get("country"),
                          "country_name": b.get("country_name"), "bursts": b["bursts"], "ebn0": b["ebn0"],
                          "relays": b["relays"], "first": b["first"], "last": b["last"], "selftest": b.get("selftest"),
                          "fields": b.get("fields", []), "hex": b.get("hex"),
                          "tac": b.get("tac"), "serial": b.get("serial"),
                          "lat": tr[-1][1] if tr else None, "lon": tr[-1][2] if tr else None,
                          "track": tr[-500:], "relay": relay.get(b["id"])})
    items.sort(key=lambda b: -b["last"])
    return {"generated": ov["generated"], "newest": ov["newest"], "hours": ov["hours"], "items": items,
            "station": where, "rx": ov["rx"], "last_rollup": ov["last_rollup"], "satellites": sats}


STATIC_TYPES = {".js": "text/javascript", ".json": "application/json", ".jpg": "image/jpeg", ".css": "text/css",
                ".md": "text/plain; charset=utf-8", "": "text/plain; charset=utf-8"}


class Handler(BaseHTTPRequestHandler):
    db_path = DB_PATH
    station = ""
    where = None
    server_version = "meosar-observatory"

    def log_message(self, *a):
        pass

    def _send(self, code, body, ctype, cache="no-store"):
        if isinstance(body, str):
            body = body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Cache-Control", cache)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, build, *a):
        try:
            self._send(200, json.dumps(build(*a)), "application/json")
        except sqlite3.OperationalError as e:   # no database yet: the rollup hasn't run
            self._send(200, json.dumps({"empty": True, "error": str(e), "generated": time.time()}), "application/json")
        except Exception as e:  # noqa: keep the server alive
            self._send(500, json.dumps({"error": str(e)}), "application/json")

    def do_GET(self):
        u = urlparse(self.path)
        q = parse_qs(u.query)
        arg = lambda k: (q.get(k, [None])[0] or None)
        if u.path in ("/", "/index.html", "/globe"):
            name = "globe.html" if u.path == "/globe" else "page.html"
            with open(os.path.join(HERE, name), encoding="utf-8") as f:
                self._send(200, f.read().replace("__STATION__", self.station), "text/html; charset=utf-8")
        elif u.path.startswith("/static/"):
            root = os.path.join(HERE, "static")
            path = os.path.realpath(os.path.join(root, u.path[len("/static/"):]))
            if not path.startswith(root + os.sep) or not os.path.isfile(path):
                self._send(404, "not found", "text/plain")
                return
            with open(path, "rb") as f:
                body = f.read()
            ext = os.path.splitext(path)[1]
            # The vendored assets never change; the page's own script may, with an update.
            cache = "max-age=86400" if "/vendor/" in path else "no-cache"
            self._send(200, body, STATIC_TYPES.get(ext, "application/octet-stream"), cache)
        elif u.path == "/api/map":
            self._json(build_map, self.db_path, arg("hours"), self.where)
        elif u.path == "/api/overview":
            self._json(build_overview, self.db_path, arg("hours"))
        elif u.path == "/api/beacon":
            self._json(build_beacon, self.db_path, arg("id"))
        elif u.path == "/healthz":
            self._send(200, "ok", "text/plain")
        else:
            self._send(404, "not found", "text/plain")

    do_HEAD = do_GET


def main():
    ap = argparse.ArgumentParser(description="MEOSAR Observatory server")
    ap.add_argument("--db", default=DB_PATH)
    ap.add_argument("--port", type=int, default=PORT)
    ap.add_argument("--host", default="0.0.0.0")
    ap.add_argument("--station", default="", help="station name for the page header")
    ap.add_argument("--where", default="", help="station position LAT,LON for the globe (empty = not shown)")
    a = ap.parse_args()
    Handler.db_path = a.db
    Handler.station = a.station
    try:
        lat, lon = (float(x) for x in a.where.split(","))
        Handler.where = {"lat": lat, "lon": lon}
    except ValueError:
        Handler.where = None
    srv = ThreadingHTTPServer((a.host, a.port), Handler)
    print(f"meosar observatory: http://{a.host}:{a.port}/  db={a.db}", flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
