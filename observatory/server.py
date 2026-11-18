#!/usr/bin/env python3
# Copyright (c) 2026 Mike Brown
# SPDX-License-Identifier: GPL-3.0-or-later
"""
MEOSAR Observatory: the dashboard server

A small read-only HTTP server over the database rollup.py builds from meosar-rx's archive.
Stdlib only, no external assets. Follows inmarsat-rx's observatory server.

  GET /                 -> the dashboard
  GET /api/overview     -> its data (?hours=N, 0 = everything)
  GET /api/beacon?id=B  -> one beacon: decoded fields and its bursts
  GET /healthz          -> ok

Usage:
    server.py [--db PATH] [--port 8896] [--host 0.0.0.0] [--station NAME]
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
        n_bursts = only_comb = 0
        for rel, co, bid in con.execute("SELECT relays, combined_only, bid FROM fgb_bursts WHERE t>=?", (t0,)):
            if bid in ok_ids:
                relays[min(rel, 8)] += 1
                n_bursts += 1
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
            "totals": {"beacons": len(beacons), "bursts": n_bursts, "only_combined": only_comb,
                       "sgb_beacons": len(sgb), "sgb_bursts": sum(s["bursts"] for s in sgb),
                       "countries": len(countries - {None}), "types": dict(bt.most_common())},
            "beacons": beacons, "sgb": sgb, "activity": activity, "relays": dict(sorted(relays.items())),
            "bands": dict(bands), "rx": rx, "last_rollup": float(meta.get("last_run", 0)),
            "dec406": meta.get("dec406") == "yes",
        }
    finally:
        con.close()


def build_beacon(db, bid):
    con = connect_ro(db)
    try:
        if "-" in (bid or ""):
            rows = con.execute("SELECT b.t, b.band, b.hex, b.ebn0, b.relays, b.alone=0, b.f_hz, b.ebn0s, b.how, "
                               "f.lat, f.lon FROM sgb_bursts b LEFT JOIN frames f ON f.hex=b.hex "
                               "WHERE b.bid=? ORDER BY b.t DESC LIMIT 1000", (bid,)).fetchall()
            bursts = [{"t": r[0], "band": r[1], "hex": r[2], "ebn0": r[3], "relays": r[4], "combined_only": r[5],
                       "f_hz": r[6], "detail": f"{r[8]} · relays {r[7]} dB", "lat": r[9], "lon": r[10]} for r in rows]
        else:
            rows = con.execute("SELECT b.t, b.band, b.hex, b.ebn0, b.relays, b.combined_only, b.f_hz, b.cn0, b.single, "
                               "f.lat, f.lon FROM fgb_bursts b LEFT JOIN frames f ON f.hex=b.hex "
                               "WHERE b.bid=? ORDER BY b.t DESC LIMIT 1000", (bid,)).fetchall()
            bursts = [{"t": r[0], "band": r[1], "hex": r[2], "ebn0": r[3], "relays": r[4], "combined_only": r[5],
                       "f_hz": r[6], "detail": ("strict" if r[8] else "weak") + (f" · C/N0 {r[7]:.1f}" if r[7] else ""),
                       "lat": r[9], "lon": r[10]} for r in rows]
        if not bursts:
            return {"id": bid, "bursts": []}
        # The latest frame's fields; a beacon that updates its position (or an orbitography
        # beacon's time) changes its frame, and the bursts carry each frame's position.
        nframes = len({b["hex"] for b in bursts})
        return {"id": bid, "bursts": bursts, "latest": {"hex": bursts[0]["hex"], **frame_info(con, bursts[0]["hex"])},
                "frames": nframes}
    finally:
        con.close()


class Handler(BaseHTTPRequestHandler):
    db_path = DB_PATH
    station = ""
    server_version = "meosar-observatory"

    def log_message(self, *a):
        pass

    def _send(self, code, body, ctype):
        if isinstance(body, str):
            body = body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Cache-Control", "no-store")
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
        if u.path in ("/", "/index.html"):
            with open(os.path.join(HERE, "page.html"), encoding="utf-8") as f:
                self._send(200, f.read().replace("__STATION__", self.station), "text/html; charset=utf-8")
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
    a = ap.parse_args()
    Handler.db_path = a.db
    Handler.station = a.station
    srv = ThreadingHTTPServer((a.host, a.port), Handler)
    print(f"meosar observatory: http://{a.host}:{a.port}/  db={a.db}", flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
