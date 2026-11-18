#!/usr/bin/env python3
# Copyright (c) 2026 Mike Brown
# SPDX-License-Identifier: GPL-3.0-or-later
"""
MEOSAR Observatory — archive rollup

Folds the JSONL archive written by meosar-rx --archive into a SQLite database the dashboard
queries: one row per decoded burst (relayed copies of a burst merged into one), the decoded
fields of every distinct frame, candidates per hour, and the receiver's status records.
Follows inmarsat-rx's and vdl2-observatory's rollup.py.

  * Incremental: per-file byte offsets in `ingest_state`; each run reads only appended bytes,
    and only whole lines. Re-running is safe. Deleting the database re-ingests everything.
  * Times: every record's "utc" becomes epoch seconds (REAL). Records without one (a raw file
    with no start time) are skipped, since their times can't be placed.
  * A burst relayed by several satellites gives several records within ~40 ms, and the same
    beacon repeats every ~50 s, so records of one beacon ID within 1.5 s are one burst.
  * A beacon is shown once it is confirmed: a strict decode, or bursts in two separate
    transmissions. For SGB, strict is a hard decode with at most 5 of BCH's 6 correctable errors:
    single Chase decodes and 6-error decodes are where noise gives false frames.
  * Repeat combining (meosar-rx fgb_match / fgb_repeat records) adds bursts that never decoded
    alone: each becomes a burst of its beacon, marked in `via`.
  * Frame fields (protocol, country, identification, position) come from dec406's hex decoder
    when it is installed (--dec406), once per distinct frame; without it, the basics are
    decoded here.

Usage:
    rollup.py [--archive DIR] [--db PATH] [--dec406 PATH] [-v] [FILE ...]

FILE arguments are ingested as well as the archive: results from recordings, say.
"""

import argparse
import fcntl
import glob
import json
import os
import re
import sqlite3
import subprocess
import sys
import time
from datetime import datetime

SCHEMA_VERSION = 1
GROUP_S = 1.5


def init_db(con):
    con.executescript(
        """
        PRAGMA journal_mode=WAL;
        PRAGMA synchronous=NORMAL;
        CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
        CREATE TABLE IF NOT EXISTS ingest_state (
            file TEXT PRIMARY KEY, offset INTEGER NOT NULL, size INTEGER NOT NULL, updated INTEGER NOT NULL);

        -- One row per first-generation burst, relayed copies merged. bid: the 15-hex beacon ID.
        -- relays: satellites it came through (distinct carriers, or the number combined).
        -- combined_only: no single copy decoded alone. single: some copy passed the strict test.
        CREATE TABLE IF NOT EXISTS fgb_bursts (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            t REAL, band TEXT, bid TEXT, hex TEXT, ebn0 REAL, relays INTEGER, fbins TEXT,
            comb INTEGER, combined_only INTEGER, single INTEGER, selftest INTEGER, f_hz REAL, cn0 REAL);
        CREATE INDEX IF NOT EXISTS ix_fgb_t ON fgb_bursts(t);
        CREATE INDEX IF NOT EXISTS ix_fgb_bid ON fgb_bursts(bid, t);

        -- One row per second-generation burst. bid: country-TAC-serial.
        CREATE TABLE IF NOT EXISTS sgb_bursts (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            t REAL, tq INTEGER, band TEXT, bid TEXT, hex TEXT, tac INTEGER, serial INTEGER,
            country INTEGER, test INTEGER, relays INTEGER, ebn0 REAL, ebn0s TEXT, how TEXT,
            alone INTEGER, single INTEGER, f_hz REAL,
            UNIQUE(hex, tq));
        CREATE INDEX IF NOT EXISTS ix_sgb_t ON sgb_bursts(t);
        CREATE INDEX IF NOT EXISTS ix_sgb_bid ON sgb_bursts(bid, t);

        -- Decoded fields of each distinct frame.
        CREATE TABLE IF NOT EXISTS frames (
            hex TEXT PRIMARY KEY, gen INTEGER, protocol TEXT, btype TEXT, country INTEGER,
            country_name TEXT, ident TEXT, lat REAL, lon REAL, fields TEXT);

        -- Candidates examined and accepted per hour, band and generation.
        CREATE TABLE IF NOT EXISTS cand_hourly (
            hour INTEGER, band TEXT, gen INTEGER, cand INTEGER, acc INTEGER,
            PRIMARY KEY (hour, band, gen));

        -- The receiver's start and status records (counts are totals since that start).
        CREATE TABLE IF NOT EXISTS rx_events (
            t REAL, kind TEXT, src TEXT, cand INTEGER, acc INTEGER, sgb INTEGER, dropped_s REAL, cpu REAL);
        CREATE INDEX IF NOT EXISTS ix_rx_t ON rx_events(t);
        """
    )
    # Added after schema 1: how a burst was decoded when not alone ('', 'match', 'repeat').
    if "via" not in [c[1] for c in con.execute("PRAGMA table_info(fgb_bursts)")]:
        con.execute("ALTER TABLE fgb_bursts ADD COLUMN via TEXT DEFAULT ''")
    con.execute("INSERT INTO meta(key,value) VALUES('schema',?) "
                "ON CONFLICT(key) DO UPDATE SET value=excluded.value", (str(SCHEMA_VERSION),))
    con.commit()


def epoch(utc):
    try:
        return datetime.fromisoformat(utc).timestamp()
    except (TypeError, ValueError):
        return None


# ---------- frame fields ----------

def bits_of(hexs, n):
    return [int(c) for c in bin(int(hexs, 16))[2:].zfill(n)]


def num(b, a, n):
    return int("".join(map(str, b[a:a + n])), 2)


# First-generation protocol codes (C/S T.001): user protocols (protocol flag 1) by their 3-bit
# code, location protocols (flag 0) by their 4-bit code.
FGB_USER = {0: ("Orbitography", "orbitography"), 1: ("ELT aviation user", "ELT"), 2: ("Maritime user", "EPIRB"),
            3: ("Serial user", None), 4: ("National user", None), 5: ("Second-generation spare", None),
            6: ("Radio call sign", "EPIRB"), 7: ("Test user", "test")}
FGB_LOC = {2: ("Standard location, EPIRB MMSI", "EPIRB"), 3: ("Standard location, ELT 24-bit address", "ELT"),
           4: ("Standard location, ELT serial", "ELT"), 5: ("Standard location, ELT operator", "ELT"),
           6: ("Standard location, EPIRB serial", "EPIRB"), 7: ("Standard location, PLB serial", "PLB"),
           8: ("National location, ELT", "ELT"), 9: ("ELT-DT location", "ELT"), 10: ("National location, EPIRB", "EPIRB"),
           11: ("National location, PLB", "PLB"), 12: ("Ship security", "SSAS"), 13: ("RLS location", "RLS"),
           14: ("Standard test", "test"), 15: ("National test", "test")}
SERIAL_TYPE = {0: "ELT", 1: "ELT", 2: "EPIRB", 3: "ELT", 4: "EPIRB", 6: "PLB"}


def fgb_basic(hexs):
    """(protocol, beacon type, country) from the bits alone. hexs: frame bits 25.. as hex."""
    b = [0] * 24 + bits_of(hexs, len(hexs) * 4)   # index = T.001 bit number - 1
    country, pf = num(b, 26, 10), b[25]
    if pf:
        code = num(b, 36, 3)
        name, btype = FGB_USER.get(code, (f"User protocol {code}", None))
        if code == 3:
            btype = SERIAL_TYPE.get(num(b, 39, 3))
    else:
        code = num(b, 36, 4)
        name, btype = FGB_LOC.get(code, (f"Location protocol {code}", None))
    return name, btype or "other", country


def sgb_basic(hexs):
    b = bits_of(hexs, 252)[2:]
    return {"tac": num(b, 0, 16), "serial": num(b, 16, 14), "country": num(b, 30, 10), "test": b[42]}


KEEP = ("Protocol", "Country", "Hex ID", "23 Hex ID", "Identification", "Type", "TAC", "Serial", "MMSI",
        "Aircraft", "Vessel", "Radio call sign", "Position", "Orbitography ID", "UTC time", "Test Protocol",
        "Homing", "Frame Mode", "Battery", "GNSS Status", "Elapsed time", "Last position", "Altitude")
POS_RE = re.compile(r"Position[^:]*:\s*([\d.]+)°?\s*([NS]),\s*([\d.]+)°?\s*([EW])")


def dec406_fields(tool, hexs):
    """Field lines, country name, identification and position from dec406's hex decoder."""
    if not tool:
        return None
    try:
        out = subprocess.run([tool, hexs], capture_output=True, text=True, timeout=10).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    lines, cname, ident, pos = [], None, None, None
    for l in out.splitlines():
        l = l.strip()
        if not l.startswith(KEEP) or l in lines or l.startswith("Coordinates"):
            continue
        lines.append(l)
        k, _, v = l.partition(":")
        v = v.strip()
        if k.strip() in ("Country", "Country Code") and "(" in v:
            cname = v[v.index("(") + 1:v.rindex(")")]
        elif k.strip() in ("Identification", "Type") and not ident:
            ident = v
        m = POS_RE.search(l)
        if m and pos is None:
            lat = float(m.group(1)) * (1 if m.group(2) == "N" else -1)
            lon = float(m.group(3)) * (1 if m.group(4) == "E" else -1)
            if (lat, lon) != (0.0, 0.0) and abs(lat) <= 90 and abs(lon) <= 180:
                pos = (lat, lon)
    return {"lines": lines[:16], "country_name": cname, "ident": ident, "pos": pos}


def ensure_frame(con, hexs, gen, tool, seen):
    if hexs in seen:
        return
    seen.add(hexs)
    if con.execute("SELECT 1 FROM frames WHERE hex=?", (hexs,)).fetchone():
        return
    if gen == 1:
        protocol, btype, country = fgb_basic(hexs)
    else:
        s = sgb_basic(hexs)
        protocol, btype, country = "Second generation", "test" if s["test"] else "SGB", s["country"]
    d = dec406_fields(tool, hexs) or {}
    if gen == 2 and d.get("ident"):
        btype = "test" if btype == "test" else d["ident"].split()[0]
    pos = d.get("pos") or (None, None)
    con.execute("INSERT INTO frames VALUES(?,?,?,?,?,?,?,?,?,?)",
                (hexs, gen, protocol, btype, country, d.get("country_name"), d.get("ident"), pos[0], pos[1],
                 json.dumps(d.get("lines", []))))


# ---------- records ----------

class Ingest:
    def __init__(self, con, tool):
        self.con, self.tool, self.seen = con, tool, set()
        self.cand = {}   # (hour, band, gen) -> [cand, acc]
        self.lines = self.bad = self.skipped = 0

    def count(self, t, band, gen, acc, cand=1):
        a = self.cand.setdefault((int(t) // 3600 * 3600, band or "?", gen), [0, 0])
        a[0] += cand
        a[1] += 1 if acc else 0

    def record(self, r):
        kind = r.get("kind")
        t = epoch(r.get("utc"))
        if t is None:
            self.skipped += 1
            return
        if kind in ("start", "stat"):
            self.con.execute("INSERT INTO rx_events VALUES(?,?,?,?,?,?,?,?)",
                             (t, kind, r.get("src"), r.get("cand"), r.get("acc"), r.get("sgb"),
                              r.get("dropped_s"), r.get("cpu")))
        elif kind in ("fgb_match", "fgb_repeat"):
            self.repeat(r, t)
        elif kind == "sgb" or ("relays" in r and "t_cw" not in r):   # (the Python scan's have no kind)
            self.sgb(r, t)
        elif "t_cw" in r:
            self.fgb(r, t)

    def fgb(self, r, t):
        combined = "combined" in r
        acc = r.get("accept")
        self.count(t, r.get("band"), 1, acc, cand=0 if combined else 1)
        hexs = r.get("hex")
        if not acc or not hexs:
            return
        self.burst(t, r.get("band"), hexs, None if combined else r.get("ebn0"), r.get("ebn0"),
                   "" if combined else str(round(r.get("f_hz", 0) / 100)), r.get("combined", 0), int(combined),
                   1 if acc == "single" else 0, int(bool(r.get("selftest_sync"))), r.get("f_hz"), r.get("cn0"), "")

    def repeat(self, r, t):
        """Repeat combining: a weak burst attributed to a known beacon (fgb_match), or a frame decoded
        from several weak bursts (fgb_repeat), each member a burst of that beacon."""
        hexs = r.get("hex")
        if not hexs:
            return
        if r.get("kind") == "fgb_match":
            members, via = [r], "match"
        else:
            members, via = r.get("members") or [], "repeat"
        for m in members:
            tm = epoch(m.get("utc")) or t
            self.burst(tm, m.get("band"), hexs, m.get("ebn0"), m.get("ebn0"), str(round(m.get("f_hz", 0) / 100)),
                       0, 1, 0, 0, m.get("f_hz"), None, via)

    def burst(self, t, band, hexs, ebn0, ebn0_any, fb, comb, combined, single, selftest, f_hz, cn0, via):
        """One copy of a burst: merged with the beacon's burst within GROUP_S, or a new burst row.
        combined: this copy was not decoded alone (relay-combined, matched or repeat-decoded)."""
        bid = hexs[:15]
        ensure_frame(self.con, hexs, 1, self.tool, self.seen)
        row = self.con.execute("SELECT id, fbins, comb, ebn0, combined_only, single FROM fgb_bursts "
                               "WHERE bid=? AND t BETWEEN ? AND ? ORDER BY ABS(t-?) LIMIT 1",
                               (bid, t - GROUP_S, t + GROUP_S, t)).fetchone()
        if row:
            bid_id, fbins, c0, e0, co, s0 = row
            bins = set(filter(None, (fbins or "").split(","))) | ({fb} if fb else set())
            c = max(c0 or 0, comb)
            e = max(x for x in (e0, ebn0, -99) if x is not None)
            self.con.execute("UPDATE fgb_bursts SET fbins=?, comb=?, relays=?, ebn0=?, combined_only=?, single=? "
                             "WHERE id=?", (",".join(sorted(bins)), c, max(len(bins), c, 1),
                                            e if e > -99 else None, co and combined, max(s0, single), bid_id))
        else:
            self.con.execute("INSERT INTO fgb_bursts(t,band,bid,hex,ebn0,relays,fbins,comb,combined_only,single,"
                             "selftest,f_hz,cn0,via) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                             (t, band, bid, hexs, ebn0 if ebn0 is not None else ebn0_any,
                              max(1 if fb else 0, comb, 1), fb, comb, combined, single, selftest, f_hz, cn0, via))

    def sgb(self, r, t):
        acc = r.get("accept", "single" if r.get("ok") else "")   # Python scan: ok means accepted
        self.count(t, r.get("band"), 2, acc and r.get("ok"))
        hexs = r.get("hex")
        if not (acc and r.get("ok") and hexs):
            return
        s = sgb_basic(hexs)
        if not any(bits_of(hexs, 252)[2:204]):
            return
        ensure_frame(self.con, hexs, 2, self.tool, self.seen)
        rel = r.get("relays") or []
        eb = [x.get("ebn0") for x in rel if x.get("ebn0") is not None]
        self.con.execute(
            "INSERT OR IGNORE INTO sgb_bursts(t,tq,band,bid,hex,tac,serial,country,test,relays,ebn0,ebn0s,how,"
            "alone,single,f_hz) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (t, round(t * 10), r.get("band"), f"{s['country']}-{s['tac']}-{s['serial']}", hexs, s["tac"],
             s["serial"], s["country"], s["test"], len(rel) or 1, max(eb) if eb else None,
             " ".join(f"{e:.1f}" for e in eb), r.get("how", ""), int(any(x.get("alone") for x in rel)),
             int(acc == "single" and "6 errors" not in r.get("how", "")), r.get("f_hz")))

    def flush(self):
        self.con.executemany(
            "INSERT INTO cand_hourly VALUES(?,?,?,?,?) ON CONFLICT(hour,band,gen) DO UPDATE SET "
            "cand=cand+excluded.cand, acc=acc+excluded.acc", [(*k, *v) for k, v in self.cand.items()])
        self.cand = {}


def ingest_file(con, path, ing, verbose):
    size = os.path.getsize(path)
    row = con.execute("SELECT offset,size FROM ingest_state WHERE file=?", (path,)).fetchone()
    offset = row[0] if row else 0
    if row and size < row[1]:
        offset = 0   # replaced or truncated: read it again
    if offset >= size:
        return 0
    with open(path, "rb") as f:
        f.seek(offset)
        buf = f.read()
    *complete, tail = buf.split(b"\n")
    consumed = len(buf) - len(tail)
    n = 0
    for raw in complete:
        if not raw.strip():
            continue
        try:
            ing.record(json.loads(raw))
            n += 1
        except (ValueError, KeyError, TypeError):
            ing.bad += 1
    ing.lines += n
    ing.flush()
    con.execute("INSERT INTO ingest_state(file,offset,size,updated) VALUES(?,?,?,?) "
                "ON CONFLICT(file) DO UPDATE SET offset=excluded.offset, size=excluded.size, "
                "updated=excluded.updated", (path, offset + consumed, size, int(time.time())))
    con.commit()   # per file, so a long first ingest shows progress and survives an interruption
    if verbose:
        print(f"  {os.path.basename(path)}: +{n} records ({offset}->{offset + consumed} of {size} bytes)")
    return n


def main():
    ap = argparse.ArgumentParser(description="MEOSAR Observatory rollup")
    ap.add_argument("--archive", default="archive", help="directory of meosar-YYYY-MM-DD.jsonl files")
    ap.add_argument("--db", default="observatory.db")
    ap.add_argument("--dec406", default=os.path.expanduser("~/dec406/build/dec406_hex"),
                    help="dec406's hex decoder, for frame fields and positions ('' = none)")
    ap.add_argument("--verbose", "-v", action="store_true")
    ap.add_argument("files", nargs="*", help="more JSONL files to ingest (recordings' results)")
    a = ap.parse_args()

    extra = []
    for f in a.files:
        if os.path.isfile(f):
            extra.append(os.path.abspath(f))
        else:
            print(f"rollup: skipping {f!r}: no such file", file=sys.stderr)
    files = sorted(glob.glob(os.path.join(a.archive, "meosar-*.jsonl"))) + extra
    if not files:
        print(f"rollup: no archive files in {a.archive}", file=sys.stderr)
        return 0
    tool = a.dec406 if a.dec406 and os.access(a.dec406, os.X_OK) else None
    os.makedirs(os.path.dirname(os.path.abspath(a.db)), exist_ok=True)
    lock = open(a.db + ".lock", "w")
    fcntl.flock(lock, fcntl.LOCK_EX)   # one rollup at a time, or two would add the same lines
    con = sqlite3.connect(a.db, timeout=30)
    init_db(con)
    ing = Ingest(con, tool)
    total = sum(ingest_file(con, p, ing, a.verbose) for p in files)
    con.execute("INSERT INTO meta(key,value) VALUES('last_run',?) ON CONFLICT(key) DO UPDATE SET "
                "value=excluded.value", (str(int(time.time())),))
    con.execute("INSERT INTO meta(key,value) VALUES('dec406',?) ON CONFLICT(key) DO UPDATE SET "
                "value=excluded.value", ("yes" if tool else "no",))
    con.commit()
    con.close()
    extra = [f"{ing.bad} malformed"] * bool(ing.bad) + [f"{ing.skipped} without a UTC time"] * bool(ing.skipped)
    print(f"rollup: {total} new records from {len(files)} file(s)" + (", " + ", ".join(extra) if extra else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
