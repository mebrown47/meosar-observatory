#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""meosar-rx dashboard: results and build progress from the files the pipeline writes.

    server.py [--config dash.json] [--port 8896] [--host 0.0.0.0]
    server.py --full [--port 8897]        full detail, loopback only (reach it by SSH tunnel)

The page is meant to be reachable from outside, so it shows full detail only for beacons that
are not operational: orbitography (calibration), test protocols and self-test transmissions.
Every other beacon (it may be a real distress beacon) is counted, never identified.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import glob
import json
import os
import re
import subprocess
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

HERE = Path(__file__).resolve().parent
CACHE_S = 15
_cache = {"t": 0.0, "state": None}
_lock = threading.Lock()
MIDS: set[int] = set()
FULL = False
LOOPBACK = {"127.0.0.1", "::1", "localhost"}
_decoded: dict[str, list[str]] = {}
KEEP = ("Protocol", "Country", "Hex ID", "23 Hex ID", "Identification", "Type", "TAC", "Serial",
        "MMSI", "Aircraft", "Vessel", "Position", "Latitude", "Longitude", "Orbitography ID",
        "UTC time", "Test Protocol", "Homing", "Frame Mode")


def decode_fields(hexs: str, cfg) -> list[str]:
    """Field text from dec406's hex decoder (full mode only; cached per frame)."""
    if hexs in _decoded:
        return _decoded[hexs]
    tool = exp(cfg.get("dec406_hex", "~/dec406/build/dec406_hex"))
    try:
        out = subprocess.run([tool, hexs], capture_output=True, text=True, timeout=10).stdout
    except (OSError, subprocess.SubprocessError):
        out = ""
    lines = []
    for l in out.splitlines():
        l = l.strip()
        if l.startswith(KEEP) and l not in lines:
            lines.append(l)
    _decoded[hexs] = lines[:12]
    return _decoded[hexs]


def exp(p: str) -> str:
    return os.path.expanduser(p)


def load_mids():
    src = HERE.parent / "src" / "mids.c"
    try:
        body = src.read_text().split("{", 1)[1].split("}", 1)[0]
        MIDS.update(int(x) for x in body.split(",") if x.strip())
    except (OSError, IndexError, ValueError):
        pass


# ---------- beacon classes (what may be shown publicly) ----------

def fgb_bits(hexs: str):
    v = bin(int(hexs, 16))[2:].zfill(len(hexs) * 4)
    return [1] * 15 + [0, 0, 0, 1, 0, 1, 1, 1, 1] + [int(c) for c in v]


def fgb_class(hexs: str) -> tuple[str, int]:
    """('orbitography' | 'test' | 'operational', country), from the protocol code only. A
    self-test transmission of an operational beacon still identifies a real person's beacon, so
    self-test does not make a beacon public."""
    b = fgb_bits(hexs)
    country = int("".join(map(str, b[26:36])), 2)
    pf = b[25]
    if pf:
        code = int("".join(map(str, b[36:39])), 2)
        if code == 0:
            return "orbitography", country
        if code == 7:
            return "test", country
    else:
        code = int("".join(map(str, b[36:40])), 2)
        if code in (14, 15):
            return "test", country
    return "operational", country


def sgb_fields(hexs: str):
    b = [int(c) for c in bin(int(hexs, 16))[2:].zfill(252)][2:]
    num = lambda a, n: int("".join(map(str, b[a:a + n])), 2)
    return {"tac": num(0, 16), "serial": num(16, 14), "country": num(30, 10), "test": bool(b[42]),
            "allzero": not any(b[:202])}


# ---------- readers ----------

def read_jsonl(patterns):
    out = []
    for pat in patterns:
        for f in sorted(glob.glob(exp(pat))):
            try:
                with open(f) as fh:
                    for line in fh:
                        if line.startswith("{"):
                            try:
                                out.append(json.loads(line))
                            except ValueError:
                                pass
            except OSError:
                pass
    return out


def fgb_summary(recs):
    acc = [r for r in recs if r.get("accept") and r.get("hex") and "t_cw" in r]
    by_id = collections.defaultdict(list)
    for r in acc:
        by_id[r["hex"][:15]].append(r)
    beacons, bursts_t, relays = [], [], collections.Counter()
    only_combined = 0
    for hid, rs in by_id.items():
        rs.sort(key=lambda r: r["t_cw"])
        groups, cur = [], []
        for r in rs:
            if cur and r["t_cw"] - cur[0]["t_cw"] > 1.5:
                groups.append(cur)
                cur = []
            cur.append(r)
        groups.append(cur)
        if not any(r["accept"] == "single" for r in rs) and len(groups) < 2:
            continue
        hx = collections.Counter(r["hex"] for r in rs).most_common(1)[0][0]
        selftest = sum(1 for r in rs if r.get("selftest_sync")) > 0.8 * len(rs)
        cls, country = fgb_class(hx)
        for g in groups:
            n = max([len({round(r["f_hz"], -2) for r in g if "combined" not in r})] + [r.get("combined", 0) for r in g])
            relays[n] += 1
            only_combined += all("combined" in r for r in g)
            bursts_t.append(g[0].get("utc"))
        eb = [r["ebn0"] for r in rs if "combined" not in r] or [rs[0]["ebn0"]]
        beacons.append({"id": hid, "hex": hx, "class": cls, "selftest": selftest, "country": country, "bursts": len(groups),
                        "first": rs[0].get("utc"), "last": rs[-1].get("utc"), "ebn0_max": round(max(eb), 1),
                        "ebn0_med": round(sorted(eb)[len(eb) // 2], 1),
                        "relays_max": max(max(len({round(r['f_hz'], -2) for r in g}) for g in groups), 1)})
    return beacons, bursts_t, relays, only_combined, len(recs)


def sgb_summary(recs):
    out = []
    for r in recs:
        # SGB records only: the C receiver tags them; the Python scan's have relays and no t_cw
        if not (r.get("kind") == "sgb" or ("relays" in r and "t_cw" not in r)):
            continue
        if not r.get("ok") or not r.get("hex"):
            continue
        if r.get("kind") == "sgb" and not r.get("accept"):
            continue
        f = sgb_fields(r["hex"])
        if f["allzero"] or (MIDS and f["country"] not in MIDS):
            continue
        out.append({"utc": r.get("utc"), "t": r.get("t"), "band": r.get("band"), "hex": r["hex"], **f,
                    "relays": len(r.get("relays", [])) or 1, "how": r.get("how", ""),
                    "ebn0": [x.get("ebn0") for x in r.get("relays", [])],
                    "alone": any(x.get("alone") for x in r.get("relays", []))})
    # one entry per burst (several scans can cover the same second)
    uniq = {}
    for s in out:
        k = (round(s["t"] or 0, 1), s["hex"])
        uniq.setdefault(k, s)
    return list(uniq.values())


def recording(cfg):
    log = exp(cfg.get("record_log", "~/meosar-p0/record.log"))
    rec_dir = exp(cfg.get("record_dir", "~/meosar-p0/rec"))
    st = {"chunks": [], "waiting_until": None, "active": None, "log_tail": []}
    try:
        lines = [l.rstrip() for l in open(log) if "Streaming at" not in l]
    except OSError:
        lines = []
    st["log_tail"] = lines[-6:]
    for l in lines:
        m = re.search(r"waiting until (\S+)", l)
        if m:
            try:
                st["waiting_until"] = dt.datetime.fromisoformat(m.group(1)).astimezone(dt.timezone.utc).isoformat()
            except ValueError:
                st["waiting_until"] = m.group(1)
    for f in sorted(glob.glob(os.path.join(rec_dir, "*.sigmf-data"))):
        size = os.path.getsize(f)
        st["chunks"].append({"name": os.path.basename(f).replace(".sigmf-data", ""), "gb": round(size / 1e9, 2),
                             "minutes": round(size / 4 / 2.5e6 / 60, 1), "mtime": os.path.getmtime(f)})
    if st["chunks"] and time.time() - st["chunks"][-1]["mtime"] < 10:
        st["active"] = st["chunks"][-1]["name"]
    return st


def jobs():
    try:
        out = subprocess.run(["ps", "-eo", "etimes,pcpu,args"], capture_output=True, text=True, timeout=5).stdout
    except (OSError, subprocess.SubprocessError):
        return []
    pats = [("airspy_rx", "recording"), ("meosar-rx", "C receiver"), ("sgb_scan.py", "SGB search (Python)"),
            ("run.py", "FGB run (Python)"), ("scoreboard.sh", "scoreboard"), ("dec406_fgb_iq", "dec406 baseline"),
            ("fgbwin", "dec406 demod (scoreboard)"), ("iridium-sniffer", "iridium-sniffer")]
    res = []
    for line in out.splitlines()[1:]:
        parts = line.split(None, 2)
        if len(parts) < 3:
            continue
        et, cpu, args = parts
        argv = args.split()
        prog = os.path.basename(argv[0])
        if prog.startswith("python") and len(argv) > 1:
            prog = os.path.basename(argv[1])
        elif prog in ("bash", "sh") and len(argv) > 1 and not argv[1].startswith("-"):
            prog = os.path.basename(argv[1])
        for p, label in pats:
            if prog == p:
                m = re.search(r"--start (\S+)", args)
                res.append({"what": label, "secs": int(et), "cpu": float(cpu),
                            "detail": (f"from {m.group(1)} s" if m else "") + (" " + os.path.basename(args.split()[-1]) if p in ("meosar-rx", "run.py", "sgb_scan.py") else "")})
                break
    return res


def satellites(cfg):
    p = exp(cfg.get("passes", "~/meosar-p0/passes.json"))
    try:
        j = json.load(open(p))
    except (OSError, ValueError):
        return None
    t0 = dt.datetime.fromisoformat(j["start"])
    now = dt.datetime.now(dt.timezone.utc)
    m = int((now - t0).total_seconds() // 60 // 5 * 5)
    rows = j["rows"].get(str(m))
    if rows is None:
        return {"stale": True}
    sats = [{"con": r[0], "id": r[1], "name": r[2], "el": r[4], "az": r[5], "dop_hz": r[6]} for r in rows]
    sats.sort(key=lambda s: -s["el"])
    return {"slot_utc": (t0 + dt.timedelta(minutes=m)).isoformat(timespec="minutes"), "sats": sats}


SCORE_RE = re.compile(r"C/N0\s+([\d.]+)\s+ours: decoded (\d+)/(\d+).*?single ok (\d+) wrong (\d+)\s+weak ok (\d+) wrong (\d+)"
                      r"(?:.*?dec406: decoded (\d+)/\d+)?(?:.*?true onset: (\d+)/\d+)?")


def scoreboard(cfg):
    out = []
    for name, pat in cfg.get("scoreboards", {}).items():
        rows = []
        for f in sorted(glob.glob(exp(pat)), key=lambda s: float(re.findall(r"c(\d+)", os.path.basename(s))[0])):
            m = SCORE_RE.search(open(f).read())
            if m:
                g = m.groups()
                rows.append({"cn0": float(g[0]), "ours": int(g[1]), "n": int(g[2]), "wrong_single": int(g[4]),
                             "dec406": int(g[7]) if g[7] else None, "dec406_demod": int(g[8]) if g[8] else None})
        if rows:
            out.append({"name": name, "rows": rows})
    return out


def build_state(cfg):
    datasets = []
    all_beacons = {}
    timeline = collections.Counter()
    sgb_all = []
    for ds in cfg.get("datasets", []):
        fgb = read_jsonl(ds.get("fgb", []))
        sgb = read_jsonl(ds.get("sgb", []))
        beacons, bt, relays, only_c, ncand = fgb_summary(fgb)
        sg = sgb_summary(sgb)
        for t in bt:
            if t:
                timeline[t[:15] + ("0" if int(t[15]) < 5 else "5")] += 1   # 5-minute bins
        for b in beacons:
            a = all_beacons.setdefault(b["id"], {**b, "datasets": 0})
            if a is not b:
                a["bursts"] += b["bursts"]
                a["last"] = max(a["last"] or "", b["last"] or "")
                a["first"] = min(a["first"] or "~", b["first"] or "~")
                a["ebn0_max"] = max(a["ebn0_max"], b["ebn0_max"])
                a["relays_max"] = max(a["relays_max"], b["relays_max"])
            a["datasets"] += 1
        sgb_all += sg
        datasets.append({"name": ds["name"], "processed": bool(fgb), "candidates": ncand,
                         "beacons": len(beacons), "bursts": len(bt), "only_combined": only_c,
                         "relays": dict(sorted(relays.items())), "sgb": len(sg),
                         "note": ds.get("note", "")})
    public, hidden = [], collections.Counter()
    for b in all_beacons.values():
        if b["class"] == "operational":
            hidden["fgb"] += 1
        if b["class"] != "operational" or FULL:
            if FULL:
                b["fields"] = decode_fields(b["hex"], cfg)
            public.append(b)
    public.sort(key=lambda b: -b["bursts"])
    sgb_public, sgb_hidden = [], 0
    for s in sorted(sgb_all, key=lambda s: s["utc"] or ""):
        if not s["test"]:
            sgb_hidden += 1
        if s["test"] or FULL:
            if FULL:
                s["fields"] = decode_fields(s["hex"], cfg)
            sgb_public.append(s)
    try:
        progress = json.load(open(HERE / "progress.json"))
    except (OSError, ValueError):
        progress = []
    return {
        "generated": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "station": cfg.get("station", "ESG-KATL"),
        "datasets": datasets,
        "totals": {"beacons": len(all_beacons), "public": len(public), "operational": hidden["fgb"],
                   "bursts": sum(d["bursts"] for d in datasets), "sgb": len(sgb_all),
                   "only_combined": sum(d["only_combined"] for d in datasets),
                   "countries": len({b["country"] for b in all_beacons.values()})},
        "classes": dict(collections.Counter(b["class"] for b in all_beacons.values())),
        "beacons": public,
        "sgb": sgb_public, "sgb_operational": sgb_hidden,
        "timeline": sorted(timeline.items()),
        "recording": recording(cfg),
        "jobs": jobs(),
        "satellites": satellites(cfg),
        "scoreboards": scoreboard(cfg),
        "progress": progress,
        "full": FULL,
    }


def state(cfg):
    with _lock:
        if time.time() - _cache["t"] > CACHE_S or _cache["state"] is None:
            _cache["state"] = build_state(cfg)
            _cache["t"] = time.time()
        return _cache["state"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(HERE / "dash.json"))
    ap.add_argument("--port", type=int, default=8896)
    ap.add_argument("--host", default=None)
    ap.add_argument("--full", action="store_true",
                    help="every beacon with its decoded fields; refuses any host but loopback")
    a = ap.parse_args()
    global FULL
    FULL = a.full
    if a.host is None:
        a.host = "127.0.0.1" if FULL else "0.0.0.0"
    if FULL and a.host not in LOOPBACK:
        raise SystemExit("--full shows operational beacons and positions: it only listens on 127.0.0.1 or ::1")
    load_mids()
    page = (HERE / "page.html").read_bytes()

    class H(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def do_GET(self):
            if self.path.split("?")[0] in ("/", "/index.html"):
                body, ctype = page, "text/html; charset=utf-8"
            elif self.path.startswith("/api/state"):
                cfg = json.load(open(a.config))
                body, ctype = json.dumps(state(cfg)).encode(), "application/json"
            else:
                self.send_error(404)
                return
            self.send_response(200)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

    srv = ThreadingHTTPServer((a.host, a.port), H)
    print(f"meosar-rx dashboard on http://{a.host}:{a.port}/", flush=True)
    srv.serve_forever()


if __name__ == "__main__":
    main()
