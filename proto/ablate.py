"""Ablation on one run: how much each receiver stage contributes, from the same records.

usage: ablate.py run.jsonl [...] [--dec406 fgbwin_results.jsonl]
Configurations (all use the same confirmation rule: a 'single' decode, or decodes of the same
ID in >= 2 different bursts):
  full           relays + relay combining + repeat pass (combine.py)
  no repeat      relays + relay combining
  no relay comb  relays only
  hard decoding  relays only, BCH decoded from hard decisions (Chase off, L=0) on the same soft bits
  dec406 demod   our detector, dec406's FGB demodulator + BCH on each burst-like candidate
"""
from __future__ import annotations

import collections
import json
import subprocess
import sys

import numpy as np

import fgb
from run import accept


def confirmed(decodes):
    """decodes: list of (id15, t_cw, tier). Returns (n_beacons, n_bursts, set of ids)."""
    by = collections.defaultdict(list)
    for hid, t, tier in decodes:
        by[hid].append((t, tier))
    ids, nb = set(), 0
    for hid, xs in by.items():
        xs.sort()
        bursts, last = 0, None
        for t, _ in xs:
            if last is None or t - last > 1.5:
                bursts += 1
                last = t
        if any(tier == "single" for _, tier in xs) or bursts >= 2:
            ids.add(hid)
            nb += bursts
    return len(ids), nb, ids


def main():
    args = sys.argv[1:]
    dec406_file = None
    if "--dec406" in args:
        i = args.index("--dec406")
        dec406_file = args[i + 1]
        args = args[:i] + args[i + 2:]
    recs = [json.loads(l) for f in args for l in open(f) if l.startswith('{"t"')]
    relays = [r for r in recs if "llr" in r]
    combined = [r for r in recs if "combined" in r]
    rows = []

    d_rel = [(r["hex"][:15], r["t_cw"], r["accept"]) for r in relays if r.get("accept")]
    d_comb = [(r["hex"][:15], r["t_cw"], r["accept"]) for r in combined if r.get("accept")]
    no_relay = confirmed(d_rel)
    no_rep = confirmed(d_rel + d_comb)

    # repeat pass: reuse combine.py and read its new IDs
    out = subprocess.run([sys.executable, "combine.py", *args], capture_output=True, text=True).stdout
    rep = [json.loads(l) for l in out.splitlines() if '"kind": "repeat"' in l]
    new_ids = {r["id"] for r in rep} - no_rep[2]
    rep_bursts = sum(r["bursts"] for r in rep if r["id"] in new_ids)
    full = (no_rep[0] + len(new_ids), no_rep[1] + rep_bursts)

    # hard decisions on the same soft bits
    d_hard = []
    for r in relays:
        d = fgb.decode_soft(np.array(r["llr"]), L=0)
        res = {"ok": d.ok, "flips": [d.bch1_flips, d.bch2_flips], "ebn0": r["ebn0"], "sync": r["sync"]}
        a = accept(res)
        if a:
            d_hard.append((d.hex[:15], r["t_cw"], a))
    hard = confirmed(d_hard)

    rows = [("full chain", full[0], full[1]),
            ("no repeat pass", no_rep[0], no_rep[1]),
            ("no relay combining", no_relay[0], no_relay[1]),
            ("hard-decision BCH (no Chase)", hard[0], hard[1])]
    if dec406_file:
        dd = [json.loads(l) for l in open(dec406_file)]
        d406 = [(fgb.frame_hex([int(c) for c in x["bits"]])[:15], x["t_cw"], "single") for x in dd if x.get("bits")]
        r406 = confirmed(d406)
        rows.append((f"dec406 demod + BCH on our {len(dd)} burst-like candidates", r406[0], r406[1]))
    print(f"relay records {len(relays)}, combined records {len(combined)}")
    print(f"{'configuration':58s} {'beacons':>8s} {'bursts':>8s}")
    for name, b, n in rows:
        print(f"{name:58s} {b:8d} {n:8d}")


if __name__ == "__main__":
    main()
