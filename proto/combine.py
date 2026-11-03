"""Second pass over a run's records: combine repeats of the same beacon across bursts.

1. Known beacons: each relay record with soft bits (llr) is tested against the frames of the
   confirmed beacons (normalised correlation over bits 24..105, z-score); z >= Z_KNOWN marks it
   as a burst of that beacon, even when it is far too weak to decode alone.
2. Unknown beacons: relays not explained by a known beacon are clustered by the correlation of
   their soft bits with each other; each cluster's LLRs are summed and decoded. A member must
   also match the decoded frame on its own (z >= Z_MEMBER), so two similar beacons are not merged.

usage: combine.py run.jsonl [more.jsonl ...] > combined.jsonl   (summary on stderr)
"""
from __future__ import annotations

import collections
import json
import sys

import numpy as np

import fgb

Z_KNOWN, Z_PAIR, Z_MEMBER = 6.0, 5.0, 3.0
SEG = slice(24, 106)          # BCH-1 codeword: PDF-1 + parity; identical in every burst of a beacon


def zscore(L: np.ndarray, s: np.ndarray) -> float:
    """Correlation of soft values L with +-1 bits s, ~N(0,1) when unrelated."""
    den = np.sqrt(np.sum(L * L))
    return float(np.sum(L * s) / den) if den > 0 else 0.0


def main():
    recs = [json.loads(l) for f in sys.argv[1:] for l in open(f) if l.startswith('{"t"')]
    relays = [r for r in recs if "llr" in r]
    acc = [r for r in recs if r.get("accept")]
    # confirmed beacons and their most common BCH-1 codeword
    by_id = collections.defaultdict(list)
    for r in acc:
        by_id[r["hex"][:15]].append(r)
    known = {}
    for hid, rs in by_id.items():
        if any(r["accept"] == "single" for r in rs) or len({round(r["t_cw"]) for r in rs}) >= 2:
            hx = collections.Counter(r["hex"] for r in rs).most_common(1)[0][0]
            bits = [1] * 15 + list(fgb.SYNC_NORMAL) + [int(c) for c in bin(int(hx, 16))[2:].zfill(len(hx) * 4)]
            known[hid] = 2 * np.array(bits[SEG]) - 1
    K = list(known)
    S = np.array([known[k] for k in K]) if K else np.zeros((0, 82))
    L = np.array([np.array(r["llr"])[SEG] for r in relays])
    norms = np.sqrt((L * L).sum(1))
    Z = (L @ S.T) / norms[:, None] if K else np.zeros((len(relays), 0))
    per_beacon = collections.Counter()
    matched = np.zeros(len(relays), dtype=bool)
    for i, r in enumerate(relays):
        if K and Z[i].max() >= Z_KNOWN:
            k = int(Z[i].argmax())
            matched[i] = True
            per_beacon[K[k]] += 1
            print(json.dumps({"kind": "known", "t_cw": r["t_cw"], "f_hz": r["f_hz"], "ebn0": r["ebn0"],
                              "id": K[k], "z": round(float(Z[i, k]), 1), "decoded_alone": bool(r.get("accept"))}))
    # unknown: cluster unmatched relays that look like bursts
    U = [i for i in range(len(relays)) if not matched[i] and relays[i]["sync"] >= 0.8]
    new = []
    if len(U) >= 2:
        Lu = L[U] / norms[U][:, None]
        C = (Lu @ Lu.T) * np.sqrt(82.0)
        np.fill_diagonal(C, 0)
        seen = set()
        for a in range(len(U)):
            if a in seen:
                continue
            grp = {a} | {b for b in np.where(C[a] >= Z_PAIR)[0] if b not in seen}
            # one level of expansion
            for b in list(grp):
                grp |= {c for c in np.where(C[b] >= Z_PAIR)[0] if c not in seen}
            if len(grp) < 2:
                continue
            members = [relays[U[g]] for g in sorted(grp)]
            # drop duplicate detections of the same relay (same time and carrier)
            uniq = []
            for m in members:
                if not any(abs(m["t_cw"] - q["t_cw"]) < 0.002 and abs(m["f_hz"] - q["f_hz"]) < 20 for q in uniq):
                    uniq.append(m)
            if len(uniq) < 2:
                continue
            llr = np.sum([np.array(m["llr"]) for m in uniq], axis=0)
            d = fgb.decode_soft(llr, L=12)
            if not d.ok or d.bch1_flips > 3 or d.bch2_flips > 2:
                continue
            s = 2 * np.array(d.bits[SEG]) - 1
            z_m = [zscore(np.array(m["llr"])[SEG], s) for m in uniq]
            if min(z_m) < Z_MEMBER:
                continue
            seen |= grp
            n_bursts = len({round(m["t_cw"], 0) for m in uniq})
            rec = {"kind": "repeat", "id": d.hex[:15], "hex": d.hex, "members": len(uniq), "bursts": n_bursts,
                   "ebn0_each": [m["ebn0"] for m in uniq], "z_each": [round(z, 1) for z in z_m],
                   "flips": [d.bch1_flips, d.bch2_flips], "t_first": uniq[0].get("utc", uniq[0]["t_cw"]),
                   "new": d.hex[:15] not in known}
            new.append(rec)
            print(json.dumps(rec))
    n_alone = sum(1 for r in relays if r.get("accept"))
    print(f"relays {len(relays)}  decoded alone {n_alone}  matched to known beacons {int(matched.sum())}  "
          f"known beacons {len(K)}  repeat clusters decoded {len(new)} "
          f"(new IDs {len({r['id'] for r in new if r['new']})})", file=sys.stderr)
    for hid, n in per_beacon.most_common(12):
        print(f"  {hid}: {n} relays matched", file=sys.stderr)


if __name__ == "__main__":
    main()
