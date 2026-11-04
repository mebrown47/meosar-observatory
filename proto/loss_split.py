"""Split the remaining loss at threshold: demodulator BER (vs oracle) and which decoding stage
rejects each frame, using the C receiver's records (with llr) for a synthetic file."""
import json, sys
from collections import Counter
import numpy as np
import fgb
base, recs = sys.argv[1], sys.argv[2]
T = json.load(open(base + ".truth.json"))
R = [json.loads(l) for l in open(recs) if l.startswith('{"t"') and '"llr"' in l]
stage = Counter(); bers = []
for b in T["bursts"]:
    bits = np.array([1] * 15 + list(fgb.SYNC_NORMAL) + [int(c) for c in bin(int(b["hex"], 16))[2:].zfill(120)])
    m = [r for r in R if abs(r["t_cw"] - b["t_cw"]) < 0.03]
    if not m:
        stage["not detected"] += 1; continue
    r = max(m, key=lambda r: r["sync"]); L = np.array(r["llr"]); hard = (L > 0).astype(int)
    bers.append(np.mean(hard[24:] != bits[24:]))
    ds = int(np.sum(hard[15:24] != bits[15:24]))
    c1, cost1, f1, _ = fgb.BCH1.chase(L[24:106], L=10)
    ok1 = c1 is not None and list(c1) == list(bits[24:106])
    c2, cost2, f2, _ = fgb.BCH2.chase(L[106:144], L=8)
    ok2 = c2 is not None and list(c2) == list(bits[106:144])
    if not ok1: stage["BCH-1 wrong/none"] += 1
    elif not ok2: stage["BCH-1 right, BCH-2 wrong/none"] += 1
    elif ds > 1: stage["BCH right, rejected by hard sync check (>1 of 9 sync bits)"] += 1
    elif r["ebn0"] < 0 or r["sync"] < 1.0: stage["right, but below the weak gate (Eb/N0 est < 0 or sync < 1)"] += 1
    else: stage["right and passes"] += 1
print(f"C/N0 {T['cn0']}: detected {40 - stage['not detected']}/40, BER median {np.median(bers):.3f} "
      f"(theory {0.5 * __import__('math').erfc((10 ** ((T['cn0'] - 27.02) / 10)) ** 0.5):.3f})")
for k, v in stage.most_common(): print(f"  {v:2d}  {k}")
