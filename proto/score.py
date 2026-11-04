"""Score decodes against synthetic truth. usage: score.py BASE ours.jsonl [dec406.txt] [fgbwin.txt]"""
import json, re, sys
truth = json.load(open(sys.argv[1] + ".truth.json"))
T = truth["bursts"]
ours = [json.loads(l) for l in open(sys.argv[2]) if l.startswith('{"t"')]
def match(t, tol=0.08):
    return next((b for b in T if abs(b["t_cw"] - t) < tol), None)
res = {"single": [0, 0], "weak": [0, 0]}
got, full = set(), set()
for r in ours:
    if not r.get("accept") or "t_cw" not in r: continue
    b = match(r["t_cw"])
    good = b is not None and r["hex"][:15] == b["hex"][:15]      # beacon ID (BCH-1 protected)
    res[r["accept"]][0 if good else 1] += 1
    if good:
        got.add(b["t_cw"])
        if r["hex"] == b["hex"]: full.add(b["t_cw"])
n = len(T)
line = f"C/N0 {truth['cn0']:4.1f}  ours: decoded {len(got)}/{n} (full frame {len(full)})  single ok {res['single'][0]} wrong {res['single'][1]}  weak ok {res['weak'][0]} wrong {res['weak'][1]}"
if len(sys.argv) > 3:
    txt = open(sys.argv[3]).read()
    # dec406_fgb_iq: "[t s] BURST ..." followed by a decode or a failure line
    blocks = re.split(r"\n(?=\[\s*[\d.]+ s\] BURST)", txt)
    ok = wrong = 0; seen = set()
    for blk in blocks:
        m = re.match(r"\[\s*([\d.]+) s\] BURST", blk)
        if not m or "BEACON DECODE" not in blk: continue
        b = match(float(m.group(1)), 0.25)
        c = re.search(r"Country: (\d+)", blk)
        if b and c and int(c.group(1)) == b["country"] and b["t_cw"] not in seen:
            ok += 1; seen.add(b["t_cw"])
        else: wrong += 1
    line += f"   dec406: decoded {ok}/{n} wrong {wrong}"
if len(sys.argv) > 4:
    ok = sum(1 for l in open(sys.argv[4]) if l.startswith("== BCH OK"))
    line += f"   dec406 demod given true onset: {ok}/{n}"
print(line)
