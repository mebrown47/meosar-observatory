"""Group accepted decodes into bursts and beacons.
A burst = decodes of one ID within 1.5 s (relays, duplicates, combined records).
A beacon ID counts when it has a 'single' decode, or 'weak' decodes in >= 2 different bursts.
usage: group.py out.jsonl [more.jsonl ...]"""
import json, sys, collections
recs = [json.loads(l) for f in sys.argv[1:] for l in open(f) if l.startswith('{"t"')]
acc = [r for r in recs if r.get("accept")]
by_id = collections.defaultdict(list)
for r in acc:
    by_id[r["hex"][:15]].append(r)
beacons, n_bursts, relays, combined_only = {}, 0, collections.Counter(), 0
for hid, rs in by_id.items():
    rs.sort(key=lambda r: r["t_cw"])
    bursts, cur = [], []
    for r in rs:
        if cur and r["t_cw"] - cur[0]["t_cw"] > 1.5:
            bursts.append(cur); cur = []
        cur.append(r)
    bursts.append(cur)
    if not any(r["accept"] == "single" for r in rs) and len(bursts) < 2:
        continue
    beacons[hid] = (rs, bursts)
    n_bursts += len(bursts)
    for b in bursts:
        relays[max([len({round(r['f_hz'], -2) for r in b if 'combined' not in r})] + [r.get('combined', 0) for r in b])] += 1
        if all("combined" in r for r in b):
            combined_only += 1
print(f"candidates {len(recs)}, accepted decodes {len(acc)}, beacons {len(beacons)}, bursts {n_bursts} "
      f"(of which decoded only after combining relays: {combined_only})")
print("relays per decoded burst:", dict(sorted(relays.items())))
for hid, (rs, bursts) in sorted(beacons.items(), key=lambda kv: -len(kv[1][1])):
    eb = sorted(r["ebn0"] for r in rs)
    tiers = collections.Counter(r["accept"] for r in rs)
    print(f"  {hid}  bursts {len(bursts):3d}  single {tiers['single']:3d} weak {tiers['weak']:3d}  Eb/N0 {eb[0]:5.1f}..{eb[-1]:4.1f} dB  first {rs[0].get('utc', rs[0]['t_cw'])}")
