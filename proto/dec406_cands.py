"""Run dec406's FGB demodulator + BCH (via fgbwin, narrow grid around our estimates) on every
burst-like candidate of a run. usage: dec406_cands.py FILE.sigmf-data run.jsonl [...] > out.jsonl"""
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

from run import meta_of
from pathlib import Path

FGBWIN = os.path.expanduser("~/meosar-p0/tools/fgbwin")


def main():
    data = sys.argv[1]
    fs, fc, _ = meta_of(Path(data))
    recs = [json.loads(l) for f in sys.argv[2:] for l in open(f) if l.startswith('{"t"')]
    cands = []
    for r in recs:
        if "llr" not in r or r["sync"] < 0.8:
            continue
        if any(abs(r["t_cw"] - c["t_cw"]) < 0.002 and abs(r["f_hz"] - c["f_hz"]) < 20 for c in cands[-20:]):
            continue
        cands.append(r)
    env = dict(os.environ, FAST="1", QUIET="1")

    def one(r):
        p = subprocess.run([FGBWIN, data, str(int(fs)), f"{r['t_cw']:.4f}", f"{r['f_hz'] - fc:.1f}"],
                           capture_output=True, text=True, env=env)
        bits = ""
        for line in p.stdout.splitlines():
            if line.startswith("== BCH OK"):
                bits = line.rsplit(" ", 1)[-1]
        return {"t_cw": r["t_cw"], "f_hz": r["f_hz"], "bits": bits}

    with ThreadPoolExecutor(int(os.environ.get("JOBS", "3"))) as ex:
        for res in ex.map(one, cands):
            print(json.dumps(res), flush=True)


if __name__ == "__main__":
    main()
