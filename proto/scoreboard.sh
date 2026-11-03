#!/bin/bash
# Synthetic scoreboard on the feeder: ours vs dec406 at several C/N0.
set -u
cd ~/meosar-rx/proto
SRC=$(ls ~/meosar-p0/rec/*am1.sigmf-data); D=${D:-~/meosar-rx/syn2cal}; mkdir -p $D; LEVELS=${LEVELS:-"22 24 26 28 30 32 34 36 38 40"}
for c in $LEVELS; do python3 synth.py $SRC $D/c$c --cn0 $c --n 40 --spacing 1.6 --src-start 100 --fixed-f; done
for c in $LEVELS; do
  (
    nice ~/meosar-rx/c/build/meosar-rx -i $D/c$c.cs16 -F cs16 --bands quiet --no-sgb --no-llr > $D/c$c.ours.jsonl 2>/dev/null
    timeout 3000 nice ~/dec406/build/dec406_fgb_iq $D/c$c.cs16 -s 2500000 -i -f -750000 > $D/c$c.dec406.txt 2>&1
    : > $D/c$c.fgbwin.txt
    python3 -c "import json;[print(b['t_cw'], b['f_hz']-1544.5e6) for b in json.load(open('$D/c$c.truth.json'))['bursts']]" |
      while read t f; do nice ../../meosar-p0/tools/fgbwin $D/c$c.cs16 2500000 $t $f 2>/dev/null | grep -E "^== BCH OK|^no decode" >> $D/c$c.fgbwin.txt; done
    python3 score.py $D/c$c $D/c$c.ours.jsonl $D/c$c.dec406.txt $D/c$c.fgbwin.txt > $D/c$c.score ) &
  [ $(jobs -r | wc -l) -ge 3 ] && wait -n
done
wait
cat $D/c*.score
