#!/bin/bash
# Synthetic scoreboard: this receiver vs dec406 at several C/N0, on bursts injected into real noise.
# usage: scoreboard.sh SRC.sigmf-data   (a recording of the downlink; its quiet band is the noise)
# env: D (output dir, default ./syn), LEVELS (C/N0 list), BIN (the receiver),
#      DEC406 (dec406_fgb_iq, skipped if missing), FGBWIN (a windowed decoder, optional)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC=$1; D=${D:-./syn}; mkdir -p "$D"; LEVELS=${LEVELS:-"22 24 26 28 30 32 34 36 38 40"}
BIN=${BIN:-$HERE/../build/meosar-rx}; DEC406=${DEC406:-$HOME/dec406/build/dec406_fgb_iq}; FGBWIN=${FGBWIN:-}
for c in $LEVELS; do python3 "$HERE/synth.py" "$SRC" "$D/c$c" --cn0 "$c" --n 40 --spacing 1.6 --src-start 100 --fixed-f; done
for c in $LEVELS; do
  (
    nice "$BIN" -i "$D/c$c.cs16" -F cs16 --bands quiet --no-sgb --no-llr > "$D/c$c.ours.jsonl" 2>/dev/null
    : > "$D/c$c.dec406.txt"; : > "$D/c$c.fgbwin.txt"
    [ -x "$DEC406" ] && timeout 3000 nice "$DEC406" "$D/c$c.cs16" -s 2500000 -i -f -750000 > "$D/c$c.dec406.txt" 2>&1
    if [ -n "$FGBWIN" ]; then
      python3 -c "import json;[print(b['t_cw'], b['f_hz']-1544.5e6) for b in json.load(open('$D/c$c.truth.json'))['bursts']]" |
        while read t f; do nice "$FGBWIN" "$D/c$c.cs16" 2500000 $t $f 2>/dev/null | grep -E "^== BCH OK|^no decode" >> "$D/c$c.fgbwin.txt"; done
    fi
    python3 "$HERE/score.py" "$D/c$c" "$D/c$c.ours.jsonl" "$D/c$c.dec406.txt" "$D/c$c.fgbwin.txt" > "$D/c$c.score" ) &
  [ $(jobs -r | wc -l) -ge 3 ] && wait -n
done
wait
cat "$D"/c*.score
