#!/bin/bash
# Run the C receiver on each finished recording chunk matching PATTERN, once.
# usage: process_chunks.sh PATTERN UNTIL_EPOCH   (e.g. '*_pm*' "$(date -d '16:00' +%s)")
set -u
PAT=$1; UNTIL=$2
BIN=~/meosar-rx/c/build/meosar-rx; OUT=~/meosar-rx/out; LOG=~/meosar-rx/process_chunks.log
while [ "$(date +%s)" -lt "$UNTIL" ]; do
  for f in ~/meosar-p0/rec/$PAT.sigmf-data; do
    [ -e "$f" ] || continue
    name=$(basename "$f" .sigmf-data); out=$OUT/pm_$name.jsonl
    [ -e "$out" ] && continue
    # finished: full size, or not written for two minutes
    size=$(stat -c %s "$f"); age=$(( $(date +%s) - $(stat -c %Y "$f") ))
    [ "$size" -ge 18000000000 ] || [ "$age" -ge 120 ] || continue
    echo "$(date -Is) start $name" >> $LOG
    nice "$BIN" -i "${f%.sigmf-data}.sigmf-meta" --sgb-threads 4 > "$out.tmp" 2>> $LOG && mv "$out.tmp" "$out"
    echo "$(date -Is) done $name: $(grep -c '"accept": "single"' "$out" 2>/dev/null) single FGB decodes, $(grep -c '"kind": "sgb".*"accept": "single"' "$out" 2>/dev/null) SGB" >> $LOG
  done
  sleep 60
done
