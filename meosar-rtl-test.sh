#!/usr/bin/env bash
# meosar-rx test run on an RTL-SDR dongle (through rtl_tcp).
#
# Builds meosar-rx from GitHub, checks the dongle, records for a few hours, then prints a summary and
# packs the results into one .tar.gz to send back. Run it as your normal user (it uses sudo only to
# install missing packages).
#
#   bash meosar-rtl-test.sh                 # 2 hours, gain 40 dB, bias-tee on
#   bash meosar-rtl-test.sh --hours 4 --gain 30 --no-bias
#
# Options:
#   --hours N    how long to record (default 2; fractions are fine, e.g. 0.5)
#   --gain DB    dongle gain in dB (default 40; try ~30 with a strong LNA in front)
#   --no-bias    don't power the LNA through the coax (bias-tee is on by default: an L-band antenna
#                needs an LNA; use this for a passive antenna, one with a DC short, or an LNA with
#                its own supply)
#   --no-sgb     first-generation beacons only: second-generation (SGB) acquisition is most of the
#                CPU load (about 85%); for a slow computer. On a Raspberry Pi it's the default.
#   --sgb        second-generation beacons too, even on a Raspberry Pi (it will fall behind)
#   --dir PATH   where to clone and build (default ~/meosar-rx-test)
#   --no-deps    don't install packages
#
# Needs: an RTL-SDR dongle (an R820T/R828D one tunes to 1544.5 MHz; a TCXO one such as the RTL-SDR
# Blog V3/V4 is best), an L-band antenna with a view of the sky and an LNA near the antenna.
# Nothing else may be using the dongle while it runs.

set -uo pipefail

HOURS=2; GAIN=40; BIAS=1; SGB=auto; DIR="$HOME/meosar-rx-test"; DEPS=1
REPO=https://github.com/mebrown47/meosar-observatory.git
RATE=2560000          # clean for the RTL's 28.8 MHz crystal; 2.5 MS/s drops samples
PORT=1234

while [ $# -gt 0 ]; do
  case "$1" in
    --hours) HOURS="$2"; shift 2;;
    --gain)  GAIN="$2"; shift 2;;
    --bias)  BIAS=1; shift;;
    --no-bias) BIAS=0; shift;;
    --no-sgb) SGB=0; shift;;
    --sgb)    SGB=1; shift;;
    --dir)   DIR="$2"; shift 2;;
    --no-deps) DEPS=0; shift;;
    -h|--help) sed -n '2,27p' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "unknown option $1 (--help for the list)"; exit 1;;
  esac
done

say()  { printf '\n==> %s\n' "$*"; }
ok()   { printf '  ok   %s\n' "$*"; }
warn() { printf '  warn %s\n' "$*"; }
die()  { printf '  error %s\n' "$*"; exit 1; }

SECS=$(awk -v h="$HOURS" 'BEGIN { s = int(h * 3600); if (s < 60) exit 1; print s }') \
  || die "--hours must be at least 0.02 (a minute)"

# A Raspberry Pi keeps up with first-generation beacons but not with second-generation acquisition too
# (measured: 0.95x real time without it, 0.16x with it, on a 4-core Pi with Debian 13).
MODEL=$(tr -d '\0' 2>/dev/null < /proc/device-tree/model || true)
if [ "$SGB" = auto ]; then
  case "$MODEL" in
    *"Raspberry Pi"*)
      SGB=0
      say "$MODEL"
      warn "second-generation (SGB) beacons off: their acquisition is most of the CPU load and a Pi"
      warn "falls far behind with it; first-generation beacons only (--sgb turns SGB on anyway)";;
    *) SGB=1;;
  esac
fi

# ---- packages ---------------------------------------------------------------
say "Packages"
need=()
command -v cmake >/dev/null      || need+=(cmake)
command -v cc >/dev/null         || need+=(build-essential)
command -v pkg-config >/dev/null || need+=(pkg-config)
command -v git >/dev/null        || need+=(git)
command -v python3 >/dev/null    || need+=(python3)
command -v rtl_tcp >/dev/null && command -v rtl_test >/dev/null || need+=(rtl-sdr)
pkg-config --exists fftw3f 2>/dev/null || need+=(libfftw3-dev)
if [ ${#need[@]} -eq 0 ]; then
  ok "everything needed is installed"
elif [ "$DEPS" -eq 0 ]; then
  die "missing: ${need[*]} (install them, or run without --no-deps)"
elif command -v apt-get >/dev/null; then
  sudo apt-get install -y -q --no-install-recommends "${need[@]}" || die "apt couldn't install: ${need[*]}"
  ok "installed ${need[*]}"
else
  die "not a Debian/Ubuntu system: install ${need[*]} (FFTW with single precision, rtl-sdr tools) and re-run"
fi

# ---- source and build -------------------------------------------------------
say "Building meosar-rx"
if [ -d "$DIR/.git" ]; then
  git -C "$DIR" pull --ff-only -q || die "git pull failed in $DIR (delete it and re-run)"
else
  git clone -q "$REPO" "$DIR" || die "git clone failed"
fi
cmake -S "$DIR" -B "$DIR/build" -DCMAKE_BUILD_TYPE=Release > "$DIR/build.log" 2>&1 \
  && cmake --build "$DIR/build" -j"$(nproc)" >> "$DIR/build.log" 2>&1 \
  || { tail -20 "$DIR/build.log"; die "build failed (full log: $DIR/build.log)"; }
BIN="$DIR/build/meosar-rx"
ok "$("$BIN" --version), commit $(git -C "$DIR" rev-parse --short HEAD)"

# ---- results folder ---------------------------------------------------------
STAMP=$(date -u +%Y%m%d-%H%M%SZ)
OUT="$HOME/meosar-rtl-$STAMP"
mkdir -p "$OUT"
{
  echo "meosar-rtl-test $STAMP"
  echo "receiver: $("$BIN" --version), commit $(git -C "$DIR" rev-parse --short HEAD)"
  echo "options: hours=$HOURS gain=$GAIN bias=$BIAS sgb=$SGB rate=$RATE"
  echo "system: $(uname -srm); $( (. /etc/os-release && echo "$PRETTY_NAME") 2>/dev/null)"
  cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ *//')
  echo "cpu: ${MODEL:-${cpu:-unknown}}, $(nproc) threads"
} > "$OUT/info.txt"

# ---- dongle -----------------------------------------------------------------
say "Dongle"
RTL_PID=""; RX_PID=""; RECORDING=0
cleanup() {
  [ "$RECORDING" -eq 0 ] && [ -n "${OUT:-}" ] && rm -rf "$OUT"   # stopped before recording: nothing to keep
  # TERM to timeout, which passes the receiver its INT: it finishes the records it has and exits.
  [ -n "$RX_PID" ] && kill -TERM "$RX_PID" 2>/dev/null && wait "$RX_PID" 2>/dev/null
  if [ -n "$RTL_PID" ]; then
    kill "$RTL_PID" 2>/dev/null; sleep 1
    kill -9 "$RTL_PID" 2>/dev/null   # rtl_tcp often ignores a plain stop and keeps the dongle busy
  fi
}
trap cleanup EXIT
trap 'echo; echo "  stopping early"; exit 130' INT TERM

# The kernel's DVB-T (TV) driver grabs an RTL dongle when it's plugged in, and then rtl_tcp can't
# open it. Blacklisted for good (the rtl-sdr packages usually do this; a build from source doesn't),
# then unloaded now: the driver first, then its helper modules, which hold it in use.
DVB_MODS=(dvb_usb_rtl28xxu rtl2832_sdr rtl2832 rtl2830)
BLACKLIST=/etc/modprobe.d/blacklist-rtl-sdr-dvb.conf
if lsmod 2>/dev/null | grep -q '^dvb_usb_rtl28xxu'; then
  warn "the kernel's DVB-T (TV) driver has the dongle; switching it off (sudo)"
  if ! grep -rqs '^[[:space:]]*blacklist[[:space:]]\+dvb_usb_rtl28xxu' /etc/modprobe.d/; then
    { echo "# Written by meosar-rtl-test.sh: keep the DVB-T driver off RTL-SDR dongles"
      printf 'blacklist %s\n' "${DVB_MODS[@]}"; } | sudo tee "$BLACKLIST" > /dev/null \
      || die "couldn't write $BLACKLIST"
    ok "blacklisted ${DVB_MODS[*]} in $BLACKLIST (stays off after a replug or reboot)"
  fi
  for m in "${DVB_MODS[@]}"; do
    lsmod | grep -q "^$m " || continue
    err=$(sudo modprobe -r "$m" 2>&1) || warn "couldn't unload $m: $err"
  done
  if lsmod | grep -q '^dvb_usb_rtl28xxu'; then
    die "the DVB-T driver is still loaded: unplug the dongle, plug it back in (the blacklist keeps the driver off now) and re-run; or reboot"
  fi
  ok "DVB-T driver unloaded"
fi
if pgrep -x rtl_tcp >/dev/null || pgrep -x rtl_sdr >/dev/null; then
  die "rtl_tcp or rtl_sdr is already running (pkill -9 rtl_tcp), or stop the program using the dongle"
fi

# Sample loss at our rate, 20 s: the dongle and USB have to keep up before anything else matters.
timeout -s INT 20 rtl_test -s "$RATE" > "$OUT/rtl_test.txt" 2>&1
grep -q 'Found .* tuner' "$OUT/rtl_test.txt" \
  || { cat "$OUT/rtl_test.txt"; die "rtl_test couldn't open the dongle (in use by another program?)"; }
tuner=$(grep -o 'Found .* tuner' "$OUT/rtl_test.txt" | sed 's/^Found //')
lost=$(grep -o 'million lost (minimum): [0-9]*' "$OUT/rtl_test.txt" | grep -o '[0-9]*$')
echo "dongle: $tuner, rtl_test at $RATE S/s: ${lost:-?} per million lost" >> "$OUT/info.txt"
ok "$tuner"
case "$tuner" in *R820T*|*R828D*) ;; *) warn "this tuner may not reach 1544.5 MHz (R820T/R828D do)";; esac
if [ -n "$lost" ] && [ "$lost" -gt 10 ]; then
  warn "rtl_test lost $lost samples per million at $RATE S/s: try another USB port (not a hub)"
else
  ok "no samples lost at $RATE S/s in 20 s"
fi
TCP_ARGS=(-a 127.0.0.1 -p "$PORT")
if [ "$BIAS" -eq 1 ]; then
  # rtl_tcp -h exits non-zero: read its help text first (pipefail would fail the grep otherwise).
  tcp_help=$(rtl_tcp -h 2>&1 || true)
  case "$tcp_help" in *"-T "*|*"[-T"*) TCP_ARGS+=(-T);; esac
  ok "bias-tee ON: the dongle will power the LNA through the coax (--no-bias for a passive antenna)"
else
  ok "bias-tee off (--no-bias): the LNA needs its own supply, or the antenna is passive"
fi
echo "bias-tee: $([ "$BIAS" -eq 1 ] && echo "on, rtl_tcp ${TCP_ARGS[*]}" || echo off)" >> "$OUT/info.txt"
sleep 1

# ---- record -----------------------------------------------------------------
say "Recording for $HOURS h (until $(date -d "+$SECS sec" '+%H:%M %Z')); Ctrl-C stops early and keeps what's recorded"
rtl_tcp "${TCP_ARGS[@]}" > "$OUT/rtl_tcp.log" 2>&1 &
RTL_PID=$!
for _ in $(seq 20); do
  ss -ltn 2>/dev/null | grep -q "127.0.0.1:$PORT " && break
  kill -0 "$RTL_PID" 2>/dev/null || { cat "$OUT/rtl_tcp.log"; die "rtl_tcp stopped"; }
  sleep 0.5
done
args=(-i "rtl_tcp:127.0.0.1:$PORT" -s "$RATE" -g "$GAIN" --stats 300)
[ "$BIAS" -eq 1 ] && args+=(-b)
[ "$SGB" -eq 0 ] && args+=(--no-sgb)
echo "command: meosar-rx ${args[*]}" >> "$OUT/info.txt"
timeout -s INT "$SECS" "$BIN" "${args[@]}" > "$OUT/records.jsonl" 2> "$OUT/receiver.log" &
RX_PID=$!
RECORDING=1
T0=$(date +%s)
if [ "$BIAS" -eq 1 ]; then
  # rtl_tcp logs it when it sets the bias-tee (from -T, or the receiver's command).
  for _ in $(seq 10); do grep -qi 'bias.*1\|bias.*on' "$OUT/rtl_tcp.log" && break; sleep 0.5; done
  if grep -qi 'bias.*1\|bias.*on' "$OUT/rtl_tcp.log"; then
    ok "rtl_tcp switched the bias-tee on"
  else
    warn "rtl_tcp didn't report switching the bias-tee on: this rtl-sdr may be too old to power the"
    warn "LNA (rtl-sdr 2.0 or the RTL-SDR Blog drivers can); otherwise give the LNA its own supply"
  fi
fi
# Ctrl-C from here on stops the recording early and still gives the summary.
trap 'echo; echo "  stopping early"; kill -TERM "$RX_PID" 2>/dev/null' INT TERM
# Show the receiver's status lines as they come (every 10 minutes).
tail -n +1 -f "$OUT/receiver.log" --pid="$RX_PID" 2>/dev/null | sed -u 's/^/  /' &
while kill -0 "$RX_PID" 2>/dev/null; do wait "$RX_PID"; done
RX_PID=""
T1=$(date +%s)
echo "recorded: $(( (T1 - T0) / 60 )) min" >> "$OUT/info.txt"
trap 'echo; echo "  stopping early"; exit 130' INT TERM
sleep 1

# ---- summary ----------------------------------------------------------------
say "Summary"
python3 - "$OUT/records.jsonl" "$((T1 - T0))" <<'PY' | tee "$OUT/summary.txt"
import json, sys
recs = []
for line in open(sys.argv[1]):
    try: recs.append(json.loads(line))
    except ValueError: pass
stats = [r for r in recs if r.get("kind") in ("start", "stat")]
fgb = [r for r in recs if "t_cw" in r and r.get("accept") and r.get("hex")]
rep = [r for r in recs if r.get("kind") in ("fgb_match", "fgb_repeat")]
sgb = [r for r in recs if r.get("kind") == "sgb" and r.get("ok") and r.get("accept") and r.get("hex")]
tried_sgb = [r for r in recs if r.get("kind") == "sgb"]
cand = [r for r in recs if "t_cw" in r and "combined" not in r]
last = stats[-1] if stats else {}
# The receiver writes a status record every 5 minutes (not when it is stopped): the totals are
# up to the last one.
print(f"  recorded {int(sys.argv[2]) / 3600:.2f} h; up to the last status, {last.get('dropped_s', 0):.2f} s "
      f"of samples dropped, {last.get('cpu', 0):.2f} CPU cores")
print(f"  FGB: {len(cand)} candidates, {len(fgb)} bursts decoded, {len(rep)} more from repeat combining")
print(f"  SGB: {len(tried_sgb)} attempts, {len(sgb)} decoded")
def ebn0(r):   # SGB keeps it per relay
    return r["ebn0"] if "ebn0" in r else max((q.get("ebn0", -99) for q in r.get("relays", [])), default=-99)
def beacon(h):   # the beacon ID the MEOSAR Observatory uses
    if len(h) <= 40:
        return "FGB " + h[:15]
    b = bin(int(h, 16))[2:].zfill(252)[2:]
    return f"SGB {int(b[30:40], 2)}-{int(b[0:16], 2)}-{int(b[16:30], 2)}"   # country-TAC-serial

heard = {}
for r in fgb + sgb:
    heard.setdefault(beacon(r["hex"]), []).append(r)
# A strict decode stands alone; a weak one (soft decoding with several bit flips) can come from noise,
# so a beacon heard only weakly counts once it is decoded in a second burst.
beacons = {k: rs for k, rs in heard.items() if len(rs) > 1 or any(r.get("accept") == "single" for r in rs)}
lone = len(heard) - len(beacons)
if beacons:
    best = max(ebn0(r) for rs in beacons.values() for r in rs)
    print(f"  {len(beacons)} beacon(s) confirmed; strongest Eb/N0 {best:.1f} dB")
    ranked = sorted(beacons.items(), key=lambda kv: -len(kv[1]))
    for bid, rs in ranked[:12]:
        print(f"    {bid:<22} {len(rs):4d} burst(s), best Eb/N0 {max(ebn0(r) for r in rs):5.1f} dB")
    if len(ranked) > 12:
        print(f"    ... and {len(ranked) - 12} more")
else:
    print("  no beacon confirmed")
if lone:
    print(f"  {lone} beacon ID(s) from a single weak decode: likely noise, not counted")
if not beacons:
    print("  if this was with an antenna: check its view of the sky and that the LNA is powered")
PY

tar -czf "$OUT.tar.gz" -C "$(dirname "$OUT")" "$(basename "$OUT")"
say "Done"
echo "  results: $OUT"
echo "  send this file back: $OUT.tar.gz ($(du -h "$OUT.tar.gz" | cut -f1))"
