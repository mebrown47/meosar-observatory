# meosar-rx

A receiver for Cospas-Sarsat 406 MHz distress-beacon bursts relayed on the MEOSAR L-band
downlink (Galileo SAR 1544.05–1544.15 MHz, GLONASS-K 1544.85–1544.95 MHz), built for weak
signals on a small fixed antenna. It decodes both beacon generations, combines the copies of a
burst relayed by different satellites, and comes with a station installer and a web dashboard,
the MEOSAR Observatory: beacons, bursts, which satellite relayed each copy, and beacon positions
estimated from the relays where the beacons don't report one.

**This is not part of Cospas-Sarsat.** Alerts reach rescue services through the system's own
ground stations, not through this program. Don't use it for alerting, and check that receiving
and decoding these signals is allowed where you are.

Status: running continuously on one station since November 2026, on an Airspy. The HackRF,
rtl_tcp and spyserver inputs come from inmarsat-rx, where they are in use, but haven't been tested
with meosar-rx yet. `proto/` holds the Python prototype used to design and measure the receiver
(see `proto/README.md`); a station doesn't need it.

## What you need

- An SDR that covers 1544.5 MHz at 2.5 MS/s (the default; both SAR bands fit in it). Airspy is
  tested; HackRF, rtl_tcp and spyserver sources are supported.
- An L-band antenna with a view of the sky and a low-noise amplifier. `-b` powers an LNA over the
  coax (bias-tee).
- Linux with systemd for the station install (Debian or Ubuntu: the installer fetches the
  packages). The receiver alone builds anywhere with CMake, a C compiler and FFTW3.

## FGB chain

1. Coherent detector: 131 ms FFTs find the 160 ms unmodulated carrier at the start of each burst.
2. Coherent demodulator: carrier frequency and drift over the whole burst, timing from the 24
   known bits, phase smoothed forwards and backwards, soft bits.
3. Chase-II soft decoding of BCH-1 and BCH-2.
4. Relay combining: copies of one burst relayed by different satellites are summed before decoding.
5. Acceptance: a strong single decode, or the same beacon in two different bursts.
6. Repeat combining: bursts too weak to decode alone are matched to known beacons by correlating
   their soft bits with each beacon's BCH-1 codeword, or summed with earlier weak bursts that
   correlate with them and decoded together (`fgb_match` and `fgb_repeat` records). On the 30-minute
   Phase 0 recording at 13:17 UTC this adds 40% more bursts and two beacons that never decoded alone.

## SGB chain

1. Acquisition: the 166.7 ms spread-spectrum preamble, split into 25 one-bit segments, is searched
   at every lag and over the whole band by FFT correlation, the segments added coherently across
   fine frequency (PMF-FFT), in threads (`--sgb-threads`).
2. Coherent despreading: timing and frequency refined on the preamble, per-bit correlations on I
   and Q, carrier drift fitted over the whole burst, soft bits.
3. Relay combining: copies of one burst from different satellites are summed before decoding.
4. Soft decoding of the BCH(250,202) code.

## Build and run

```bash
sudo apt install build-essential cmake pkg-config libfftw3-dev libairspy-dev libhackrf-dev
cmake -S . -B build && cmake --build build -j
tests/run_tests.sh build
./build/meosar-rx -i recording.sigmf-meta > records.jsonl
./build/meosar-rx -i airspy -g 18 -b --archive archive --stats 600     # live
```

The build optimises for the CPU it runs on; `-DNATIVE=OFF` gives a portable binary. Input is a
SigMF recording (or a raw file with `-F cu8|cs8|cs16|cf32`), or live: `airspy[:SERIAL]`,
`hackrf[:SERIAL]`, `rtl_tcp:HOST:PORT`, `spyserver:HOST:PORT`. Output is one JSON record per burst
candidate on stdout (and over UDP with `-u HOST:PORT`). With `--archive DIR` records are also
appended to one file per UTC day, and `--stats 600` adds a status record every 10 minutes.
`meosar-rx` with no arguments lists every option.

## Beacon fields and positions: dec406 (optional)

The observatory decodes each frame's fields (protocol, country, identification, the position the
beacon reports) with the hex decoder of [dec406](https://github.com/moricef/dec406). Without it the
dashboard shows the beacon ID, protocol, type and country only. To build it:

```bash
git clone https://github.com/moricef/dec406 ~/dec406
make -C ~/dec406 build/dec406_hex
```

The installer finds `~/dec406/build/dec406_hex` on its own.

## Install as a station (systemd)

```bash
git clone <this repository> meosar-rx && cd meosar-rx
./install.sh            # asks for source, gain, bias-tee, station name and position; --help for the rest
./install.sh --update   # git pull, rebuild, test, reinstall, restart
```

This installs the receiver (`meosar-rx.service`, daily archive in `/var/lib/meosar-rx/archive`),
the MEOSAR Observatory dashboard on :8896 (`meosar-observatory`), a rollup that folds the archive
into `/var/lib/meosar-rx/observatory.db` every minute (`meosar-rollup.timer`), an hourly job that
estimates beacon positions (`meosar-geoloc.timer`, idle priority) and a daily timer that compresses
and expires the raw archive. Settings are in `/etc/meosar-rx.conf`. Run it as your normal user; it
uses sudo when it needs to. `./install.sh --uninstall` removes the services and programs and keeps
the settings, archive and database.

The station position (`STATION_POS`, to ~1 km) is what satellite attribution, beacon positions and
the globe need; without it the dashboard still shows beacons and bursts.

One SDR serves one program. If another service on the station (iridium-sniffer, say) uses the same SDR,
the installer can make them take turns: `sudo systemctl start meosar-rx` stops the other, and
starting the other stops meosar-rx.

**The dashboard is meant for your own network.** It listens on every interface, has no login, and
shows every beacon heard with its identification and position. Don't expose port 8896 to the
internet; put it behind a VPN or an authenticating proxy if you want it from outside.

## MEOSAR Observatory

`observatory/rollup.py` merges the relayed copies of each burst into one row, decodes every
distinct frame once (with dec406's hex decoder when it is installed: protocol, identification,
position), and counts candidates per hour. `observatory/server.py` serves the dashboard: beacons
of both generations with their fields and positions, bursts over time, satellites per burst, the
receiver's status, and each beacon's burst history. A beacon is shown once it is confirmed: a
strict decode, or bursts from two separate transmissions.

Satellite attribution (`observatory/satattr.py`): for beacons whose frames carry a position, the
rollup works out which satellite relayed each copy of a burst. Under the right assignment of copies
to satellites (orbits from CelesTrak's TLEs, propagated with SGP4), every copy implies the same
beacon frequency and transmit time once each satellite's Doppler and path delay are removed. A
burst heard through one satellite only is attributed by frequency, when the same beacon's frequency
is known from a multi-copy burst within the hour. The dashboard lists the satellites (where each is
now, copies relayed, beacons, the elevations they relay at) and the satellites of every burst, and
shows the receiver's frequency reference error, which is the common offset of the implied beacon
frequency (second-generation beacons hold 406.05 MHz tightly). This needs `python3-sgp4` and the
station's position to ~1 km in `STATION_POS`; a wrong position shows as bursts that fit no
assignment (`proto/stationfit.py` finds the position from the bursts themselves). TLEs are fetched
once a day into `/var/lib/meosar-rx/tle/`, the only internet access the observatory makes; offline,
attribution waits for TLEs at most 30 days from each burst.

Estimated positions (`observatory/geoloc.py`): the same relay geometry run the other way. A beacon's
copies only agree on one frequency and transmit time at the beacon's true position, so a search over
trial positions, scored over its multi-copy bursts from the last 72 hours as the satellites move,
finds where it is, whether or not its frames say. Checked against beacons that report their
precise position: second generation within ~2 km, first generation within 3-13 km (the stated 95%
radius is never smaller than 3 and 15 km). Many first-generation frames carry only a coarse
position (0.5 degrees, about 50 km), where the estimate is the better of the two. A beacon needs 10 multi-copy bursts over an hour at least; moving
beacons fit nowhere and get no estimate. The dashboard and the globe show an estimate (marked ≈,
circled at 95% on the globe) for beacons whose frames carry no position, and for the others how far
it is from their own position. One beacon takes 5-80 s of CPU; the hourly job runs at idle priority
and re-estimates beacons still heard every 6 hours.

`/globe` lists the beacons on the left and maps their positions on a 3D globe on the right: colour
by beacon type, first generation round and second generation diamond, tracks joining a beacon's
positions (broken where the implied speed is impossible, as with test beacons that alternate a
placeholder position), and the station itself when `STATION_POS` is set in `/etc/meosar-rx.conf`.
With satellite attribution it also marks the SAR satellites in view now, and links a selected
beacon to the satellites that relayed its latest burst and on to the station.
The globe is rtl-sdr-workbench's renderer with its vendored assets in `observatory/static/`
(three.js, MIT; NASA Black Marble, public domain; world-atlas coastlines, ISC), so it needs no
internet.

Results from recordings can be added at any time:

```bash
python3 /opt/meosar-rx/observatory/rollup.py --archive /var/lib/meosar-rx/archive \
    --db /var/lib/meosar-rx/observatory.db results/*.jsonl
```
