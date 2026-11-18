# meosar-rx

A receiver for Cospas-Sarsat 406 MHz distress-beacon bursts relayed on the MEOSAR L-band
downlink (Galileo SAR 1544.05–1544.15 MHz, GLONASS-K 1544.85–1544.95 MHz), built for weak
signals on a small fixed antenna.

Status: Phase 0 (research). `proto/` holds the Python prototype used to design and measure the
receiver; `src/` is the C receiver, which so far implements the first-generation (FGB) chain.

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

## Build and run

```bash
sudo apt install build-essential cmake pkg-config libfftw3-dev libairspy-dev libhackrf-dev
cmake -S . -B build && cmake --build build -j
./build/test_fgb
./build/meosar-rx -i recording.sigmf-meta > records.jsonl
```

One JSON record per burst candidate. With `--archive DIR` records are also appended to one file
per UTC day, and `--stats 600` adds a status record every 10 minutes.

## Install as a station (systemd)

```bash
./install.sh            # asks for source, gain, bias-tee, station name; --help for the rest
./install.sh --update   # pull, rebuild, test, reinstall, restart
```

This installs the receiver (`meosar-rx.service`, daily archive in `/var/lib/meosar-rx/archive`),
the MEOSAR Observatory dashboard on :8896 (`meosar-observatory`), a rollup that folds the archive
into `/var/lib/meosar-rx/observatory.db` every minute (`meosar-rollup.timer`) and a daily timer
that compresses and expires the raw archive. Settings are in `/etc/meosar-rx.conf`.

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

One SDR serves one program. If another service here (iridium-sniffer, say) uses the same SDR,
the installer can make them take turns: `sudo systemctl start meosar-rx` stops the other, and
starting the other stops meosar-rx.
