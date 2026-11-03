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

## Build and run

```bash
sudo apt install build-essential cmake pkg-config libfftw3-dev libairspy-dev libhackrf-dev
cmake -S . -B build && cmake --build build -j
./build/test_fgb
./build/meosar-rx -i recording.sigmf-meta > records.jsonl
```

One JSON record per burst candidate. Records include beacon IDs and positions, which can belong
to real distress beacons: keep them on the station unless the beacon is a test, system or
orbitography beacon.
