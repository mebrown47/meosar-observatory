# proto/: the research prototype

The Python prototype used to design and measure the receiver before and alongside the C code in
`src/`. **A station doesn't need any of it**: `install.sh` builds and runs the C receiver and the
observatory in `observatory/`. It is kept because it records how each part of the receiver was
chosen and measured, and the synthetic tests are still the way to compare changes near threshold.

Needs Python 3 with NumPy; `satid.py` and `stationfit.py` also need `python3-sgp4`. Recordings are
SigMF, complex int16 (`ci16`) at 2.5 MS/s centred on 1544.5 MHz unless a script says otherwise.

## Reference implementation (Python versions of what `src/` does)

| Script | What it is |
|---|---|
| `dsp.py` | Signal helpers: reading ci16, channelising, FGB modulation for synthetic tests |
| `fgb.py` | First-generation frame: layout, BCH codes, Chase soft decoding |
| `demod.py` | Coherent FGB demodulator for weak bursts, on a whole burst in memory |
| `run.py` | The FGB pipeline over a recording: carrier detector, demodulator, soft decoder |
| `group.py` | Groups accepted decodes into bursts and beacons |
| `combine.py` | Second pass over a run's records: combines repeats of one beacon across bursts |
| `sgb.py` | Second-generation beacon: spreading codes, frame, BCH(250,202), soft decoding |
| `sgb_rx.py` | SGB receiver for weak bursts: preamble acquisition, coherent despreading |
| `sgb_scan.py` | Searches a recording for SGB bursts in one band |
| `sgb_pair.py` | Demodulates given SGB relays of one burst, combines their soft bits, decodes |

## Measurement

| Script | What it measures |
|---|---|
| `synth.py`, `synth_sgb.py` | Inject synthetic FGB / SGB bursts at a known C/N0 into real recorded noise (the quiet band between the SAR bands) |
| `score.py` | Scores decodes against the synthetic truth |
| `scoreboard.sh` | Synthetic scoreboard across C/N0: this receiver against dec406 |
| `ablate.py` | How much each receiver stage contributes, from one run's records |
| `oracle_fgb.py`, `oracle_snr.py`, `loss_split.py` | Where the FGB chain loses near threshold, with oracle timing and frequency |
| `dec406_cands.py`, `win2457.py` | Run dec406's FGB decoder on our candidates (it is validated at 2.4576 MS/s) |
| `try_real.py` | Two bursts of the first Phase 0 recording, step by step |
| `process_chunks.sh` | Runs the C receiver on each finished chunk of a long recording, once |

## Satellites and positions

| Script | What it does |
|---|---|
| `satid.py` | The first version of satellite attribution (now `observatory/satattr.py`) |
| `stationfit.py` | Fits the receiving station's position from relayed SGB bursts; useful to check `STATION_POS` |

`observatory/geoloc.py` grew out of this directory too and is now part of the observatory.
