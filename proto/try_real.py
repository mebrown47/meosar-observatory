import sys, numpy as np
import dsp, demod, fgb
FS, FC = 2.5e6, 1544.5e6
fn = sys.argv[1]
for t_win, t_cw, f_mhz in [(0.0, 0.85, 1544.0941), (4.6, 5.83, 1544.0700)]:
    x = dsp.read_ci16(fn, t_win, 2.0, FS)
    # coarse: 100 ms FFT around the guess for a +-80 Hz start
    xb = dsp.channelise(x, (f_mhz * 1e6 - FC), FS, 125, 8000)
    fs = FS / 125
    seg = xb[int((t_cw - t_win) * fs): int((t_cw - t_win + 0.15) * fs)]
    S = np.abs(np.fft.fft(seg, 1 << 16)); fr = np.fft.fftfreq(1 << 16, 1 / fs); m = np.abs(fr) < 600
    f1 = fr[m][np.argmax(S[m])]
    xb = dsp.channelise(x, (f_mhz * 1e6 - FC) + f1, FS, 125, 8000)
    d = demod.demodulate(xb, fs, t_cw - t_win)
    print(f"== burst near {t_cw}s {f_mhz}: coarse {f1:+.1f} Hz -> f {d.f_hz:+.2f} Hz fdot {d.fdot:+.2f} Hz/s  data at {t_win + d.t_data:.4f}s  pol {d.polarity} selftest {d.selftest_sync}  sync score {d.sync_score:.2f}  Eb/N0 {d.ebn0_db:.1f} dB  C/N0 {d.cn0_dbhz:.1f}")
    if d.ok:
        hard = (d.llr > 0).astype(int)
        print("   hard:", "".join(map(str, hard)))
        for L in (8, 12, 14):
            r = fgb.decode_soft(d.llr, L=L)
            print(f"   chase L={L}: ok={r.ok} long={r.long} cost {r.bch1_cost:.1f}+{r.bch2_cost:.1f} flips {r.bch1_flips}+{r.bch2_flips} reasons {r.reasons} hex {r.hex}")
        print("   |llr| median", np.median(np.abs(d.llr)).round(2), " preamble+sync errors:", int(np.sum(hard[:24] != (demod.KNOWN_ST if d.selftest_sync else demod.KNOWN))))
