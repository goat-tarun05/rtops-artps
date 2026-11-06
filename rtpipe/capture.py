"""SDR capture that writes the raw formats run_all.py reads.  No hardware? add --fake to test the whole chain.
  python -m rtpipe.capture solar    --day 1 --az 185.4 --el 42.5 --el-offset 1.0 --seconds 3600 --fake
  python -m rtpipe.capture starlink --az 180 --el 55 --seconds 2400 --fake
  python -m rtpipe.capture hi       --az 0 --el 78.6 --target 10.6847 41.269 --seconds 43200 --fake
Real hardware: pip install pyrtlsdr (RTL-SDR).  Ku: tune to RF minus the LNB LO.  HI needs >= 5 MHz bandwidth."""
import argparse, json, os, time
import numpy as np

class Source:
    def __init__(self, center_hz, rate, gain, fake):
        self.fake, self.rate = fake, rate
        if not fake:
            from rtlsdr import RtlSdr
            self.sdr = RtlSdr(); self.sdr.sample_rate = rate; self.sdr.center_freq = center_hz; self.sdr.gain = gain
    def read(self, n):
        if self.fake:
            return (np.random.randn(n) + 1j * np.random.randn(n)) / np.sqrt(2)
        return self.sdr.read_samples(n)
    def close(self):
        if not self.fake: self.sdr.close()

def total_power(src, seconds, cadence_s, nread=262144):
    unix, cnt, t_end = [], [], time.time() + seconds
    while time.time() < t_end:
        t0, acc = time.time(), []
        while time.time() - t0 < cadence_s: acc.append(np.mean(np.abs(src.read(nread)) ** 2))
        unix.append(t0 + cadence_s / 2); cnt.append(np.mean(acc))
    return np.array(unix), np.array(cnt)

def spectra(src, seconds, dump_s, nch=8192, frames=64):
    unix, spec, t_end = [], [], time.time() + seconds
    while time.time() < t_end:
        t0, acc, n = time.time(), np.zeros(nch), 0
        while time.time() - t0 < dump_s:
            x = src.read(nch * frames).reshape(frames, nch) * np.hanning(nch)
            acc += (np.abs(np.fft.fftshift(np.fft.fft(x, axis=1), axes=1)) ** 2).sum(0); n += frames
        unix.append(t0 + dump_s / 2); spec.append(acc / n)
    return np.array(unix), np.array(spec, dtype=np.float32)

def _wait(msg, noprompt):
    if not noprompt: input(msg + "  [Enter when ready] ")

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("mode", choices=["solar", "starlink", "hi"])
    ap.add_argument("--day", type=int, default=1); ap.add_argument("--az", type=float, required=True); ap.add_argument("--el", type=float, required=True)
    ap.add_argument("--el-offset", type=float, default=0.0); ap.add_argument("--seconds", type=float, default=3600); ap.add_argument("--hot", type=float, default=120)
    ap.add_argument("--cold", type=float, default=120); ap.add_argument("--cadence", type=float, default=1.0); ap.add_argument("--dump", type=float, default=60)
    ap.add_argument("--center-hz", type=float, default=1.5e9); ap.add_argument("--rate", type=float, default=2.4e6); ap.add_argument("--gain", type=float, default=30)
    ap.add_argument("--target", type=float, nargs=2, default=[10.6847, 41.269]); ap.add_argument("--freq-ghz", type=float, default=11.2)
    ap.add_argument("--fake", action="store_true"); ap.add_argument("--no-prompt", action="store_true"); ap.add_argument("--out", default="data/raw")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); src = Source(a.center_hz, a.rate, a.gain, a.fake)
    if a.mode == "solar":
        _wait("Point the feed at the GROUND (hot load)", a.no_prompt); _, hot = total_power(src, a.hot, a.cadence)
        _wait("Point at EMPTY SKY (cold)", a.no_prompt); _, cold = total_power(src, a.cold, a.cadence)
        _wait(f"Set the dish to az {a.az}, el {a.el} (offset {a.el_offset:+.1f})", a.no_prompt); unix, cnt = total_power(src, a.seconds, a.cadence)
        np.savez(f"{a.out}/solar_transit_{a.day}.npz", unix=unix, counts=cnt, hot=hot, cold=cold)
        json.dump(dict(az_deg=a.az, el_deg=a.el, el_offset_set_deg=a.el_offset, freq_ghz=a.freq_ghz), open(f"{a.out}/solar_transit_{a.day}.json", "w"))
    elif a.mode == "starlink":
        unix, cnt = total_power(src, a.seconds, 0.1)
        np.savez(f"{a.out}/starlink_scan.npz", unix=unix, counts=cnt)
        json.dump(dict(az_deg=a.az, el_deg=a.el, freq_ghz=a.freq_ghz, gain_k_per_count=1.0), open(f"{a.out}/starlink_scan.json", "w"))
        print("edit gain_k_per_count in the json (kelvin per count from the solar Y-factor) before running run_all.py")
    else:
        _wait("Point at the GROUND (hot load)", a.no_prompt); _, hot = spectra(src, 10 * a.dump, a.dump)
        _wait(f"Set the antenna to az {a.az}, el {a.el} (fixed)", a.no_prompt); unix, spec = spectra(src, a.seconds, a.dump)
        np.savez(f"{a.out}/hi_m31.npz", freq_hz=a.center_hz + (np.arange(spec.shape[1]) - spec.shape[1] / 2) * a.rate / spec.shape[1], unix=unix, spec=spec, hot=hot[:10])
        json.dump(dict(mode="drift", target=dict(ra_deg=a.target[0], dec_deg=a.target[1]), az_deg=a.az, el_deg=a.el, rest_mhz=1420.40575), open(f"{a.out}/hi_m31.json", "w"))
    src.close(); print("saved to", a.out)

if __name__ == "__main__":
    main()
