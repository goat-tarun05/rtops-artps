"""Starlink: orbit propagation from published elements/TLEs, pass prediction through the beam, burst detection,
identification and an RFI mask for the other datasets (orbital tracking)."""
import math, numpy as np, astropy.units as u
from astropy.time import Time
from sgp4.api import Satrec, WGS72
from skyfield.api import load, EarthSatellite, wgs84
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from .common import load_meta, ang_sep

def satrec_from_elements(e):
    s = Satrec()
    s.sgp4init(WGS72, 'i', e["satnum"], e["epoch_jd"] - 2433281.5, 0.0, 0.0, 0.0, e["ecc"], math.radians(e["argp"]),
               math.radians(e["inc"]), math.radians(e["ma"]), e["n_rad_min"], math.radians(e["raan"]))
    return s

def load_catalogue(path):
    """accepts the JSON element list, or a real 3-line TLE file from CelesTrak"""
    if path.endswith(".json"):
        return [(e["name"], satrec_from_elements(e)) for e in load_meta(path)]
    lines = [l.rstrip() for l in open(path) if l.strip()]
    return [(lines[i].strip(), Satrec.twoline2rv(lines[i + 1], lines[i + 2])) for i in range(0, len(lines) - 2, 3)]

def run(cfg, beam_fwhm_deg, raw="data/raw", out="outputs"):
    d = np.load(f"{raw}/starlink_scan.npz"); m = load_meta(f"{raw}/starlink_scan.json")
    ts = load.timescale(builtin=True); obs = wgs84.latlon(cfg["site"]["lat_deg"], cfg["site"]["lon_deg"], cfg["site"]["height_m"])
    unix = d["unix"]; T = d["counts"] * m["gain_k_per_count"]; dt = np.median(np.diff(unix))
    # 1. predict every catalogue satellite's track through the beam
    t_c = ts.from_astropy(Time(unix[::100], format="unix")); t_f = Time(unix[::5], format="unix"); tf = ts.from_astropy(t_f)
    preds, tracks = [], []
    for name, sr in load_catalogue(f"{raw}/starlink_elements.json"):
        es = EarthSatellite.from_satrec(sr, ts)
        alt, az, _ = (es - obs).at(t_c).altaz()
        if np.min(ang_sep(m["az_deg"], m["el_deg"], az.degrees, alt.degrees)) > 12: continue
        alt, az, dist = (es - obs).at(tf).altaz()
        sep = ang_sep(m["az_deg"], m["el_deg"], az.degrees, alt.degrees); i = int(np.argmin(sep))
        if sep[i] < 2.0 * beam_fwhm_deg:
            rate = np.gradient(sep, unix[::5]); preds.append(dict(name=name, unix=float(unix[::5][i]), min_sep_deg=float(sep[i]),
                         range_km=float(dist.km[i]), exp_rel_amp=float(np.exp(-4 * np.log(2) * sep[i] ** 2 / beam_fwhm_deg ** 2))))
            tracks.append((name, az.degrees, alt.degrees))
    # 2. detect bursts: robust baseline + 6 sigma
    from scipy.ndimage import median_filter
    base = median_filter(T, size=int(60 / dt)); r = T - base
    sig = 1.4826 * np.median(np.abs(r - np.median(r))); hot = r > 6 * sig
    bursts, i = [], 0
    while i < len(hot):
        if hot[i]:
            j = i
            while j < len(hot) and hot[j]: j += 1
            if j - i >= 3:
                k = i + int(np.argmax(r[i:j])); bursts.append(dict(unix=float(unix[k]), peak_k=float(r[k]), dur_s=float((j - i) * dt), i0=i, i1=j))
            i = j
        else: i += 1
    # 3. identify: nearest predicted pass within +-5 s
    used = set()
    for b in bursts:
        c = [(abs(p["unix"] - b["unix"]), n) for n, p in enumerate(preds) if n not in used and abs(p["unix"] - b["unix"]) < 5]
        if c:
            dtm, n = min(c); used.add(n); b.update(id=preds[n]["name"], dt_s=float(b["unix"] - preds[n]["unix"]), min_sep_deg=preds[n]["min_sep_deg"])
        else: b.update(id="UNIDENTIFIED (terrestrial RFI?)", dt_s=None, min_sep_deg=None)
    missed = [p for n, p in enumerate(preds) if n not in used]
    mask = [[Time(unix[b["i0"]] - 3, format="unix").isot, Time(unix[min(b["i1"], len(unix) - 1)] + 3, format="unix").isot] for b in bursts]
    ident = [b for b in bursts if b["dt_s"] is not None]
    res = dict(n_predicted=len(preds), n_detected=len(bursts), n_identified=len(ident), n_unidentified=len(bursts) - len(ident),
               n_missed=len(missed), id_rate=float(len(ident) / max(1, len([b for b in bursts]))),
               timing_rms_s=float(np.sqrt(np.mean([b["dt_s"] ** 2 for b in ident]))) if ident else None,
               masked_fraction=float(sum(b["i1"] - b["i0"] + 60 for b in bursts) / len(unix)),
               detections=[{k: v for k, v in b.items() if k not in ("i0", "i1")} for b in bursts],
               missed=missed, rfi_mask_utc=mask)
    strong = [n for n, p in enumerate(preds) if p["exp_rel_amp"] > 0.3]
    res["strong_predicted"] = len(strong); res["strong_identified"] = len([n for n in strong if n in used])
    res["checks"] = {">=80% of strong predicted passes detected and identified": len(strong) == 0 or res["strong_identified"] / len(strong) >= 0.8,
                     "timing rms < 2 s": res["timing_rms_s"] is not None and res["timing_rms_s"] < 2}
    mins = (unix - unix[0]) / 60
    fig, ax = plt.subplots(2, 1, figsize=(8, 6.5))
    ax[0].plot(mins, r, lw=0.5, label="power above baseline (K)")
    for b in bursts:
        lab = b["id"].replace("STARLINK-SIM-", "SL-") if b["dt_s"] is not None else "RFI?"
        ax[0].annotate(lab, ((b["unix"] - unix[0]) / 60, b["peak_k"]), fontsize=7, rotation=60, color="g" if b["dt_s"] is not None else "r")
    for p in preds: ax[0].axvline((p["unix"] - unix[0]) / 60, color="orange", lw=0.6, alpha=0.6)
    ax[0].set_xlabel("minutes"); ax[0].set_ylabel("K"); ax[0].set_title(f"Starlink: {len(ident)}/{len(bursts)} bursts identified from orbits, orange = predicted passes")
    th = np.linspace(0, 2 * np.pi, 100)
    for name, az, el in tracks:
        xs = (az - m["az_deg"]) * np.cos(np.radians(m["el_deg"])); ys = el - m["el_deg"]; ok = (np.abs(xs) < 10) & (np.abs(ys) < 10)
        ax[1].plot(xs[ok], ys[ok], lw=0.8)
    ax[1].plot(beam_fwhm_deg / 2 * np.cos(th), beam_fwhm_deg / 2 * np.sin(th), "k-", lw=2, label="beam (FWHM from solar fit)")
    ax[1].set_aspect("equal"); ax[1].set_xlim(-10, 10); ax[1].set_ylim(-10, 10); ax[1].legend(fontsize=8)
    ax[1].set_xlabel("sky offset in azimuth (deg)"); ax[1].set_ylabel("elevation offset (deg)"); ax[1].set_title("predicted satellite tracks across the beam")
    fig.tight_layout(); fig.savefig(f"{out}/starlink.png", dpi=130); plt.close(fig)
    return res
