"""Science reduction from the ARCHIVE only (works identically for simulated and ingested real recordings): converts archive blocks to the
rtpipe input formats and runs the rtpipe science modules; adds a satellite light-curve analysis."""
import os, json, tempfile, numpy as np
from rtpipe import solar as rt_solar, hi as rt_hi
from .archive import Archive
from rtpipe.common import ang_sep

def _cal(arc, e):
    c = arc.entries(obid=e["cal_ref"]); return arc.read(c[0]["path"]) if c else (None, None)

def reduce_solar(arc, cfg, out, tmp=None):
    ents = [e for e in arc.entries(obtype="solar_drift") if e["status"] in ("DONE", "ABORTED_WIND")]; tmp = tmp or tempfile.mkdtemp(); os.makedirs(out, exist_ok=True)
    if len(ents) < 2: return None
    for k, e in enumerate(ents[:2], 1):
        hdr, tabs = arc.read(e["path"]); _, ct = _cal(arc, e); np.savez(f"{tmp}/solar_transit_{k}.npz", unix=tabs["data"]["t"], counts=tabs["data"]["counts"], hot=ct["hot"]["counts"], cold=ct["cold"]["counts"])
        json.dump(dict(az_deg=float(hdr["AZ_CMD"]), el_deg=float(hdr["EL_CMD"]), el_offset_set_deg=float(hdr["EL_OFFS"]), freq_ghz=float(hdr["FREQ_GHZ"])), open(f"{tmp}/solar_transit_{k}.json", "w"))
    r = rt_solar.run(cfg, raw=tmp, out=out); r["blocks"] = [e["obid"] for e in ents[:2]]; r["origin"] = ents[0]["origin"]; return r

def reduce_hi(arc, cfg, base, out, tmp=None):
    ents = [e for e in arc.entries() if e["obtype"].startswith("hi_") and (e["obid"] == base or e["obid"].startswith(base + "_r"))]; tmp = tmp or tempfile.mkdtemp(); os.makedirs(out, exist_ok=True)
    if not ents: return None
    parts = [arc.read(e["path"]) for e in ents]; hdr0 = parts[0][0]; mode = "switch" if ents[0]["obtype"] == "hi_track_nod" else "drift"; _, ct = _cal(arc, ents[0])
    freq = parts[0][1]["freq"]["freq_hz"]; unix = np.concatenate([p[1]["data"]["t"] for p in parts]); spec = np.concatenate([p[1]["data"]["spec"] for p in parts]); on = np.concatenate([p[1]["data"]["is_on"] for p in parts]).astype(bool)
    sv = dict(freq_hz=freq, unix=unix, spec=spec, hot=ct["hot"]["spec"]); meta = dict(mode=mode, target=cfg["targets"]["m31"], rest_mhz=cfg["hi"]["rest_mhz"])
    if mode == "switch": sv["is_on"] = on; meta["off"] = dict(ra_deg=float(hdr0["OFF_RA"]), dec_deg=float(hdr0["OFF_DEC"]))
    else: meta["az_deg"], meta["el_deg"] = float(hdr0["AZ_CMD"]), float(hdr0["EL_CMD"])
    np.savez(f"{tmp}/hi_m31.npz", **sv); json.dump(meta, open(f"{tmp}/hi_m31.json", "w"))
    r = rt_hi.run(cfg, raw=tmp, out=out); r["blocks"] = [e["obid"] for e in ents]; r["mode_used"] = mode; r["n_dumps"] = int(len(unix)); r["origin"] = ents[0]["origin"]
    r["on_source_hours"] = float(on.sum() * float(hdr0["DUMP_S"]) / 3600) if mode == "switch" else float(len(unix) * float(hdr0["DUMP_S"]) / 3600); return r

def reduce_sat(arc, cfg, e):
    hdr, tabs = arc.read(e["path"]); _, ct = _cal(arc, e); ku, rc = cfg["ku"], cfg["receiver"]; Y = ct["hot"]["counts"].mean() / ct["cold"]["counts"].mean()
    t_rx = (ku["t_hot_k"] - Y * ku["t_cold_k"]) / (Y - 1); g = ct["cold"]["counts"].mean() / (t_rx + ku["t_cold_k"]); t, c = tabs["data"]["t"], tabs["data"]["counts"]; eph = tabs["ephem"]
    az, el, rng = np.interp(t, eph["t"], eph["az"]), np.interp(t, eph["t"], eph["el"]), np.interp(t, eph["t"], eph["range_km"])
    sky = 2.7 + rc["t_atm_zenith"] / np.sin(np.radians(np.clip(el, 10, 90))); T = c / g - t_rx - sky; ideal = cfg["satellites"]["t_sat_k"] * (cfg["satellites"]["ref_range_km"] / rng) ** 2
    dt = float(np.median(np.diff(t))); half = T > 0.5 * ideal; eff = T / ideal; i = int(np.argmax(ideal)); rate = np.degrees(np.abs(np.gradient(np.unwrap(np.radians(az)), t))) * np.cos(np.radians(el)); rate_el = np.abs(np.gradient(el, t)); ang_rate = np.hypot(rate, rate_el)
    fw = cfg["ku"]["fwhm_theory_deg"]; noise = float(np.std(np.diff(T)) / np.sqrt(2))
    return dict(name=hdr["SATNAME"], block=e["obid"], t=t.tolist(), T=T.tolist(), ideal=ideal.tolist(), range_km=rng.tolist(), el=el.tolist(), max_el=float(hdr["MAX_EL"]), duration_s=float(t[-1] - t[0]),
                peak_k=float(T.max()), mean_eff=float(np.mean(eff)), time_above_half_s=float(half.sum() * dt), snr_peak=float(T.max() / max(noise, 1e-6)), noise_k=noise,
                fixed_beam_equiv_s=float(fw / ang_rate[i]), gain_k_per_count=float(1 / g), origin=e["origin"])

def reduce_all(arc, cfg, out):
    """every product is reduced independently; a product that cannot be reduced is reported in res["errors"] instead of stopping the run"""
    res = {"hi": {}, "errors": {}}
    try:
        sd = reduce_solar(arc, cfg, os.path.join(out, "solar"))
        if sd: res["solar"] = sd
    except Exception as e: res["errors"]["solar"] = f"{type(e).__name__}: {e}"
    for base, label in (("HITRK", "tracked ON/OFF"), ("HIDFT", "fixed-dish drift")):
        try:
            r = reduce_hi(arc, cfg, base, os.path.join(out, base.lower()))
            if r: res["hi"][label] = r
        except Exception as e: res["errors"][f"hi {label}"] = f"{type(e).__name__}: {e}"
    try:
        sats = [e for e in arc.entries(obtype="sat_track") if e["status"] == "DONE"]
        if sats: res["sat"] = reduce_sat(arc, cfg, sats[0])
    except Exception as e: res["errors"]["sat"] = f"{type(e).__name__}: {e}"
    return res
