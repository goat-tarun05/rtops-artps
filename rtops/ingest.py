"""Adapters that bring REAL recordings into the same archive, so the reduction and reporting run unchanged.
  solar : unix, counts (total power), hot, cold arrays  +  nominal az/el, elevation offset     (rtpipe 'solar_transit' format)
  hi    : freq_hz, unix, spec (time x channel), hot, [is_on]  +  mode, target, az/el or OFF position   (rtpipe 'hi_m31' format)"""
import json, numpy as np
from .archive import Archive

def ingest_solar(arc, cfg, npz_path, meta_path, obid, origin="REAL"):
    d = np.load(npz_path); m = json.load(open(meta_path)); t0, t1 = float(d["unix"][0]), float(d["unix"][-1]); cal = f"CAL_KU_{obid}"
    arc.write(cal, "cal_ku", dict(TSTART=t0 - 300, TSTOP=t0 - 60), dict(cold=dict(t=np.arange(len(d["cold"]), dtype=float) + t0 - 300, counts=d["cold"]), hot=dict(t=np.arange(len(d["hot"]), dtype=float) + t0 - 200, counts=d["hot"])), origin=origin, target="calibration")
    return arc.write(obid, "solar_drift", dict(TSTART=t0, TSTOP=t1, AZ_CMD=m["az_deg"], EL_CMD=m["el_deg"], EL_OFFS=m.get("el_offset_set_deg", 0.0), FREQ_GHZ=m.get("freq_ghz", cfg["ku"]["freq_ghz"])), dict(data=dict(t=d["unix"], counts=d["counts"])), origin=origin, target="Sun", cal_ref=cal)

def ingest_hi(arc, cfg, npz_path, meta_path, obid, origin="REAL"):
    d = np.load(npz_path); m = json.load(open(meta_path)); t0, t1 = float(d["unix"][0]), float(d["unix"][-1]); cal = f"CAL_HI_{obid}"; mode = m["mode"]; n = len(d["unix"])
    arc.write(cal, "cal_hi", dict(TSTART=t0 - 300, TSTOP=t0 - 60), dict(hot=dict(t=np.arange(len(d["hot"]), dtype=float) + t0 - 300, spec=d["hot"]), freq=dict(freq_hz=d["freq_hz"])), origin=origin, target="calibration")
    hdr = dict(TSTART=t0, TSTOP=t1, BASEID=obid, DUMP_S=float(np.median(np.diff(d["unix"]))), MODE=mode)
    if mode == "switch": hdr.update(OFF_RA=m["off"]["ra_deg"], OFF_DEC=m["off"]["dec_deg"])
    else: hdr.update(AZ_CMD=m["az_deg"], EL_CMD=m["el_deg"])
    on = d["is_on"].astype(bool) if "is_on" in d.files else np.ones(n, bool)
    return arc.write(obid, "hi_track_nod" if mode == "switch" else "hi_stare", hdr, dict(data=dict(t=d["unix"], spec=d["spec"], resp=np.ones(n), theta=np.zeros(n), is_on=on), freq=dict(freq_hz=d["freq_hz"])), origin=origin, target="M31", cal_ref=cal)
