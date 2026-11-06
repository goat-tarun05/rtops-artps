"""python -m rtops.run_session [--config config/session.yaml] [--out session] [--no-ab]"""
import sys, os, json, time, copy, shutil, argparse, numpy as np
from . import config
from .tcs import TCS
from .reduce import reduce_all, reduce_solar, reduce_hi
from .ingest import ingest_solar, ingest_hi
from .archive import Archive
from .pointing import PointingModel
from . import report

def metrics(tcs, R): return dict(solar=R.get("solar"), sat=R.get("sat"), tracking_rms=[t["rms_true_deg"] for t in tcs.results["tracking"]] or [float("nan")])

def roundtrip(arc, cfg, out):
    """export archived blocks to plain raw files, ingest them as REAL into a second archive, reduce again, compare"""
    rows = []; arc2 = Archive(os.path.join(out, "archive_real")); exp = os.path.join(out, "export"); os.makedirs(exp, exist_ok=True)
    try: ref = reduce_solar(arc, cfg, os.path.join(out, "rt_a"), tmp=exp)
    except Exception: ref = None
    if ref:
        for k in (1, 2): ingest_solar(arc2, cfg, f"{exp}/solar_transit_{k}.npz", f"{exp}/solar_transit_{k}.json", f"REALSOL{k}")
        new = reduce_solar(arc2, cfg, os.path.join(out, "rt_b"))
        for q, lab in (("fwhm_deg", "beam FWHM (°)"), ("pointing_dy_deg", "pointing offset in elevation (°)"), ("sun_peak_k", "Sun height (K)")): rows.append(["Sun cross-scan", lab, f"{ref[q]:.4f}", f"{new[q]:.4f}", f"{new[q] - ref[q]:+.1e}"])
    exp2 = os.path.join(out, "export_hi"); os.makedirs(exp2, exist_ok=True)
    try: refh = reduce_hi(arc, cfg, "HIDFT", os.path.join(out, "rt_c"), tmp=exp2)
    except Exception: refh = None
    if refh:
        ingest_hi(arc2, cfg, f"{exp2}/hi_m31.npz", f"{exp2}/hi_m31.json", "REALHI"); newh = reduce_hi(arc2, cfg, "REALHI", os.path.join(out, "rt_d"))
        for q, lab in (("vsys_helio_kms", "v_sys (km/s)"), ("line_snr", "line S/N")): rows.append(["M31 H I (drift)", lab, f"{refh[q]:.3f}", f"{newh[q]:.3f}", f"{newh[q] - refh[q]:+.1e}"])
    return rows

def main(argv=None):
    ap = argparse.ArgumentParser(); ap.add_argument("--config", default=None); ap.add_argument("--out", default="session"); ap.add_argument("--no-ab", action="store_true"); a = ap.parse_args(argv)
    cfg = config.load(a.config); os.makedirs(a.out, exist_ok=True); t0 = time.time(); tcs = TCS(cfg, a.out); res = tcs.run_session(); res["t0"] = tcs.t0; res["events"] = tcs.mon.events; wall = time.time() - t0
    print(f"session simulated in {wall:.1f} s; {len(tcs.arc.entries())} archive entries"); R = reduce_all(tcs.arc, cfg, os.path.join(a.out, "reduced")); ab = None
    if not a.no_ab:
        cfg2 = copy.deepcopy(cfg)
        for o in cfg2["observations"]:
            if o["id"] == "PR1": o["params"]["update_model"] = False
        tcs2 = TCS(cfg2, os.path.join(a.out, "ab_nocal")); tcs2.run_session(); R2 = reduce_all(tcs2.arc, cfg2, os.path.join(a.out, "ab_nocal", "reduced")); ab = metrics(tcs2, R2)
    ing = roundtrip(tcs.arc, cfg, a.out); summary = dict(wall_s=wall, pointing=[{k: v for k, v in p.items() if k != "meas"} for p in res["pointing_runs"]], tracking=res["tracking"], stow=res["stow_events"], ob_log=res["ob_log"],
        reduced={k: ({kk: {q: w for q, w in vv.items() if q not in ("t", "T", "ideal", "range_km", "el")} for kk, vv in v.items()} if k == "hi" else {q: w for q, w in (v or {}).items() if q not in ("t", "T", "ideal", "range_km", "el")}) for k, v in R.items()})
    json.dump(summary, open(os.path.join(a.out, "results.json"), "w"), indent=1, default=float)
    p = report.build(a.out, cfg, res, None, tcs.mon.hk, tcs.arc, R, ab, ing, wall, tcs.truth); print("report ->", p); return res, R

if __name__ == "__main__": main()
