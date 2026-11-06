"""Starter tests (run: pytest -q). Add your own for WP4 and WP6."""
import sys, os, math, tempfile; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from rtops import config
from rtops.pointing import PointingModel
from rtops.sky import geo_azel, refraction_deg
from rtops.tcs import TCS
from rtops.reduce import reduce_solar
from rtops.archive import Archive
from rtops.ingest import ingest_solar
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "experiments"))
import servo_lib

def test_pointing_model_inverse_is_exact():
    m = PointingModel([0.35, -0.2, 0.08, 0.1, -0.15]); az, el = np.array([123.0, 250.0]), np.array([40.0, 60.0]); ea, ee = m.sky_to_enc(az, el); da, de = m.delta(ea, ee)
    assert np.allclose(ea + da, az, atol=1e-6) and np.allclose(ee + de, el, atol=1e-6)

def test_pointing_fit_recovers_terms():
    truth = PointingModel([0.35, -0.2, 0.08, 0.1, -0.15]); r = np.random.default_rng(1); az, el = r.uniform(60, 300, 40), r.uniform(20, 70, 40); da, de = truth.delta(az, el)
    fit, err, rms = PointingModel.fit(az, el, da, az, el, de); assert np.allclose(fit.vector(), truth.vector(), atol=1e-6)

def test_geo_satellite_elevation_is_sensible_for_roorkee():
    az, el = geo_azel(29.8649, 77.8966, 268, 77.9)             # satellite on the observer's meridian: due south
    assert abs(az - 180) < 0.5 and 55 < el < 60

def test_refraction_size():
    assert 0.012 < float(refraction_deg(45.0)) < 0.02 and float(refraction_deg(10.0)) > 0.08

def test_slew_time_matches_trapezoid_theory_roughly():
    cfg = config.load("config/mini.yaml"); t = servo_lib.slew_time(cfg, 30.0, 0.0); th = servo_lib.slew_theory(30.0, cfg["mount"]["vmax_az"], cfg["mount"]["accel"])
    assert th < t < th * 2.0

def test_mini_session_aborts_in_wind_and_resumes():
    cfg = config.load("config/mini.yaml"); out = tempfile.mkdtemp(); tcs = TCS(cfg, out); res = tcs.run_session(); ids = [e["obid"] for e in tcs.arc.entries()]
    assert "HITRK" in ids and "HITRK_r1" in ids and len(res["stow_events"]) == 1
    assert res["pointing_runs"][0]["map_rms_true_after_deg"] < 0.2 * res["pointing_runs"][0]["map_rms_true_before_deg"]

def test_real_data_adapter_round_trip_is_identical():
    cfg = config.load("config/mini_solar.yaml"); out = tempfile.mkdtemp(); tcs = TCS(cfg, out); tcs.run_session(); exp = os.path.join(out, "exp"); os.makedirs(exp)
    ref = reduce_solar(tcs.arc, cfg, os.path.join(out, "a"), tmp=exp); arc2 = Archive(os.path.join(out, "real"))
    for k in (1, 2): ingest_solar(arc2, cfg, f"{exp}/solar_transit_{k}.npz", f"{exp}/solar_transit_{k}.json", f"REAL{k}")
    new = reduce_solar(arc2, cfg, os.path.join(out, "b")); assert abs(new["fwhm_deg"] - ref["fwhm_deg"]) < 1e-9 and new["origin"] == "REAL"
