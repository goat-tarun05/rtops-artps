import sys, os, json
sys.path.insert(0, ".")
import numpy as np
from rtops.experiment import run_variant

os.makedirs("results/WP6", exist_ok=True)
print("--- Running Monte Carlo (10 seeds) ---")
rows = []
for seed in range(21, 31):
    print(f"Running seed {seed}...")
    r = run_variant({"session.seed": seed}, out=f"runs/mc/{seed}")
    R, pr = r["reduced"], r["results"]["pointing_runs"]
    
    hi_trk = R["hi"].get("tracked ON/OFF", {})
    hi_dft = R["hi"].get("fixed-dish drift", {})
    
    rows.append(dict(
        seed=seed, 
        map_after=pr[0]["map_rms_true_after_deg"],
        pr2=pr[1]["rms_resid_before_deg"],
        sat_eff=R.get("sat", {}).get("mean_eff", 0.0), 
        fwhm=R.get("solar", {}).get("fwhm_deg", 0.0), 
        sun=R.get("solar", {}).get("sun_peak_k", 0.0),
        v_trk=hi_trk.get("vsys_helio_kms", 0.0),
        e_trk=hi_trk.get("vsys_err_kms", 0.0),
        snr_trk=hi_trk.get("line_snr", 0.0),
        v_dft=hi_dft.get("vsys_helio_kms", 0.0),
        e_dft=hi_dft.get("vsys_err_kms", 0.0),
        snr_dft=hi_dft.get("line_snr", 0.0)
    ))

json.dump(rows, open("results/WP6/montecarlo.json", "w"), indent=1)
print("Monte Carlo complete. Saved to results/WP6/montecarlo.json")

print("\n--- Running Stress Test ---")
base = [0.35, -0.20, 0.08, 0.10, -0.15]
harsh_config = {
    "telescope.dish_m": 3.0, 
    "pointing_truth": [3 * b for b in base], 
    "mount.wind_sigma_per_ms": 0.0075,
    "satellites.catalogue_err_km": 5.0, 
    "weather.gusts.0.speed": 20.0, 
    "session.seed": 7
}
r_stress = run_variant(harsh_config, out="runs/stress")
print(f"Harsh PR1 Sky Error: {r_stress['results']['pointing_runs'][0]['map_rms_true_after_deg']:.3f}°")
print(f"Harsh Sat Efficiency: {r_stress['reduced'].get('sat', {}).get('mean_eff', 0.0):.3f}")

print("\nUndoing dish size constraint (reverting to 1.5m)...")
harsh_config["telescope.dish_m"] = 1.5
r_stress_fix = run_variant(harsh_config, out="runs/stress_fix")
print(f"Fixed PR1 Sky Error: {r_stress_fix['results']['pointing_runs'][0]['map_rms_true_after_deg']:.3f}°")
print(f"Fixed Sat Efficiency: {r_stress_fix['reduced'].get('sat', {}).get('mean_eff', 0.0):.3f}")