import sys
sys.path.insert(0, ".")
from rtops.experiment import run_variant

experiments = {
    "baseline": {},
    "low_gain": {"mount.kp": 0.5},
    "high_gain": {"mount.kp": 6.0},
    "heavy_wind": {"weather.base_ms": 9.0, "mount.wind_sigma_per_ms": 0.015}
}

for name, overrides in experiments.items():
    print(f"Running {name}...")
    run_variant(overrides, out=f"results/WP2/{name}", config_path="config/mini.yaml", reduce=False)

print("All WP2 Servo experiments generated successfully!")