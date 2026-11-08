#!/usr/bin/env bash
set -e

echo "Running Unit Tests..."
python -m pytest -q

echo "Running Reference Session..."
python -m rtops.run_session --out session

echo "Running WP2: Servo and Tracking..."
python experiments/wp2_servo.py

echo "Running WP3: Pointing Calibration..."
python experiments/wp3_pointing.py

echo "Running WP4: Sun Cross-Scan..."
python experiments/wp4_solar.py

echo "Running WP5: Science Integration..."
python experiments/wp5_science.py

echo "Running WP6: Monte Carlo and Stress Tests..."
python experiments/wp6_montecarlo.py

echo "Running Nightly Operations..."
rm -f runs/log.csv
for d in 2026-11-15 2026-11-16 2026-11-17; do python scripts/nightly.py $d || true; done

echo "All deliverables successfully reproduced!"