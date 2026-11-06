# WP1 Report Notes (Sections 3-5)

* **True sky pointing error (0.384° → 0.016°):** Measures the absolute pointing error before and after the calibration run; it drops significantly because the newly fitted pointing model corrects the physical hardware offsets.
* **Verification run PR2 (0.016°):** Measures the residual error on independent sources without refitting the model; it confirms the pointing model remains accurate over time.
* **Servo following error (0.0015° rms):** Measures the difference between the encoder reading and the commanded trajectory; it sits at the absolute noise floor of the encoder because feed-forward is enabled.
* **Satellite true tracking error (0.19°) & efficiency (0.94):** Measures how far the boresight lags the fast-moving satellite; the dish follows accurately enough to capture 94% of the maximum possible signal.
* **Satellite time above half power (244 s vs 1.9 s):** Compares observation lengths; active tracking keeps the satellite in the beam roughly 130 times longer than letting it drift through a stationary dish.
* **Sun cross-scan (FWHM 1.328°, 872 K):** Measures the effective beam width and the measured antenna temperature of the Sun; validates the dish diameter and the receiver's calibration.
* **M31 v_sys line S/N (Tracked 5.6, Drift 5.8):** Measures the signal-to-noise ratio of the galactic hydrogen line; shows that both observation modes successfully detect the galaxy with similar confidence.
* **Weather stow (1 event, 21 min):** Measures automated safety interruptions; the dish parked itself for 21 minutes because the wind speed exceeded the 10 m/s safety limit.