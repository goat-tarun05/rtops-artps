# WP4: Sun Cross-Scan & Data Reduction

* **Pipeline Architecture (`reduce_solar`):** 
  1. **Extraction:** It filters the session archive for the two `solar_drift` observation blocks (`SOL1` and `SOL2`, which scan slightly above and below the Sun).
  2. **Calibration:** It calls the `_cal()` helper function to find the associated reference block. This allows the software to convert raw ADC power counts into physical Kelvin temperatures using the known values of the internal Hot and Cold loads.
  3. **Formatting:** It packages the time series, encoder angles, and calibrated temperatures into standard NumPy `.npz` arrays.
  4. **Fitting:** It passes the clean arrays to the `rtpipe.solar` module, which performs a least-squares Gaussian fit to find the exact peak temperature and Full-Width Half-Maximum (FWHM).
* **Physics Verification:** The theoretical beamwidth of a 1.5m dish at 11.2 GHz is 1.25°. Because the Sun has its own angular diameter of 0.45°, the signal convolves. $\sqrt{1.25^2 + 0.45^2} = 1.328^\circ$, which perfectly matches the digital observatory's baseline output.