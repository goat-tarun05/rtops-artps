# Code Tour

| Function | Question Answered |
| :--- | :--- |
| `TCS.slew` | It decides the telescope has settled by checking if the pointing error stays below a threshold for a continuous `hold_s` duration. |
| `TCS._drive` | Target sky position → inverse pointing model → encoder command → mount servo → add hidden error and wind → real boresight. |
| `Mount.run` | It enforces rate (`vmax_az`, `vmax_el`) and acceleration (`accel`) limits, and adds wind by injecting normally distributed noise scaled by `wind_sigma_per_ms`. |
| `PointingModel.sky_to_enc` | It is iterated four times because the pointing corrections (like NPAE) depend on the elevation, creating a non-linear set of equations that must converge. |
| `Dispatcher.next` | Tests fixed blocks within the lead time first, then flexible blocks by priority (if they fit the duration and constraints), then defaults to idle. |
| `TCS.run_pointing_run` | Measures the azimuth and elevation offsets by driving a cross-scan over a source to fit a Gaussian peak, which are then used to fit the 5 pointing model terms. |
| `Receiver.ku_temperature` | It feeds the beam response using the *real* commanded position (which includes hidden errors and wind), not the ideal target position. |
| `Archive.write` | Observation metadata (parameters, target, times) goes in the FITS header, and the time-series measurements (counts, t, az/el) go into the binary tables. |
| `reduce_solar` | It reads the `solar_drift` science block and its associated `cal_ku` reference block (linked via `cal_ref`). |