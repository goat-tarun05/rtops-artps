Markdown
# WP2: Servo and Tracking Experiments

* **Baseline (`kp=2.0`):** The motor perfectly follows the commanded trajectory. The servo tracking error is 0.0015° rms, which precisely matches the encoder noise floor (0.001°) combined with standard atmospheric wind buffering.
* **Low Gain (`kp=0.5`):** The telescope acts sluggish and heavily damped. It takes nearly twice as long to settle into position after a fast slew, causing the `TCS.slew` function to delay observation starts.
* **High Gain (`kp=6.0`):** The proportional response is too aggressive, causing the telescope to overshoot the target and oscillate (ring) before settling. This degrades the tracking efficiency.
* **Heavy Wind:** By increasing base wind to 9 m/s and wind noise to 0.015, the physical dish gets buffered constantly. The servo following error spikes because the motor torque cannot instantly counter the random aerodynamic forces.