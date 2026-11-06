"""Helpers for the servo experiments (WP2). They drive the virtual mount directly, without a whole session."""
import sys, os, math; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from rtops import config
from rtops.mount import Mount

def slew_theory(dist_deg, vmax, accel):
    """ideal bang-bang / trapezoid slew time with rate limit vmax and acceleration limit accel"""
    return dist_deg / vmax + vmax / accel if dist_deg >= vmax ** 2 / accel else 2 * math.sqrt(dist_deg / accel)

def slew_time(cfg, d_az, d_el, wind_ms=4.0, seed=0, tol=None, hold_s=2.0):
    """time for the mount to arrive within `tol` degrees of the target on both axes and stay there for `hold_s` seconds"""
    m = cfg["mount"]; tol = tol or m["settle_tol_deg"]; dt = 1.0 / m["rate_hz"]; mt = Mount(cfg, np.random.default_rng(seed), az0=100.0, el0=30.0)
    n = int(300 / dt); caz = np.full(n, 100.0 + d_az); cel = np.full(n, 30.0 + d_el); z = np.zeros(n); a, e, *_ = mt.run(caz, cel, z, z, np.full(n, wind_ms))
    ok = (np.abs(a - caz) * math.cos(math.radians(40)) < tol) & (np.abs(e - cel) < tol); need = int(hold_s / dt); run = np.convolve(ok, np.ones(need), mode="valid") >= need
    return float((np.argmax(run) + need) * dt) if run.any() else float("nan")

def tracking_error(cfg, rate_deg_s, wind_ms=4.0, duration_s=120.0, seed=0, feedforward=True):
    """follow a target moving at `rate_deg_s` (azimuth sky-plane; feedforward=False removes the rate feed-forward term); returns rms servo following error and rms wind deflection (deg)"""
    m = cfg["mount"]; dt = 1.0 / m["rate_hz"]; n = int(duration_s / dt); mt = Mount(cfg, np.random.default_rng(seed), az0=100.0, el0=40.0); t = np.arange(n) * dt
    caz = 100.0 + rate_deg_s * t / math.cos(math.radians(40)); ffa = np.full(n, rate_deg_s / math.cos(math.radians(40))) if feedforward else np.zeros(n); cel = np.full(n, 40.0); z = np.zeros(n)
    mt.run(caz[:50], cel[:50], ffa[:50], z[:50], np.full(50, wind_ms)); a, e, da, de, *_ = mt.run(caz[50:], cel[50:], ffa[50:], z[50:], np.full(n - 50, wind_ms))
    follow = np.hypot((a - caz[50:]) * math.cos(math.radians(40)), e - 40.0); return float(np.sqrt(np.mean(follow[100:] ** 2))), float(np.sqrt(np.mean(da ** 2 + de ** 2)))
