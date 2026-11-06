"""Alt-az mount servo simulation: rate/acceleration limits, position loop with feed-forward, actuator lag, encoder quantisation,
wind-induced pointing jitter. Runs at a fixed control rate (default 10 Hz)."""
import math, numpy as np

class Mount:
    def __init__(self, cfg, rng, az0=0.0, el0=60.0):
        m = cfg["mount"]; self.m = m; self.rng = rng; self.dt = 1.0 / m["rate_hz"]
        self.az, self.el, self.vaz, self.vel, self.d_az, self.d_el = az0, el0, 0.0, 0.0, 0.0, 0.0
    def run(self, caz, cel, ffaz, ffel, wind):
        """follow the commanded encoder angles; returns encoder az/el (quantised), wind deflection az/el (sky-plane deg), true axis az/el"""
        m = self.m; dt = self.dt; kp, acc, tau = m["kp"], m["accel"], m["tau_act"]; vmz, vme = m["vmax_az"], m["vmax_el"]
        azl, azh, ell, elh = m["az_range"][0], m["az_range"][1], m["el_range"][0], m["el_range"][1]
        n = len(caz); out = np.empty((n, 4)); noise = self.rng.normal(size=(n, 3)); sig = m["wind_sigma_per_ms"]; tw = m["wind_tau_s"]; k_ou = math.sqrt(2 * dt / tw)
        az, el, vaz, vel, da, de = self.az, self.el, self.vaz, self.vel, self.d_az, self.d_el; sq = math.sqrt; ad = acc * dt; g = dt / tau
        for i in range(n):
            out[i, 0], out[i, 1], out[i, 2], out[i, 3] = az, el, da, de                    # state AT the sample time; the command for this sample then drives the step to the next one
            e = caz[i] - az; vp = math.copysign(min(kp * abs(e), sq(2 * acc * abs(e))), e); vd = max(-vmz, min(vmz, ffaz[i] + vp))
            vaz += max(-ad, min(ad, (vd - vaz) * g)); az = max(azl, min(azh, az + vaz * dt))
            e = cel[i] - el; vp = math.copysign(min(kp * abs(e), sq(2 * acc * abs(e))), e); vd = max(-vme, min(vme, ffel[i] + vp))
            vel += max(-ad, min(ad, (vd - vel) * g)); el = max(ell, min(elh, el + vel * dt))
            s = sig * wind[i] * k_ou; da += -da * dt / tw + s * noise[i, 0]; de += -de * dt / tw + s * noise[i, 1]
        self.az, self.el, self.vaz, self.vel, self.d_az, self.d_el = az, el, vaz, vel, da, de
        res, en = m["enc_res_deg"], m["enc_noise_deg"]
        return (np.round((out[:, 0] + en * noise[:, 2]) / res) * res, np.round((out[:, 1] + en * noise[::-1, 2]) / res) * res, out[:, 2], out[:, 3], out[:, 0], out[:, 1])
