"""Virtual receivers: Ku-band total power and L-band (H I) spectrometer. The signal is computed from where the beam REALLY points
(commanded + hidden pointing errors + wind), so tracking and pointing errors show up in the data exactly as they would on a real dish."""
import math, numpy as np, astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from rtpipe.common import ang_sep, helio_correction_kms, lsr_correction_kms, site as rt_site, C_KMS
from rtpipe import sensitivity as S

def beam_resp(az_b, el_b, az_s, el_s, fwhm): return np.exp(-4 * math.log(2) * (ang_sep(az_b, el_b, az_s, el_s) / fwhm) ** 2)

class Scene:
    """the sky as the Ku receiver sees it: list of (kind, track, amplitude K, extra)"""
    def __init__(self): self.src = []
    def add(self, kind, track, amp_k, extra=None): self.src.append((kind, track, float(amp_k), extra)); return self

class Receiver:
    def __init__(self, cfg, rng, t_ref):
        self.cfg, self.rng, self.t_ref = cfg, rng, t_ref; r = cfg["receiver"]; self.r = r; self.loc = rt_site(cfg)
        lam = 2.998e8 / (cfg["ku"]["freq_ghz"] * 1e9); self.F = cfg["ku"]["fwhm_theory_deg"]; self.aeff = cfg["telescope"]["aeff_m2"]
        omega = math.pi * math.radians(r["sun_fwhm_deg"] / 1.1774 * 1.0) ** 2        # unused detail kept simple
        flux_sun = r["sun_flux_sfu"] * 1e4; self.T_sun = S.t_a_from_flux(flux_sun, self.aeff)
        self.Fs = r["sun_fwhm_deg"]; self.F_sun_eff = math.hypot(self.F, self.Fs); self.sun_dilution = self.F ** 2 / (self.F ** 2 + self.Fs ** 2)
        self.m31 = SkyCoord(ra=cfg["targets"]["m31"]["ra_deg"] * u.deg, dec=cfg["targets"]["m31"]["dec_deg"] * u.deg); self.freq = None; self.lo_log = []
        self.F_hi = S.fwhm_deg(cfg["hi"]["rest_mhz"] * 1e6, cfg["telescope"]["dish_m"]); self._shape = None
    def gain(self, t, g0): return g0 * (1 + self.r["drift_per_h"] * (np.asarray(t) - self.t_ref) / 3600.0 + 0.002 * np.sin(2 * np.pi * np.asarray(t) / 5400.0))
    def sky_temp(self, el): return 2.7 + self.r["t_atm_zenith"] / np.sin(np.radians(np.clip(el, 10, 90)))
    # ---------------- Ku total power ----------------
    def ku_temperature(self, t, az_b, el_b, scene):
        T = self.sky_temp(el_b).copy()
        for kind, trk, amp, extra in scene.src:
            az_s, el_s = trk.at(t)
            if kind == "sun": T += amp * self.sun_dilution * np.exp(-4 * math.log(2) * (ang_sep(az_b, el_b, az_s, el_s) / self.F_sun_eff) ** 2)
            elif kind == "sat": T += amp * (extra / trk.range_at(t)) ** 2 * beam_resp(az_b, el_b, az_s, el_s, self.F) * (el_s > 0)
            else: T += amp * beam_resp(az_b, el_b, az_s, el_s, self.F)
        return T
    def ku_counts(self, t, az_b, el_b, scene, cadence_s, dt):
        """integrate to the requested cadence (block mean) and add radiometer noise"""
        k = max(1, int(round(cadence_s / dt))); n = len(t) // k * k
        T = self.ku_temperature(t[:n], az_b[:n], el_b[:n], scene); tc = t[:n].reshape(-1, k).mean(1); Tm = T.reshape(-1, k).mean(1)
        c = self.gain(tc, self.r["gain_ku"]) * (self.r["t_rx_ku"] + Tm) * (1 + self.rng.normal(0, self.r["noise_rel_1s"] / math.sqrt(cadence_s), len(tc)))
        return tc, c, el_b[:n].reshape(-1, k).mean(1)
    def ku_load(self, t, cadence_s, load_k):
        c = self.gain(t, self.r["gain_ku"]) * (self.r["t_rx_ku"] + load_k) * (1 + self.rng.normal(0, self.r["noise_rel_1s"] / math.sqrt(cadence_s), len(t))); return c
    # ---------------- L-band spectrometer ----------------
    def hi_tune(self, t):
        h = self.cfg["hi"]; f0 = h["rest_mhz"] * 1e6; corr = float(helio_correction_kms(self.m31, Time(t, format="unix"), self.loc)); bw = h["bandwidth_mhz"] * 1e6; nch = h["nchan"]
        fc = f0 * (1 - (h["center_vhelio_kms"] - corr) / C_KMS); self.freq = fc + (np.arange(nch) - nch / 2) * bw / nch; self.lo_log.append((float(t), float(fc)))
        self.bp = 1 + 0.08 * np.sin(2 * np.pi * np.arange(nch) / 900) - 0.25 * ((np.arange(nch) - nch / 2) / (nch / 2)) ** 6
        self.rfi = self.rng.choice(nch, 6, replace=False); return self.freq
    def _m31_shape(self):
        if self._shape is None:
            m = self.r["m31"]; g = np.arange(-1200.0, 600.0, 2.0); phi = np.linspace(0, 2 * np.pi, 20000, endpoint=False)
            vlos = m["vsys_helio_kms"] + m["vrot_kms"] * math.sin(math.radians(m["inc_deg"])) * np.cos(phi); h, _ = np.histogram(vlos, bins=np.append(g - 1.0, g[-1] + 1.0))
            k = np.exp(-0.5 * (np.arange(-30, 31) * 2.0 / 12.0) ** 2); k /= k.sum(); p = np.convolve(h, k, mode="same"); self._shape = (g, p / (p.sum() * 2.0))
        return self._shape
    @staticmethod
    def mw_profile(v): return 27.0 * np.exp(-0.5 * ((v + 0.8) / 6.2) ** 2) + 9.3 * np.exp(-0.5 * ((v + 11.5) / 15.5) ** 2)
    def hi_spectra(self, t_c, resp, theta_deg, dump_s):
        """one spectrum per dump: beam response to M31 `resp`, mean angular distance `theta_deg` (drives the Milky Way foreground)"""
        h = self.cfg["hi"]; f0 = h["rest_mhz"] * 1e6; v_topo = C_KMS * (f0 - self.freq) / f0; nch = len(self.freq); bw = h["bandwidth_mhz"] * 1e6
        corr = helio_correction_kms(self.m31, Time(t_c, format="unix"), self.loc); lsr = float(lsr_correction_kms(self.m31)); g, pg = self._m31_shape()
        sc = 1e-26 * self.aeff / (2 * S.K_B) * self.r["m31_flux_jy_kms"]; sig = 1 / math.sqrt(bw / nch * dump_s); out = np.empty((len(t_c), nch), np.float32)
        for i in range(len(t_c)):
            vh = v_topo + corr[i]; vl = vh + lsr; mwf = max(0.3, 1 + 0.01 * (theta_deg[i] - 30.0))
            Tsky = h["t_cold_k"] + mwf * self.mw_profile(vl) + resp[i] * sc * np.interp(vh, g, pg)
            c = self.gain(t_c[i], self.r["gain_hi"]) * self.bp * (self.r["t_rx_hi"] + Tsky) * (1 + self.rng.normal(0, sig, nch)); c[self.rfi] *= 1.6; out[i] = c
        return out
    def hi_hot(self, t_c, dump_s):
        h = self.cfg["hi"]; nch = len(self.freq); bw = h["bandwidth_mhz"] * 1e6; sig = 1 / math.sqrt(bw / nch * dump_s)
        return np.array([self.gain(t, self.r["gain_hi"]) * self.bp * (self.r["t_rx_hi"] + h["t_hot_k"]) * (1 + self.rng.normal(0, sig, nch)) for t in t_c], np.float32)
