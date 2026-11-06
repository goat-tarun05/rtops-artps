"""Ephemeris service: apparent az/el tracks (with atmospheric refraction) for sidereal, solar, GEO-beacon and satellite targets."""
import math, numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, get_sun
from rtpipe.common import site as rt_site, ang_sep
from rtpipe.satellites import satrec_from_elements
from skyfield.api import load, EarthSatellite, wgs84

def refraction_deg(el_true):
    """atmospheric refraction (Saemundsson), degrees, to ADD to the true elevation"""
    h = np.clip(np.asarray(el_true, float), 1.0, 90.0)
    return (1.02 / np.tan(np.radians(h + 10.3 / (h + 5.11)))) / 60.0

class Track:
    """tabulated apparent position of a target; at(t) interpolates (az mod 360, el)"""
    def __init__(self, name, t, az, el, rng_km=None):
        self.name = name; self.t = np.asarray(t, float); self.el = np.asarray(el, float)
        self.azu = np.degrees(np.unwrap(np.radians(np.asarray(az, float)))); self.rng = None if rng_km is None else np.asarray(rng_km, float)
    def at(self, t):
        t = np.asarray(t, float); return np.interp(t, self.t, self.azu) % 360.0, np.interp(t, self.t, self.el)
    def range_at(self, t): return np.interp(np.asarray(t, float), self.t, self.rng)

class FixedTrack:
    def __init__(self, name, az, el): self.name, self.az, self.el = name, float(az), float(el)
    def at(self, t):
        t = np.asarray(t, float); return np.full(t.shape, self.az), np.full(t.shape, self.el)

def geo_azel(lat, lon, h_m, sat_lon):
    la, lo, sl = map(math.radians, (lat, lon, sat_lon)); a, f = 6378.137, 1 / 298.257; e2 = f * (2 - f); h = h_m / 1000.0
    N = a / math.sqrt(1 - e2 * math.sin(la) ** 2)
    obs = np.array([(N + h) * math.cos(la) * math.cos(lo), (N + h) * math.cos(la) * math.sin(lo), (N * (1 - e2) + h) * math.sin(la)])
    sat = np.array([42164.17 * math.cos(sl), 42164.17 * math.sin(sl), 0.0]); d = sat - obs
    e = np.array([-math.sin(lo), math.cos(lo), 0]); n = np.array([-math.sin(la) * math.cos(lo), -math.sin(la) * math.sin(lo), math.cos(la)])
    up = np.array([math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la)])
    E, Nn, U = d @ e, d @ n, d @ up
    return math.degrees(math.atan2(E, Nn)) % 360, math.degrees(math.atan2(U, math.hypot(E, Nn)))

def shell_elements(cfg, epoch_jd):
    tr = cfg["satellites"]; mu, re = 398600.4418, 6378.135; a = re + tr["alt_km"]; n = math.sqrt(mu / a ** 3) * 60.0; els = []; k = 0
    for p in range(tr["planes"]):
        for s in range(tr["per_plane"]):
            k += 1
            els.append(dict(satnum=70000 + k, name=f"STARLINK-SIM-{k:04d}", inc=tr["inc_deg"], raan=360.0 * p / tr["planes"], argp=0.0, ecc=1e-4,
                            ma=(360.0 * s / tr["per_plane"] + 360.0 * p / (tr["planes"] * tr["per_plane"]) * 17) % 360, n_rad_min=n, epoch_jd=epoch_jd))
    return els

class Ephemeris:
    def __init__(self, cfg):
        self.cfg = cfg; self.loc = rt_site(cfg); s = cfg["site"]; self.lat, self.lon, self.h = s["lat_deg"], s["lon_deg"], s["height_m"]
        self._ts = load.timescale(builtin=True); self._obs = wgs84.latlon(self.lat, self.lon, self.h)
    def _frame(self, t): return AltAz(obstime=Time(t, format="unix"), location=self.loc)
    def radec(self, name, ra, dec, t0, t1, step=10.0):
        t = np.arange(t0, t1 + step, step); aa = SkyCoord(ra=ra * u.deg, dec=dec * u.deg).transform_to(self._frame(t))
        return Track(name, t, aa.az.deg, aa.alt.deg + refraction_deg(aa.alt.deg))
    def sun(self, t0, t1, step=10.0):
        t = np.arange(t0, t1 + step, step); aa = get_sun(Time(t, format="unix")).transform_to(self._frame(t))
        return Track("Sun", t, aa.az.deg, aa.alt.deg + refraction_deg(aa.alt.deg))
    def geo(self, sat_lon):
        az, el = geo_azel(self.lat, self.lon, self.h, sat_lon); return FixedTrack(f"GEO{sat_lon:+.0f}E", az, el + float(refraction_deg(el)))
    def satellite(self, name, satrec, t0, t1, step=1.0):
        t = np.arange(t0, t1 + step, step); alt, az, dist = (EarthSatellite.from_satrec(satrec, self._ts) - self._obs).at(self._ts.from_astropy(Time(t, format="unix"))).altaz()
        el = alt.degrees; return Track(name, t, az.degrees, el + refraction_deg(el), dist.km)
    def find_passes(self, elements, t0, t1, step=20.0, min_el=25.0):
        """coarse pass search over a catalogue: max elevation, AOS/LOS, peak az rate"""
        t = np.arange(t0, t1, step); tsf = self._ts.from_astropy(Time(t, format="unix")); out = []
        for e in elements:
            alt, az, dist = (EarthSatellite.from_satrec(satrec_from_elements(e), self._ts) - self._obs).at(tsf).altaz(); el = alt.degrees
            up = el > min_el
            if not up.any(): continue
            idx = np.where(np.diff(np.r_[0, up.astype(int), 0]) != 0)[0]
            for a, b in zip(idx[::2], idx[1::2]):
                if b - a < 3: continue
                seg = slice(a, b); azs = np.unwrap(np.radians(az.degrees[seg])); rate = np.degrees(np.abs(np.gradient(azs, step))) * np.cos(np.radians(el[seg]))
                out.append(dict(name=e["name"], t_aos=float(t[a]), t_los=float(t[b - 1]), max_el=float(el[seg].max()), max_az_rate=float(rate.max()), min_range_km=float(dist.km[seg].min()), az_span=float(np.degrees(np.ptp(azs)))))
        return out
