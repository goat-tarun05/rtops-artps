import json, warnings
import numpy as np, yaml
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_sun
from astropy.time import Time
from astropy.utils import iers
iers.conf.auto_download = False            # offline: astropy's bundled tables are accurate to << 1 arcmin
iers.conf.auto_max_age = None
warnings.filterwarnings("ignore")
C_KMS = 299792.458

def load_config(path=None):
    import os
    path = path or os.environ.get("RTPIPE_CONFIG", "config.yaml")          # RTPIPE_CONFIG=scenarios/x.yaml python run_all.py --simulate
    with open(path) as f:
        return yaml.safe_load(f)

def site(cfg):
    s = cfg["site"]
    return EarthLocation(lat=s["lat_deg"] * u.deg, lon=s["lon_deg"] * u.deg, height=s["height_m"] * u.m)

def sun_altaz(t, loc):
    return get_sun(t).transform_to(AltAz(obstime=t, location=loc))

def ang_sep(az1, el1, az2, el2):
    """great-circle separation, degrees (all inputs degrees, arrays ok)"""
    a1, e1, a2, e2 = map(np.radians, (az1, el1, az2, el2))
    c = np.sin(e1) * np.sin(e2) + np.cos(e1) * np.cos(e2) * np.cos(a1 - a2)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))

def offset_pointing(az, el, dx_sky, dy):
    """apply a sky-plane offset (dx along azimuth on the sky, dy in elevation), degrees"""
    return az + dx_sky / np.cos(np.radians(el)), el + dy

def helio_correction_kms(coord, times, loc):
    """km/s to ADD to a topocentric radial velocity to get the barycentric (~heliocentric) one"""
    return coord.radial_velocity_correction(kind="barycentric", obstime=times, location=loc).to(u.km / u.s).value

_APEX = SkyCoord(ra=270.97 * u.deg, dec=30.0 * u.deg)          # standard solar motion, 20 km/s
def lsr_correction_kms(coord):
    return 20.0 * np.cos(coord.separation(_APEX).radian)

def save_meta(path, meta):
    with open(path, "w") as f:
        json.dump(meta, f, indent=1, default=float)

def load_meta(path):
    with open(path) as f:
        return json.load(f)
