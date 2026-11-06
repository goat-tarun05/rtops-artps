"""Software tracking of the beam: every sample's fixed alt/az (corrected by the fitted pointing model) -> RA/Dec and l/b."""
import os, numpy as np, astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, get_sun
from astropy.table import Table
from .common import load_meta, site, offset_pointing

def tag_fixed(cfg, unix, az, el, value, name, out, step=10):
    loc = site(cfg); i = np.arange(0, len(unix), step); t = Time(unix[i], format="unix")
    c = SkyCoord(az=np.full(len(i), az) * u.deg, alt=np.full(len(i), el) * u.deg, frame=AltAz(obstime=t, location=loc))
    icrs, gal = c.icrs, c.galactic
    tab = Table({"utc": t.isot, "az_deg": np.full(len(i), az), "el_deg": np.full(len(i), el), "ra_deg": icrs.ra.deg, "dec_deg": icrs.dec.deg,
                 "l_deg": gal.l.deg, "b_deg": gal.b.deg, "value": value[i]})
    os.makedirs(f"{out}/tables", exist_ok=True); tab.write(f"{out}/tables/{name}.ecsv", overwrite=True)
    return tab

def run(cfg, raw="data/raw", out="outputs", solar_res=None):
    res = {}; dx = solar_res["pointing_dx_sky_deg"] if solar_res else 0.0; dy = solar_res["pointing_dy_deg"] if solar_res else 0.0
    res["pointing_model_applied_deg"] = [dx, dy]
    for k in (1, 2):
        d = np.load(f"{raw}/solar_transit_{k}.npz"); m = load_meta(f"{raw}/solar_transit_{k}.json")
        az, el = offset_pointing(m["az_deg"], m["el_deg"], dx, dy)                           # true beam centre = commanded + fitted offsets
        tab = tag_fixed(cfg, d["unix"], az, el, d["counts"], f"solar_scan{k}_tagged", out)
        if solar_res:
            tpk = Time(solar_res["transit_obs_utc"][k - 1]); j = int(np.argmin(np.abs(Time(tab["utc"]).unix - tpk.unix)))
            sun = get_sun(tpk); beam = SkyCoord(ra=tab["ra_deg"][j] * u.deg, dec=tab["dec_deg"][j] * u.deg)
            res[f"scan{k}_beam_to_sun_deg"] = float(beam.separation(SkyCoord(ra=sun.ra, dec=sun.dec)).deg)
            res[f"scan{k}_commanded_offset_deg"] = abs(float(m["el_offset_set_deg"]))
    d = np.load(f"{raw}/starlink_scan.npz"); m = load_meta(f"{raw}/starlink_scan.json")
    az, el = offset_pointing(m["az_deg"], m["el_deg"], dx, dy)
    tag_fixed(cfg, d["unix"], az, el, d["counts"], "starlink_scan_tagged", out, step=100)
    res["tables"] = sorted(os.listdir(f"{out}/tables"))
    return res
