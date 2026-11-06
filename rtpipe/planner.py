"""Observation planner: Sun noon pointing (+- 1 deg cross-scan), M31 transit pointing and hours above 40 deg."""
import csv, os, numpy as np, astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, get_sun
from .common import site
from . import sensitivity as S

def plan(cfg, start=None, days=14):
    loc = site(cfg); rows = []
    t0 = Time((start or Time.now().iso[:10]) + "T00:00:00")
    m31 = SkyCoord(ra=cfg["targets"]["m31"]["ra_deg"] * u.deg, dec=cfg["targets"]["m31"]["dec_deg"] * u.deg)
    fhi = S.fwhm_deg(1420.40575e6, cfg["hi"]["dish_m"])
    for d in range(days):
        tt = t0 + d * u.day + np.arange(0, 1440) * u.min
        fr = AltAz(obstime=tt, location=loc)
        sun = get_sun(tt).transform_to(fr); i = int(np.argmax(sun.alt.deg))
        mm = m31.transform_to(fr); j = int(np.argmax(mm.alt.deg))
        rows.append(dict(date=tt[0].iso[:10], sun_noon_utc=tt[i].iso[11:16], sun_az=round(float(sun.az.deg[i]), 1), sun_el=round(float(sun.alt.deg[i]), 1),
                         scan_el_plus=round(float(sun.alt.deg[i]) + 1, 1), scan_el_minus=round(float(sun.alt.deg[i]) - 1, 1),
                         m31_transit_utc=tt[j].iso[11:16], m31_az=round(float(mm.az.deg[j]), 1), m31_el=round(float(mm.alt.deg[j]), 1),
                         m31_hours_above_40=round(float(np.sum(mm.alt.deg > 40)) / 60, 1), hi_beam_fwhm_deg=round(fhi, 1)))
    return rows

def main(cfg, days=14, out="outputs"):
    os.makedirs(out, exist_ok=True); rows = plan(cfg, days=days)
    with open(f"{out}/plan.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print("date        Sun noon(UTC) az    el   | cross-scan el  | M31 transit(UTC) az   el   h>40")
    for r in rows:
        print(f"{r['date']}  {r['sun_noon_utc']}        {r['sun_az']:5.1f} {r['sun_el']:4.1f} | {r['scan_el_plus']:4.1f} / {r['scan_el_minus']:4.1f}  | {r['m31_transit_utc']}            {r['m31_az']:5.1f} {r['m31_el']:4.1f} {r['m31_hours_above_40']}")
    print(f"written {out}/plan.csv   (HI beam {rows[0]['hi_beam_fwhm_deg']} deg: for 'switch' mode the OFF position must be >= {1.5*rows[0]['hi_beam_fwhm_deg']:.0f} deg away)")
    return rows
