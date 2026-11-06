"""Solar transit: Y-factor calibration, predicted solar track, beam + pointing fit (celestial tracking)."""
import numpy as np, astropy.units as u
from astropy.time import Time
from scipy.optimize import curve_fit
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from .common import load_meta, site, sun_altaz, ang_sep, offset_pointing

def _load_scan(raw, k, cfg, loc, mask=None):
    d = np.load(f"{raw}/solar_transit_{k}.npz"); m = load_meta(f"{raw}/solar_transit_{k}.json"); ku = cfg["ku"]
    unix, counts = d["unix"], d["counts"]; n0 = len(unix)
    if mask:                                                                                 # Starlink passes contaminate the Ku-band Sun scan
        keep = np.ones(n0, bool)
        for a, b in mask: keep &= ~((unix >= Time(a).unix) & (unix <= Time(b).unix))
        unix, counts = unix[keep], counts[keep]
    Y = d["hot"].mean() / d["cold"].mean(); t_rx = (ku["t_hot_k"] - Y * ku["t_cold_k"]) / (Y - 1)
    g = d["cold"].mean() / (t_rx + ku["t_cold_k"]); T = counts / g                          # Y-factor -> kelvin
    mins = (unix - unix[0]) / 60
    edge = (mins < 8) | (mins > mins[-1] - 8)                                                 # Sun outside the beam
    dT = T - np.polyval(np.polyfit(mins[edge], T[edge], 1), mins)                             # remove gain drift
    s = sun_altaz(Time(unix, format="unix"), loc)                                            # predicted solar track per sample
    return dict(unix=unix, mins=mins, dT=dT, saz=s.az.deg, sel=s.alt.deg, m=m, Y=Y, t_rx=t_rx, g=g, removed=n0 - len(unix))

def run(cfg, raw="data/raw", out="outputs", rfi_mask=None):
    loc = site(cfg); ku = cfg["ku"]
    sc = [_load_scan(raw, k, cfg, loc, rfi_mask) for k in (1, 2)]
    def beam(s, A, F, dx, dy):
        az, el = offset_pointing(s["m"]["az_deg"], s["m"]["el_deg"], dx, dy)
        return A * np.exp(-4 * np.log(2) * ang_sep(az, el, s["saz"], s["sel"]) ** 2 / F ** 2)
    n1 = len(sc[0]["dT"])
    def model(i, A, F, dx, dy):
        return np.concatenate([beam(sc[0], A, F, dx, dy), beam(sc[1], A, F, dx, dy)])
    y = np.concatenate([s["dT"] for s in sc]); idx = np.arange(len(y))
    p, cov = curve_fit(model, idx, y, p0=[y.max() * 1.5, ku["fwhm_theory_deg"], 0.0, 0.0], maxfev=20000)
    err = np.sqrt(np.diag(cov)); A, F, dx, dy = p; rms = float(np.sqrt(np.mean((y - model(idx, *p)) ** 2)))
    tim = []
    for s in sc:
        sep_nom = ang_sep(s["m"]["az_deg"], s["m"]["el_deg"], s["saz"], s["sel"])
        tim.append((float(s["unix"][np.argmin(sep_nom)]), float(s["unix"][np.argmax(beam(s, *p))])))
    res = dict(y_factor=float(sc[0]["Y"]), t_rx_k=float(sc[0]["t_rx"]), t_sys_k=float(sc[0]["t_rx"] + ku["t_cold_k"]),
               sun_peak_k=float(A), fwhm_deg=float(F), fwhm_err_deg=float(err[1]), fwhm_theory_deg=ku["fwhm_theory_deg"],
               pointing_dx_sky_deg=float(dx), pointing_dy_deg=float(dy), pointing_err_deg=float(np.hypot(err[2], err[3])),
               transit_pred_utc=[Time(a, format="unix").isot for a, b in tim], transit_obs_utc=[Time(b, format="unix").isot for a, b in tim],
               transit_offset_s=[b - a for a, b in tim], fit_rms_k=rms, masked_samples=int(sum(x["removed"] for x in sc)), pointing_systematic_floor_deg=0.3)
    res["checks"] = {"beam within 20% of 70 lambda/D": bool(abs(F / ku["fwhm_theory_deg"] - 1) < 0.2),
                     "fit residual < 2% of Sun peak": bool(rms < 0.02 * A),
                     "pointing fit uncertainty < 0.1 deg (statistical; manual setting adds ~0.3 deg)": bool(res["pointing_err_deg"] < 0.1)}
    fig, ax = plt.subplots(2, 1, figsize=(8, 6.5))
    for k, (s, c) in enumerate(zip(sc, ("tab:blue", "tab:green"))):
        ax[0].plot(s["mins"], s["dT"], ".", ms=1.5, color=c, label=f"scan {k+1}: dish {s['m']['el_offset_set_deg']:+.0f}° in el")
        ax[0].plot(s["mins"], beam(s, *p), "r-", lw=1)
        ax[0].axvline((tim[k][0] - s["unix"][0]) / 60, color=c, ls="--", lw=0.8)
    ax[0].set_xlabel("minutes from scan start"); ax[0].set_ylabel("K above baseline"); ax[0].legend(fontsize=8)
    ax[0].set_title(f"Solar cross-scan: FWHM {F:.2f}° (theory {ku['fwhm_theory_deg']}°), pointing ({dx:+.2f}°, {dy:+.2f}°)")
    th = np.linspace(0, 2 * np.pi, 200); m0 = sc[0]["m"]
    for k, (s, c) in enumerate(zip(sc, ("tab:blue", "tab:green"))):
        xs = (s["saz"] - s["m"]["az_deg"]) * np.cos(np.radians(s["m"]["el_deg"])); ys = s["sel"] - s["m"]["el_deg"]
        ax[1].plot(xs, ys, color=c, lw=1, label=f"Sun track rel. to dish, scan {k+1}")
    ax[1].plot(dx + F / 2 * np.cos(th), dy + F / 2 * np.sin(th), "r-", label="fitted beam (FWHM)")
    ax[1].plot(F / 2 * np.cos(th), F / 2 * np.sin(th), "k:", label="where we thought it pointed")
    ax[1].set_aspect("equal"); ax[1].set_xlim(-8, 8); ax[1].set_ylim(-4, 4); ax[1].legend(fontsize=7)
    ax[1].set_xlabel("sky offset in azimuth (deg)"); ax[1].set_ylabel("elevation offset (deg)")
    fig.tight_layout(); fig.savefig(f"{out}/solar_transit.png", dpi=130); plt.close(fig)
    return res
