"""Andromeda HI.  Two ways to get an OFF reference (config hi.mode):
  drift  : dish FIXED. The pipeline computes the beam response of every spectrum to M31 (celestial tracking);
           high response = ON, ~zero response = OFF; ON spectra are combined with beam-response weights.
  switch : dish moved ON/OFF; the OFF position's leakage of M31 (it may sit inside the beam) is corrected.
Doppler tracking: every ON spectrum is bin-averaged onto a heliocentric velocity grid at its own timestamp."""
import numpy as np, astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz
from astropy.table import Table
from scipy.ndimage import median_filter, uniform_filter1d
from scipy.optimize import least_squares, minimize_scalar
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from .common import load_meta, site, helio_correction_kms, lsr_correction_kms, C_KMS
from . import sensitivity as S

VG = np.arange(-900.0, 300.0, 1.0)
FREE = ((VG > -880) & (VG < -620)) | ((VG > 80) & (VG < 280))

def _rebin(v, y):
    """bin-average native channels onto the 1 km/s grid (np.interp would throw away most of the sensitivity)"""
    idx = np.floor(v - (VG[0] - 0.5)).astype(int); ok = (idx >= 0) & (idx < len(VG))
    s = np.bincount(idx[ok], weights=y[ok], minlength=len(VG)); n = np.bincount(idx[ok], minlength=len(VG))
    out = np.full(len(VG), np.nan); m = n > 0; out[m] = s[m] / n[m]
    bad = np.isnan(out)
    if bad.any(): out[bad] = np.interp(VG[bad], VG[~bad], out[~bad])
    return out - np.median(out[FREE])                                                   # per-dump DC (gain-drift) removal

_G = np.arange(-1500.0, 800.0, 1.0)
def template(v, vsys, vrot, sm=12.0):
    """rotating-disc double-horn profile, unit integral per km/s, Gaussian-smoothed (sigma sm km/s); exact bin masses from the arcsine CDF"""
    e = np.arange(_G[0] - 0.5, _G[-1] + 1.0, 1.0); x = np.clip((e - vsys) / vrot, -1, 1)
    mass = np.diff((np.arcsin(x) + np.pi / 2) / np.pi); k = np.exp(-0.5 * (np.arange(-40, 41) / sm) ** 2); k /= k.sum()
    return np.interp(v, _G, np.convolve(mass, k, mode="same"))

def _vp(v, y, sigma, prior, vs, vr):
    """chi2 at fixed (vsys, vrot): amplitude and quadratic baseline are linear, so they are solved exactly (variable projection)"""
    x = (v + 300.0) / 600.0; A = np.column_stack([template(v, vs, vr), np.ones_like(x), x, x ** 2]); coef = np.linalg.lstsq(A, y, rcond=None)[0]; r = (y - A @ coef) / sigma
    return float(r @ r) + ((vr - prior[0]) / prior[1]) ** 2, coef

def _best_vr(v, y, sigma, prior, vs):
    g = np.arange(80.0, 401.0, 20.0); c = [_vp(v, y, sigma, prior, vs, q)[0] for q in g]; q0 = g[int(np.argmin(c))]
    r = minimize_scalar(lambda q: _vp(v, y, sigma, prior, vs, q)[0], bounds=(max(80.0, q0 - 20), min(400.0, q0 + 20)), method="bounded", options={"xatol": 0.3}); return float(r.x), float(r.fun)

def fit_template(v, y, sigma, prior):
    """joint fit of  I x template(v; vsys, vrot) + quadratic baseline  with a Gaussian prior on vrot (the receding edge is blended with the
    Milky Way, so the data alone cannot separate v_sys from the width). The v_sys error comes from the PROFILE likelihood (global minimum over
    v_sys, error = the region where chi2 rises by less than 1 x chi2/dof), not from the local curvature: with S/N of a few the chi2 surface has
    several minima and the local curvature can be wildly optimistic.  Returns (I, vsys, vrot), their 1-sigma errors, sqrt(chi2/dof), baseline."""
    grid = np.arange(-440.0, -159.0, 5.0); prof = np.array([_best_vr(v, y, sigma, prior, g)[1] for g in grid]); c2min = float(prof.min()); dof = len(v) + 1 - 6; s2 = max(1.0, c2min / dof)
    w = np.exp(-(prof - c2min) / (2 * s2)); w /= w.sum(); vs = float(np.sum(w * grid)); vs_err = max(float(np.sqrt(np.sum(w * (grid - vs) ** 2))), 2.5)       # posterior mean and standard deviation over v_sys
    vr, c2 = _best_vr(v, y, sigma, prior, vs); coef = _vp(v, y, sigma, prior, vs, vr)[1]
    def fun(q): return np.append((y - q[0] * template(v, q[1], q[2]) - q[3] - q[4] * (v + 300.0) / 600.0 - q[5] * ((v + 300.0) / 600.0) ** 2) / sigma, (q[2] - prior[0]) / prior[1])
    q0 = [coef[0], vs, vr, coef[1], coef[2], coef[3]]
    try:
        lq = least_squares(fun, q0, x_scale=[abs(coef[0]) + 1e-3, 10, 10, 0.05, 0.05, 0.05], bounds=([-np.inf, -460, 80, -np.inf, -np.inf, -np.inf], [np.inf, -140, 400, np.inf, np.inf, np.inf]))
        cov = np.linalg.pinv(lq.jac.T @ lq.jac) * s2; I_err, vr_err = float(np.sqrt(cov[0, 0])), float(np.sqrt(cov[2, 2]))
    except Exception: I_err, vr_err = float(abs(coef[0]) / 5), 25.0
    return np.array([coef[0], vs, vr]), np.array([I_err, vs_err, vr_err]), float(np.sqrt(s2)), np.array(coef[1:])

def run(cfg, raw="data/raw", out="outputs", rfi_mask=None):
    h = cfg["hi"]; loc = site(cfg)
    d = np.load(f"{raw}/hi_m31.npz"); m = load_meta(f"{raw}/hi_m31.json")
    f0 = m["rest_mhz"] * 1e6; v_topo = C_KMS * (f0 - d["freq_hz"]) / f0
    tgt = SkyCoord(ra=m["target"]["ra_deg"] * u.deg, dec=m["target"]["dec_deg"] * u.deg)
    unix = d["unix"]; spec = d["spec"].astype(float); hot = d["hot"].astype(float); t = Time(unix, format="unix")
    F = S.fwhm_deg(f0, h["dish_m"]); mode = m["mode"]; leak = 0.0
    if mode == "drift":
        beam = SkyCoord(az=np.full(len(t), m["az_deg"]) * u.deg, alt=np.full(len(t), m["el_deg"]) * u.deg, frame=AltAz(obstime=t, location=loc)).icrs
        theta = beam.separation(tgt).deg; r = np.exp(-4 * np.log(2) * (theta / F) ** 2)
        near = np.abs(unix - unix[np.argmax(r)]) < h["off_window_hours"] * 3600                  # OFF close in time = similar Milky Way foreground
        is_on = r > h["on_response_min"]; is_off = (r < h["off_response_max"]) & near; off_sep = None
    else:
        is_on = d["is_on"].astype(bool); is_off = ~is_on
        off_sep = float(tgt.separation(SkyCoord(ra=m["off"]["ra_deg"] * u.deg, dec=m["off"]["dec_deg"] * u.deg)).deg)
        leak = float(np.exp(-4 * np.log(2) * (off_sep / F) ** 2)); r = np.where(is_on, 1.0, leak)
    if rfi_mask:
        keep = np.ones(len(unix), bool)
        for a, b in rfi_mask: keep &= ~((unix >= Time(a).unix) & (unix <= Time(b).unix))
        is_on &= keep; is_off &= keep
    n_on, n_off = int(is_on.sum()), int(is_off.sum())
    if n_on == 0 or n_off == 0: raise ValueError(f"need ON and OFF spectra (found {n_on} ON, {n_off} OFF); check hi.mode, pointing and the observing window")
    # 1. per-channel Y-factor -> Tsys; 2. ON-OFF antenna temperature
    off_ref = spec[is_off].mean(0); Y = hot.mean(0) / off_ref
    t_sys = median_filter((h["t_hot_k"] - Y * h["t_cold_k"]) / (Y - 1) + h["t_cold_k"], 101)
    ta = (spec[is_on] - off_ref) / off_ref * t_sys
    # 3. RFI: narrow spikes vs a running median, flagged and interpolated
    mean_ta = ta.mean(0); resid = mean_ta - median_filter(mean_ta, 31)
    bad = np.convolve(np.abs(resid) > 6 * 1.4826 * np.median(np.abs(resid)), np.ones(5), mode="same") > 0; ch = np.arange(ta.shape[1])
    for i in range(n_on): ta[i, bad] = np.interp(ch[bad], ch[~bad], ta[i, ~bad])
    # 4. Doppler tracking + beam-weighted combination: T_peak(v) = sum r_i T_i / sum r_i^2
    corr = helio_correction_kms(tgt, t[is_on], loc); lsr = float(lsr_correction_kms(tgt)); r_on = r[is_on]
    num = np.zeros(len(VG)); num0 = np.zeros(len(VG))
    for k in range(n_on):
        num += r_on[k] * _rebin(v_topo + corr[k], ta[k]); num0 += r_on[k] * _rebin(v_topo, ta[k])
    den = float((r_on ** 2).sum()); norm = den * ((1 - leak) if h.get("correct_off_leakage", True) else 1.0)
    spec_h = num / norm; nodop = num0 / norm
    # 5. joint line + baseline fit on the part of the band NOT blended with the Milky Way
    cub = np.polyval(np.polyfit(VG[FREE], spec_h[FREE], 3), VG); rms1 = float(np.std((spec_h - cub)[FREE]))      # noise scale for the weights
    mw_lo, mw_hi = h["mw_mask_kms"]; fitw = ((VG > -880) & (VG < mw_lo)) | ((VG > mw_hi) & (VG < 280))
    (I_fit, vsys, vrot), (I_err, vsys_err, vrot_err), chi, bc = fit_template(VG[fitw], spec_h[fitw], rms1, h["vrot_prior_kms"])
    xb = (VG + 300.0) / 600.0; spec_b = spec_h - (bc[0] + bc[1] * xb + bc[2] * xb ** 2)
    sm = uniform_filter1d(spec_b, 10); rms = float(np.std(sm[FREE])); model = I_fit * template(VG, vsys, vrot)
    mw = (VG > -80) & (VG < 80); v_mw_lsr = float(VG[mw][np.argmax(np.abs(spec_h[mw]))] + lsr)
    # 6. radiometer prediction of the noise, and the literature-flux prediction of the signal
    dump = float(np.median(np.diff(unix))); tsys_med = float(np.median(t_sys))
    sig_pred = S.sigma_pred(tsys_med, 10.0, dump, float(r_on.sum()), den, n_off) / ((1 - leak) if h.get("correct_off_leakage", True) else 1.0)
    int_exp = S.t_a_from_flux(S.M31_SINT_JY_KMS, h["aeff_m2"])
    res = dict(mode=mode, beam_fwhm_deg=F, n_on=n_on, n_off=n_off, effective_on_hours=float(den * dump / 3600), off_separation_deg=off_sep, off_leakage=leak,
               vsys_helio_kms=float(vsys), vsys_err_kms=float(vsys_err), vsys_lsr_kms=float(vsys + lsr), vrot_sini_kms=float(vrot), vrot_err_kms=float(vrot_err), width_kms=float(2 * vrot),
               rms_k=rms, rms_predicted_k=float(sig_pred), rms_ratio=float(rms / sig_pred),
               integrated_k_kms=float(I_fit), integrated_err_k_kms=float(I_err), int_ul3_k_kms=float(3 * I_err), int_ul3_jy_kms=float(3 * I_err * 2 * S.K_B / (h["aeff_m2"] * 1e-26)), integrated_expected_k_kms=float(int_exp), integrated_ratio=float(I_fit / int_exp),
               line_snr=float(I_fit / I_err), ideal_line_snr=float(int_exp / (sig_pred * np.sqrt(S.M31_W_KMS * 10.0))), fit_chi_per_dof=chi, tsys_k=tsys_med,
               fit_ranges_kms=[[-880.0, float(mw_lo)], [float(mw_hi), 280.0]], doppler_corr_range_kms=[float(corr.min()), float(corr.max())], rfi_channels_flagged=int(bad.sum()),
               milky_way_residual_lsr_kms=v_mw_lsr, literature_vsys_helio_kms=-300.0)
    res["checks"] = {"systemic velocity within 25 km/s of -300 (literature)": bool(abs(vsys + 300) < 25),
                     "M31 line detected: integrated S/N >= 5": bool(res["line_snr"] >= 5),
                     "measured noise within 40 % of the radiometer prediction": bool(0.6 < res["rms_ratio"] < 1.4),
                     "Milky Way residual found near 0 km/s LSR and kept apart from M31": bool(abs(v_mw_lsr) < 20 and abs(vsys - v_mw_lsr) > 200)}
    if mode == "switch": res["checks"]["OFF position >= 1.5 beam widths from M31 (else leakage corrected)"] = bool(off_sep >= 1.5 * F)
    Table({"v_helio_kms": VG, "v_lsr_kms": VG + lsr, "ta_k": spec_b, "ta_smoothed_k": sm}).write(f"{out}/m31_spectrum.ecsv", overwrite=True)
    fig, ax = plt.subplots(3, 1, figsize=(8, 9.5)); hrs = (unix - unix[len(unix) // 2]) / 3600
    ax[0].plot(hrs, r, "k", lw=1); ax[0].fill_between(hrs, 0, 1, where=is_on, color="tab:orange", alpha=0.35, label=f"ON ({n_on})"); ax[0].fill_between(hrs, 0, 1, where=is_off, color="tab:blue", alpha=0.2, label=f"OFF ({n_off})")
    ax[0].set_xlabel("hours from transit"); ax[0].set_ylabel("beam response to M31"); ax[0].legend(fontsize=8); ax[0].set_title(f"{mode} mode, beam {F:.1f}°: ON/OFF chosen from the computed beam response")
    order = np.argsort(v_topo); vt = v_topo[order]
    ax[1].plot(vt, (off_ref / off_ref.mean())[order], lw=0.6, label="OFF bandpass"); ax[1].plot(vt, (spec[is_on].mean(0) / off_ref.mean())[order], lw=0.6, label="ON bandpass")
    ax[1].set_xlabel("topocentric radio velocity (km/s)"); ax[1].legend(fontsize=8); ax[1].set_title("raw spectra: bandpass ripple and Milky Way dominate, M31 invisible")
    ax[2].plot(VG, sm, "k", lw=0.8, label="ON-OFF, Doppler-tracked, beam-weighted, calibrated")
    ax[2].plot(VG, uniform_filter1d(nodop - np.polyval(np.polyfit(VG[FREE], nodop[FREE], 3), VG), 10), color="0.6", lw=0.6, label="same without Doppler tracking")
    ax[2].plot(VG, model, "r--", lw=1.2, label=f"disc template fit: v_sys = {vsys:.1f} ± {vsys_err:.1f} km/s (lit. -300)")
    ax[2].axvspan(mw_lo, mw_hi, color="orange", alpha=0.2, label="blended with Milky Way: not fitted")
    ax[2].set_xlim(-750, 150); ax[2].set_ylim(-4 * rms - 0.01, max(0.08, float(model.max()) * 2.0)); ax[2].set_xlabel("heliocentric radio velocity (km/s)"); ax[2].set_ylabel("T_A at beam centre (K)")
    ax[2].legend(fontsize=7, loc="upper left"); ax[2].set_title(f"Andromeda HI: width 2V sin i = {2 * vrot:.0f} km/s, integrated S/N {res['line_snr']:.1f}")
    fig.tight_layout(); fig.savefig(f"{out}/hi_m31.png", dpi=130); plt.close(fig)
    return res
