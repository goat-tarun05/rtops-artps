"""Expected signal and noise BEFORE you observe: is the detection feasible with your antenna and time?
    python -m rtpipe.sensitivity --aeff 0.35 --tsys 100 --hours-on 1 --hours-off 4
"""
import math, argparse
K_B = 1.380649e-23; C = 2.998e8; F_HI = 1420.40575e6
M31_SINT_JY_KMS = 31400.0          # integrated HI flux: 4.5e9 Msun at 0.78 Mpc via M = 2.36e5 D^2 S
M31_W_KMS = 470.0                  # line width used for the integrated-S/N estimate

def fwhm_deg(freq_hz, dish_m):
    """beam FWHM ~ 70 lambda / D (degrees)"""
    return 70.0 * (C / freq_hz) / dish_m

def t_a_from_flux(s_jy, aeff_m2):
    """antenna temperature (K) of a source of flux density s_jy: T = S A / 2k"""
    return s_jy * 1e-26 * aeff_m2 / (2 * K_B)

def sun_ta(aeff_m2, freq_hz=11.2e9, tb_k=1.0e4, radius_deg=0.28):
    """quiet-Sun antenna temperature: S = 2 k Tb Omega / lambda^2 (Sun much smaller than the beam)"""
    lam = C / freq_hz; om = math.pi * math.radians(radius_deg) ** 2
    s = 2 * K_B * tb_k * om / lam ** 2
    return s / 1e-22, t_a_from_flux(s / 1e-26, aeff_m2)

def sigma_pred(tsys_k, dv_kms, dump_s, sum_r, sum_r2, n_off):
    """rms (K) of the beam-weighted ON-OFF estimate in a dv_kms channel (radiometer equation)"""
    bw = F_HI * dv_kms * 1e3 / C
    s1 = tsys_k / math.sqrt(bw * dump_s)
    return math.sqrt(s1 ** 2 / sum_r2 + (s1 / math.sqrt(n_off) * sum_r / sum_r2) ** 2)

def drift_effective_hours(fwhm, dec_deg):
    """hours of full-response ON time a fixed beam gives a source passing through its centre"""
    omega = 15.0 * math.cos(math.radians(dec_deg))          # deg of sky per hour
    return math.sqrt(math.pi / (4 * math.log(2) * 2)) * fwhm / omega

def m31_expected(aeff_m2, tsys_k, hours_on_eff, hours_off, dv=10.0):
    bw = F_HI * dv * 1e3 / C
    sig = tsys_k * math.sqrt(1 / (bw * hours_on_eff * 3600) + 1 / (bw * hours_off * 3600))
    i_k = t_a_from_flux(M31_SINT_JY_KMS, aeff_m2)              # K km/s
    return dict(peak_mk=1e3 * 1.3 * i_k / M31_W_KMS, int_k_kms=i_k, sigma_mk=1e3 * sig,
                line_snr=i_k / (sig * math.sqrt(M31_W_KMS * dv)))

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--aeff", type=float, required=True); ap.add_argument("--tsys", type=float, default=100)
    ap.add_argument("--hours-on", type=float, default=1.0); ap.add_argument("--hours-off", type=float, default=4.0); ap.add_argument("--dish", type=float, default=None)
    a = ap.parse_args(); r = m31_expected(a.aeff, a.tsys, a.hours_on, a.hours_off)
    print(f"M31 peak antenna temperature : {r['peak_mk']:.1f} mK   (integrated {r['int_k_kms']:.2f} K km/s)")
    print(f"noise per 10 km/s channel   : {r['sigma_mk']:.1f} mK   ({a.hours_on} h full-response ON, {a.hours_off} h OFF, Tsys {a.tsys} K)")
    print(f"expected integrated line S/N: {r['line_snr']:.1f}   (need >= 5 to claim a detection)")
    if a.dish: f = fwhm_deg(F_HI, a.dish); print(f"beam FWHM {f:.1f} deg -> {drift_effective_hours(f, 41.27):.2f} h of full-response time per transit; OFF must be >= {1.5*f:.0f} deg away")
