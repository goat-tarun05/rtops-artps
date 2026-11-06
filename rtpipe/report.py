"""Self-contained HTML report: results, PASS/FAIL checks, figures, injection-recovery comparison, assumptions and limits."""
import base64, json, os, datetime

def _img(path):
    with open(path, "rb") as f: return "data:image/png;base64," + base64.b64encode(f.read()).decode()

def _checks(c):
    return "".join(f"<li class={'ok' if v else 'bad'}>{'PASS' if v else 'FAIL'} — {k}</li>" for k, v in c.items())

def _f(x, nd=2): return "—" if x is None else f"{x:.{nd}f}"

def _mc(out):
    p = f"{out}/monte_carlo_hi.json"
    if not os.path.exists(p): return ""
    m = json.load(open(p))
    return (f"<p><b>Robustness ({m['n']} simulated noise realisations):</b> v<sub>sys</sub> = {m['vsys_mean']:.1f} ± {m['vsys_std']:.1f} km/s, "
            f"width {m['width_mean']:.0f} ± {m['width_std']:.0f} km/s, integrated-intensity ratio {m['int_ratio_mean']:.2f} ± {m['int_ratio_std']:.2f}, "
            f"mean integrated S/N {m['snr_mean']:.1f}. Error calibration: pull mean {m['pull_mean']:+.2f}, pull spread {m['pull_std']:.2f} (1.0 = quoted errors are correct).</p>")

def run(res, out="outputs", truth=None, cfg=None):
    s, st, h, tg = res["solar"], res["starlink"], res["hi"], res["tag"]; sim = truth is not None
    rec = ""
    if sim:
        ts, th, tt = truth["solar_params"], truth["hi_params"], truth["hi_truth"]
        rows = [("Beam FWHM (deg)", ts["fwhm_deg"], s["fwhm_deg"]), ("Pointing offset, azimuth on sky (deg)", ts["daz_sky_deg"], s["pointing_dx_sky_deg"]),
                ("Pointing offset, elevation (deg)", ts["del_deg"], s["pointing_dy_deg"]), ("Sun peak above sky (K)", ts["t_sun_k"], s["sun_peak_k"]),
                ("Receiver temperature (K)", ts["t_rx_k"], s["t_rx_k"]), ("M31 systemic velocity, heliocentric (km/s)", th["vsys_helio_kms"], h["vsys_helio_kms"]),
                ("M31 integrated intensity (K km/s)", tt["int_k_kms"], h["integrated_k_kms"]),
                ("Starlink passes in the beam / identified", truth["n_sat_events"], st["n_identified"]), ("Terrestrial RFI bursts / flagged unidentified", truth["n_rogue"], st["n_unidentified"])]
        rec = ("<h2>Injection–recovery test (simulated data)</h2><table><tr><th>Quantity</th><th>Injected (truth)</th><th>Recovered</th></tr>" +
               "".join(f"<tr><td>{a}</td><td>{b:.2f}</td><td>{c:.2f}</td></tr>" if isinstance(b, float) else f"<tr><td>{a}</td><td>{b}</td><td>{c}</td></tr>" for a, b, c in rows) + "</table>")
    det = "".join(f"<tr><td>{d['id']}</td><td>{d['peak_k']:.0f}</td><td>{d['dur_s']:.1f}</td><td>{_f(d['dt_s'], 1)}</td><td>{_f(d['min_sep_deg'])}</td></tr>" for d in st["detections"])
    hi_c = cfg["hi"] if cfg else {}
    limits = [
        "The demo data are <b>simulated</b>. Injected: Sun 200 K (quiet-Sun flux × A<sub>eff</sub>); M31 from its literature HI flux × the assumed antenna area; Starlink burst brightness is an <b>assumption</b>.",
        "Injection–recovery tests the <b>code</b>, not the physics: simulator and pipeline share conventions (Gaussian beam 70λ/D, SGP4, the Doppler sign, the disc template). The Doppler sign was checked independently by hand.",
        "Pointing: the quoted uncertainty is statistical. Setting the dish by hand adds a systematic of about 0.3°.",
        "M31's receding edge (−71 km/s) is blended with Milky Way hydrogen. That range is excluded and the width needs a prior (245 ± 25 km/s here, 16 km/s above the injected value), so v<sub>sys</sub> carries about +5 km/s of prior-driven systematic.",
        "The Sun's finite size (0.56°) widens the fitted beam by about 0.7 % (ignored). The beam is treated as Gaussian.",
        "Starlink: the demo constellation is a synthetic single shell, truth and analysis use the same SGP4 code, and the Ku-band mask applies only to Ku data (it is not relevant at 1420 MHz).",
        "The ideal S/N prediction assumes white noise and a known baseline; the fitted S/N is lower because baseline and blending are fitted too."]
    html = f"""<!DOCTYPE html><html lang=en><head><meta charset=utf-8><meta name=viewport content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>Radio Telescope Automation &amp; Tracking — pipeline report</title><style>
:root{{box-sizing:border-box;padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px);--bg:#f6f4ef;--ink:#1b1e21;--mute:#5d646b;--line:#d9d4c8;--card:#fff}}
@media (prefers-color-scheme:dark){{:root:not([data-theme="light"]){{--bg:#15181b;--ink:#e8e6e1;--mute:#9aa1a8;--line:#30363c;--card:#1d2125}}}}
body{{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}}
main{{max-width:900px;margin:auto;padding:18px}} h1{{font-size:22px}} h2{{font-size:17px;margin-top:26px}}
.card{{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px 14px;margin:10px 0;overflow-x:auto}}
table{{border-collapse:collapse;width:100%;font-size:13.5px}} td,th{{border-bottom:1px solid var(--line);padding:5px 7px;text-align:left}}
img{{max-width:100%;border-radius:8px;background:#fff}} .ok{{color:#2e7d4f}} .bad{{color:#b3261e}} ul{{padding-left:18px}} .mute{{color:var(--mute)}}
</style></head><body><main>
<h1>Radio Telescope Automation &amp; Tracking — pipeline report</h1>
<p class=mute>Generated {datetime.datetime.utcnow().isoformat(timespec='seconds')} UTC by <code>python run_all.py</code>. Data: <b>{'SIMULATED injection datasets' if sim else 'observed data'}</b>. HI antenna: {hi_c.get('dish_m','?')} m, A<sub>eff</sub> {hi_c.get('aeff_m2','?')} m², mode <b>{h['mode']}</b>.</p>
{rec}
<h2>1 · Solar cross-scan — celestial tracking &amp; calibration (Ku band)</h2><div class=card>
<img src="{_img(f'{out}/solar_transit.png')}" alt="solar transit">
<table><tr><td>Y-factor / receiver temp / system temp</td><td>{s['y_factor']:.2f} / {s['t_rx_k']:.0f} K / {s['t_sys_k']:.0f} K</td></tr>
<tr><td>Sun above sky (peak)</td><td>{s['sun_peak_k']:.0f} K</td></tr>
<tr><td>Beam FWHM (theory 70λ/D)</td><td>{s['fwhm_deg']:.2f} ± {s['fwhm_err_deg']:.2f}° ({s['fwhm_theory_deg']}°)</td></tr>
<tr><td>Pointing offset (az on sky, el)</td><td>{s['pointing_dx_sky_deg']:+.2f}°, {s['pointing_dy_deg']:+.2f}° ± {s['pointing_err_deg']:.2f}° (statistical) — manual setting adds ≈ {s['pointing_systematic_floor_deg']}°</td></tr>
<tr><td>Transit predicted → observed</td><td>{', '.join(f'{a[11:19]} → {b[11:19]} ({o:+.0f} s)' for a, b, o in zip(s['transit_pred_utc'], s['transit_obs_utc'], s['transit_offset_s']))} UTC</td></tr>
<tr><td>Starlink samples masked in these Ku scans</td><td>{s['masked_samples']}</td></tr>
<tr><td>Beam centre vs Sun at the peaks (pointing model applied)</td><td>{', '.join(f"{tg[f'scan{k}_beam_to_sun_deg']:.2f}° (commanded {tg[f'scan{k}_commanded_offset_deg']:.1f}°)" for k in (1, 2))}</td></tr></table>
<ul>{_checks(s['checks'])}</ul></div>
<h2>2 · Starlink — orbital tracking, identification, Ku RFI mask</h2><div class=card>
<img src="{_img(f'{out}/starlink.png')}" alt="starlink">
<p>{st['n_predicted']} passes predicted from orbital elements · {st['n_detected']} bursts detected · {st['n_identified']} identified · {st['n_unidentified']} unidentified (flagged as terrestrial RFI) · timing RMS {_f(st['timing_rms_s'])} s · {100*st['masked_fraction']:.1f} % of samples masked</p>
<table><tr><th>Identified as</th><th>Peak (K)</th><th>Duration (s)</th><th>Δt vs prediction (s)</th><th>Closest approach (°)</th></tr>{det}</table>
<ul>{_checks(st['checks'])}</ul></div>
<h2>3 · Andromeda (M31) H I — {h['mode']} mode, Doppler tracking</h2><div class=card>
<img src="{_img(f'{out}/hi_m31.png')}" alt="M31 HI">
<table><tr><td>Beam FWHM · ON / OFF spectra · effective full-response time</td><td>{h['beam_fwhm_deg']:.1f}° · {h['n_on']} / {h['n_off']} · {h['effective_on_hours']:.2f} h</td></tr>
<tr><td>Systemic velocity (heliocentric / LSR)</td><td>{h['vsys_helio_kms']:.1f} ± {h['vsys_err_kms']:.1f} / {h['vsys_lsr_kms']:.1f} km/s (literature ≈ −300 helio)</td></tr>
<tr><td>Rotation speed V sin i · width 2V sin i</td><td>{h['vrot_sini_kms']:.0f} ± {h['vrot_err_kms']:.0f} · {h['width_kms']:.0f} km/s</td></tr>
<tr><td>Integrated intensity: measured / expected from literature flux × A<sub>eff</sub></td><td>{h['integrated_k_kms']:.1f} ± {h['integrated_err_k_kms']:.1f} / {h['integrated_expected_k_kms']:.1f} K km/s (ratio {h['integrated_ratio']:.2f}; tests T<sub>sys</sub> scale and aperture efficiency)</td></tr>
<tr><td>Line S/N: fitted · ideal white-noise prediction</td><td>{h['line_snr']:.1f} · {h['ideal_line_snr']:.1f}</td></tr>
<tr><td>3σ upper limit on the integrated intensity (use it if there is no detection)</td><td>{h['int_ul3_k_kms']:.1f} K km/s = {h['int_ul3_jy_kms']:.0f} Jy km/s (M31's expected flux ≈ 31 400 Jy km/s){' — <b>weaker than the expected flux: no constraint</b>' if h['int_ul3_jy_kms'] > 31400 else ''}</td></tr>
<tr><td>Noise per 10 km/s channel: measured / radiometer prediction</td><td>{1e3*h['rms_k']:.1f} / {1e3*h['rms_predicted_k']:.1f} mK (ratio {h['rms_ratio']:.2f}) at T<sub>sys</sub> {h['tsys_k']:.0f} K</td></tr>
<tr><td>Doppler correction applied per spectrum</td><td>{h['doppler_corr_range_kms'][0]:.2f} … {h['doppler_corr_range_kms'][1]:.2f} km/s</td></tr>
<tr><td>RFI channels flagged · Milky Way residual (LSR)</td><td>{h['rfi_channels_flagged']} · {h['milky_way_residual_lsr_kms']:.1f} km/s</td></tr></table>
<ul>{_checks(h['checks'])}</ul>{_mc(out)}</div>
<h2>Assumptions and limits</h2><div class=card><ul>{''.join(f'<li>{x}</li>' for x in limits)}</ul></div>
<h2>Outputs</h2><div class=card><ul>{''.join(f'<li><code>tables/{t}</code></li>' for t in tg['tables'])}<li><code>results.json</code>, <code>rfi_mask.json</code>, <code>starlink_detections.csv</code>, <code>m31_spectrum.ecsv</code>, <code>plan.csv</code> (from <code>run_all.py plan</code>)</li></ul></div>
</main></body></html>"""
    with open(f"{out}/report.html", "w") as f: f.write(html)
    return f"{out}/report.html"
