"""Self-contained HTML report of a digital-observatory session."""
import base64, io, json, datetime, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from .pointing import PointingModel, PARAMS

COL = dict(pointing_run="#7b6cd9", solar_drift="#e5a100", sat_track="#d1495b", hi_track_nod="#2a9d8f", hi_stare="#4c78a8", cal="#9aa0a6")
def b64(fig): b = io.BytesIO(); fig.savefig(b, format="png", dpi=120, bbox_inches="tight"); plt.close(fig); return "data:image/png;base64," + base64.b64encode(b.getvalue()).decode()
def png(path):
    with open(path, "rb") as f: return "data:image/png;base64," + base64.b64encode(f.read()).decode()
def hm(t): return datetime.datetime.utcfromtimestamp(t).strftime("%H:%M")
def tbl(head, rows, cls=""): return f"<div class=tw><table class='{cls}'><tr>" + "".join(f"<th>{h}</th>" for h in head) + "</tr>" + "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows) + "</table></div>"

def fig_timeline(res, hk, cfg, arc, t0):
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 4.8), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1]))
    def bar(y, e, c, lab=None): a1.barh(y, max((e["t1"] - e["t0"]) / 3600, 0.02), left=(e["t0"] - t0) / 3600, color=c, edgecolor="white", height=0.7)
    for p in res["planned"]: bar(2, p, COL.get(p["type"], "#999"))
    for o in res["ob_log"]: bar(1, o, COL.get(o["type"], "#999"))
    for e in arc.entries():
        if e["obtype"].startswith("cal_"): bar(1, dict(t0=e["t_start"], t1=e["t_stop"]), COL["cal"])
    for s in res["stow_events"]: a1.barh(1, (s["t1"] - s["t0"]) / 3600, left=(s["t0"] - t0) / 3600, color="none", edgecolor="#c0392b", hatch="////", height=0.7)
    for o in res["ob_log"]:
        a1.text((o["t0"] - t0) / 3600 + 0.05, 1.45, o["id"], fontsize=7, rotation=45, va="bottom")
    a1.set_yticks([1, 2]); a1.set_yticklabels(["executed", "planned"]); a1.set_ylim(0.4, 2.9)
    a1.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c) for c in COL.values()] + [plt.Rectangle((0, 0), 1, 1, fill=False, hatch="////", edgecolor="#c0392b")], labels=[k.replace("_", " ") for k in COL] + ["weather stow"], ncol=7, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, 1.25), frameon=False)
    h = np.array(hk); a2.plot((h[:, 0] - t0) / 3600, h[:, 1], lw=0.8, color="#34495e"); a2.axhline(cfg["weather"]["stow_limit_ms"], color="#c0392b", ls="--", lw=1); a2.set_ylabel("wind (m/s)"); a2.set_xlabel(f"hours since {hm(t0)} UTC")
    return b64(fig)

def sky_err(params, truth, az, el):
    da = (truth.delta(az, el)[0] - PointingModel(params).delta(az, el)[0]) * np.cos(np.radians(el)); de = truth.delta(az, el)[1] - PointingModel(params).delta(az, el)[1]; return np.hypot(da, de)

def fig_pointing(res, truth):
    az, el = np.meshgrid(np.linspace(60, 300, 60), np.linspace(20, 70, 30)); prs = res["pointing_runs"]; sets = [("no pointing model", np.zeros(5))]
    if prs: sets.append((f"after {prs[0]['ob']} (fitted model)", np.array(prs[0]["params_fit"])))
    fig, axs = plt.subplots(1, len(sets), figsize=(5.4 * len(sets), 3.6), sharey=True); axs = np.atleast_1d(axs)
    for ax, (lab, p) in zip(axs, sets):
        e = sky_err(p, truth, az, el); im = ax.pcolormesh(az, el, e, cmap="viridis", vmin=0, vmax=0.5, shading="auto"); ax.set_title(f"{lab}\nrms {np.sqrt(np.mean(e ** 2)):.3f}° (true error over the sky)", fontsize=9); ax.set_xlabel("azimuth (°)")
        for pr in prs:
            for m in pr["meas"]: ax.plot(m["az_enc"], m["el_enc"], ".", color="white", ms=3)
    axs[0].set_ylabel("elevation (°)"); fig.colorbar(im, ax=axs, label="pointing error (°)", shrink=0.9); return b64(fig)

def fig_sat(q):
    t = np.array(q["t"]); fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 4.6), sharex=True, gridspec_kw=dict(height_ratios=[2, 1])); tm = (t - t[0]) / 60
    a1.plot(tm, q["ideal"], "k--", lw=1, label="ideal (perfect pointing, ∝ 1/range²)"); a1.plot(tm, q["T"], lw=0.8, color="#d1495b", label="measured, tracked from the published orbit"); a1.set_ylabel("satellite signal (K)"); a1.legend(fontsize=8)
    a2.plot(tm, q["el"], color="#4c78a8"); a2.set_ylabel("elevation (°)"); a2.set_xlabel("minutes after acquisition"); return b64(fig)

def fig_hi(R, comb):
    fig, ax = plt.subplots(figsize=(6.4, 2.8)); names = list(R["hi"]) + ["combined"]; v = [R["hi"][k]["vsys_helio_kms"] for k in R["hi"]] + [comb["v"]]; e = [R["hi"][k]["vsys_err_kms"] for k in R["hi"]] + [comb["e"]]
    ax.errorbar(v, range(len(v)), xerr=e, fmt="o", color="#2a9d8f", capsize=4); ax.axvline(-300, color="k", ls="--", lw=1, label="injected / literature −300"); ax.set_yticks(range(len(v))); ax.set_yticklabels(names); ax.set_xlabel("M31 systemic velocity (km/s)"); ax.legend(fontsize=8); return b64(fig)

SVG = """<svg viewBox="0 0 860 330" role="img" aria-label="architecture" style="width:100%;max-width:860px"><defs><marker id="a" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto"><path d="M2 1L8 5L2 9" fill="none" stroke="currentColor" stroke-width="1.5"/></marker></defs>
<g font-family="system-ui,sans-serif" font-size="12" fill="currentColor" stroke="currentColor"><g fill="none" stroke-width="1.2">
<rect x="10" y="20" width="150" height="54" rx="8"/><rect x="200" y="20" width="170" height="54" rx="8"/><rect x="410" y="20" width="170" height="54" rx="8"/><rect x="620" y="20" width="220" height="54" rx="8"/>
<rect x="200" y="130" width="170" height="54" rx="8"/><rect x="410" y="130" width="170" height="54" rx="8"/><rect x="620" y="130" width="220" height="54" rx="8"/>
<rect x="10" y="240" width="150" height="54" rx="8"/><rect x="200" y="240" width="170" height="54" rx="8"/><rect x="410" y="240" width="170" height="54" rx="8"/><rect x="620" y="240" width="220" height="54" rx="8"/>
<path d="M160 47H198" marker-end="url(#a)"/><path d="M370 47H408" marker-end="url(#a)"/><path d="M580 47H618" marker-end="url(#a)"/><path d="M285 74V128" marker-end="url(#a)"/><path d="M495 74V128" marker-end="url(#a)"/>
<path d="M370 157H408" marker-end="url(#a)"/><path d="M580 157H618" marker-end="url(#a)"/><path d="M730 184V238" marker-end="url(#a)"/><path d="M620 267H582" marker-end="url(#a)"/><path d="M410 267H372" marker-end="url(#a)"/><path d="M200 267H162" marker-end="url(#a)"/></g>
<g stroke="none"><text x="85" y="42" text-anchor="middle" font-weight="600">Observation queue</text><text x="85" y="60" text-anchor="middle" opacity=".7">blocks + constraints</text>
<text x="285" y="42" text-anchor="middle" font-weight="600">Scheduler</text><text x="285" y="60" text-anchor="middle" opacity=".7">priority, windows, replan</text>
<text x="495" y="42" text-anchor="middle" font-weight="600">Telescope control (TCS)</text><text x="495" y="60" text-anchor="middle" opacity=".7">slew / track / scan</text>
<text x="730" y="42" text-anchor="middle" font-weight="600">Monitor + weather</text><text x="730" y="60" text-anchor="middle" opacity=".7">stow, resume, alarms</text>
<text x="285" y="152" text-anchor="middle" font-weight="600">Ephemeris</text><text x="285" y="170" text-anchor="middle" opacity=".7">Sun, M31, GEO, satellites</text>
<text x="495" y="152" text-anchor="middle" font-weight="600">Pointing model</text><text x="495" y="170" text-anchor="middle" opacity=".7">fitted from pointing runs</text>
<text x="730" y="152" text-anchor="middle" font-weight="600">Servo / mount (virtual)</text><text x="730" y="170" text-anchor="middle" opacity=".7">rate limits, wind, encoders</text>
<text x="85" y="262" text-anchor="middle" font-weight="600">Report</text><text x="85" y="280" text-anchor="middle" opacity=".7">HTML + figures</text>
<text x="285" y="262" text-anchor="middle" font-weight="600">Reduction (rtpipe)</text><text x="285" y="280" text-anchor="middle" opacity=".7">solar, satellite, H I</text>
<text x="495" y="262" text-anchor="middle" font-weight="600">FITS archive</text><text x="495" y="280" text-anchor="middle" opacity=".7">+ SQLite catalogue</text>
<text x="730" y="262" text-anchor="middle" font-weight="600">Receivers (virtual)</text><text x="730" y="280" text-anchor="middle" opacity=".7">Ku power, L-band spectra, cal</text></g></g></svg>"""

def build(out, cfg, res, planned_arc, hk, arc, R, ab, ing, wall, truth, tcs_wall=None):
    t0 = res["t0"]; prs = res["pointing_runs"]; css = """:root{box-sizing:border-box;padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px);--bg:#f6f4ef;--ink:#1b1e21;--mute:#5d646b;--line:#d9d4c8;--card:#fff}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#15181b;--ink:#e8e6e1;--mute:#9aa1a8;--line:#30363c;--card:#1d2125}}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}main{max-width:980px;margin:auto;padding:18px}h1{font-size:23px}h2{font-size:18px;margin-top:30px;border-bottom:1px solid var(--line);padding-bottom:4px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px 14px;margin:10px 0}.kpi{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:10px}.kpi div{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:8px 12px}.kpi b{font-size:21px;display:block}
.tw{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:13px}td,th{border-bottom:1px solid var(--line);padding:5px 7px;text-align:left}img{max-width:100%;border-radius:8px;background:#fff}.mute{color:var(--mute)}.ok{color:#2e7d4f}.bad{color:#b3261e}code{font-size:12.5px}ul{padding-left:18px}"""
    trk = [t for t in res["tracking"]]; sat = R.get("sat"); sol = R.get("solar"); hi = R.get("hi", {})
    # combined HI
    vs = [(h["vsys_helio_kms"], h["vsys_err_kms"], h["line_snr"]) for h in hi.values()]
    if vs: w = np.array([1 / e ** 2 for _, e, _ in vs]); comb = dict(v=float(np.sum(w * np.array([v for v, _, _ in vs])) / w.sum()), e=float(1 / np.sqrt(w.sum())), snr=float(np.sqrt(sum(s ** 2 for _, _, s in vs))))
    pr1 = prs[0] if prs else None
    kpi = [("Blocks executed", f"{len([o for o in res['ob_log']])} + {len([e for e in arc.entries() if e['obtype'].startswith('cal_')])} cal"), ("Pointing error", f"{pr1['map_rms_true_before_deg']:.2f}° → {pr1['map_rms_true_after_deg']:.3f}°" if pr1 else "—"),
           ("Servo following", f"{np.mean([t['rms_following_deg'] for t in trk]) * 1000:.1f} m° rms" if trk else "—"), ("Weather stows", f"{len(res['stow_events'])} ({sum(s['t1'] - s['t0'] for s in res['stow_events']) / 60:.0f} min)"),
           ("Satellite on-target", f"{sat['time_above_half_s']:.0f} s vs {sat['fixed_beam_equiv_s']:.1f} s fixed" if sat else "—"), ("M31 v_sys (combined)", f"{comb['v']:.0f} ± {comb['e']:.0f} km/s" if vs else "—")]
    h = [f"<!DOCTYPE html><html lang=en><head><meta charset=utf-8><meta name=viewport content='width=device-width, initial-scale=1, viewport-fit=cover'><title>rtops — digital observatory session report</title><style>{css}</style></head><body><main>",
         f"<h1>Digital observatory — session report</h1><p class=mute>{cfg['telescope']['name']} · {cfg['telescope']['dish_m']} m dish · {cfg['session']['date']} {cfg['session']['start']}–{cfg['session']['end']} UTC · generated {datetime.datetime.utcnow().isoformat(timespec='seconds')}Z · <b>all data SIMULATED by a software telescope; no hardware involved</b> · full session simulated in {wall:.0f} s.</p>",
         "<div class=kpi>" + "".join(f"<div><span class=mute>{k}</span><b>{v}</b></div>" for k, v in kpi) + "</div>",
         "<h2>1 · Architecture: what a real control system does, in software</h2><div class=card>" + SVG + "</div>",
         tbl(["Layer here", "Does", "Real-observatory counterpart"], [["Scheduler / dispatcher", "blocks with priorities, windows, fixed times; fills gaps; re-queues after interruptions; picks satellite passes by constraints", "observing queue / dynamic scheduling"], ["Ephemeris", "apparent az/el with refraction for sidereal, Sun, GEO, satellites", "ephemeris + coordinate chain in the antenna controller"],
          ["Pointing model", "TPOINT-style terms (IA, IE, NPAE, AN, AW) fitted from automated pointing runs, inverted for every command", "pointing model from pointing runs"], ["Servo / mount", "rate and acceleration limits, feed-forward tracking, encoder quantisation, wind jitter, cable wrap", "drive servo + axis encoders"],
          ["Calibration blocks", "hot-load / cold-sky Y-factor before science, refreshed hourly", "noise diode / load calibration"], ["Safety + monitor", "stow above the wind limit, resume after a calm period, event log, housekeeping", "wind/safety interlocks, operator console"],
          ["Archive", "one FITS file per block + SQLite catalogue with provenance", "FITS/SDFITS archive + metadata database"], ["Reduction", "calibration, beam and pointing fits, line fitting from the archive only", "science pipeline"]]),
         "<h2>2 · Automation: the day as scheduled and as executed</h2><div class=card><img src='" + fig_timeline(res, hk, cfg, arc, t0) + "' alt=timeline></div>"]
    ev = [e for e in res["events"] if e["level"] == "WARN" or "re-queued" in e["msg"] or "skip" in e["msg"] or "resuming" in e["msg"] or "APPLIED" in e["msg"] or "BEGIN" in e["msg"] or "END" in e["msg"]]
    h.append("<div class=card><b>Operator log (selected)</b>" + tbl(["UTC", "level", "source", "message"], [[e["utc"], e["level"], e["src"], e["msg"]] for e in ev[:40]]) + "</div>")
    h.append("<h2>3 · Automatic pointing calibration</h2>")
    if pr1:
        rows = [[k, f"{t:+.3f}", f"{v:+.3f} ± {er:.3f}", f"{v - t:+.3f}"] for k, t, v, er in zip(PARAMS, cfg["pointing_truth"], pr1["params_fit"], pr1["params_err"])]
        h.append("<p>The telescope starts with an unknown pointing error (hidden truth in the simulator). The first block scans the Sun and eight geostationary beacons in azimuth and elevation, fits the peak of each scan, and fits the model. Every later command uses it.</p>")
        h.append(tbl(["term (deg)", "hidden truth", "fitted ± error", "difference"], rows) + f"<p class=mute>{pr1['n_meas']} scan peaks · residual rms {pr1['rms_resid_before_deg']:.3f}° before → {pr1['rms_resid_after_fit_deg']:.4f}° after. Terms can be partly degenerate over a limited part of the sky; what matters is the error over the whole sky:</p><div class=card><img src='" + fig_pointing(res, truth) + "' alt='pointing maps'></div>")
        if len(prs) > 1: h.append(f"<p>Verification run {prs[1]['ob']} at the end of the day, with the model applied and <i>not</i> refitted: remaining residual rms <b>{prs[1]['rms_resid_before_deg']:.3f}°</b> (true sky error {prs[1]['map_rms_true_before_deg']:.3f}°).</p>")
    if ab and ab.get("solar") and sol:
        rows = [["Sun height at the beam centre (K)", f"{ab['solar']['sun_peak_k']:.0f}", f"{sol['sun_peak_k']:.0f}"], ["pointing offset seen in the Sun cross-scan (°)", f"({ab['solar']['pointing_dx_sky_deg']:+.2f}, {ab['solar']['pointing_dy_deg']:+.2f})", f"({sol['pointing_dx_sky_deg']:+.2f}, {sol['pointing_dy_deg']:+.2f})"]]
        if ab.get("sat") and sat: rows.append(["satellite tracking efficiency", f"{ab['sat']['mean_eff']:.2f}", f"{sat['mean_eff']:.2f}"]); rows.append(["peak satellite signal (K)", f"{ab['sat']['peak_k']:.0f}", f"{sat['peak_k']:.0f}"])
        rows.append(["true pointing error while tracking (°, rms)", f"{np.mean(ab['tracking_rms']):.3f}", f"{np.mean([t['rms_true_deg'] for t in trk]):.3f}"])
        h.append("<h3>Why it matters: the same day with and without the pointing calibration</h3>" + tbl(["quantity", "without calibration", "with calibration"], rows) + "<p class=mute>The 9.8° L-band beam hides a 0.4° error; the 1.25° Ku beam does not.</p>")
    h.append("<h2>4 · Tracking performance</h2>")
    h.append(tbl(["block", "target", "true pointing error rms (°)", "95th percentile (°)", "servo following error rms (°)", "samples"], [[t["ob"], t["label"], f"{t['rms_true_deg']:.3f}", f"{t['p95_true_deg']:.3f}", f"{t['rms_following_deg']:.4f}", f"{t['n']}"] for t in trk]) + "<p class=mute>'True' error compares the real boresight (known only inside the simulator) with the target. 'Servo following' is what a real operator can see: encoder minus command.</p>")
    if sat:
        h.append(f"<div class=card><b>Satellite pass — {sat['name']}</b> (max elevation {sat['max_el']:.0f}°, selected automatically by the scheduler)<img src='{fig_sat(sat)}' alt='satellite pass'></div>")
        h.append(tbl(["quantity", "tracked from the published orbit", "fixed dish pointing at the pass centre"], [["time within half power (s)", f"{sat['time_above_half_s']:.0f}", f"{sat['fixed_beam_equiv_s']:.1f}"], ["tracking efficiency (mean measured / ideal)", f"{sat['mean_eff']:.2f}", "—"], ["peak signal (K) · S/N", f"{sat['peak_k']:.0f} · {sat['snr_peak']:.0f}", "—"]]) + "<p class=mute>Efficiency below 1 is the catalogue error: the published orbit is a few km off, which points the beam slightly away. The satellite brightness is an assumption.</p>")
    h.append("<h2>5 · Science products (reduced from the archive)</h2>")
    if R.get("errors"): h.append("<div class=card><b>Products that could not be reduced in this configuration</b><ul>" + "".join(f"<li><code>{k}</code>: {v}</li>" for k, v in R["errors"].items()) + "</ul><p class=mute>Usually the configuration leaves no usable data: for example a beam so wide that the drift scan never reaches an OFF position.</p></div>")
    if sol: h.append(f"<div class=card><b>Sun cross-scan</b> ({', '.join(sol['blocks'])}): beam FWHM {sol['fwhm_deg']:.2f}° (theory {sol['fwhm_theory_deg']:.2f}° plus the Sun's size → 1.33°) · pointing offset ({sol['pointing_dx_sky_deg']:+.3f}, {sol['pointing_dy_deg']:+.3f})° · T_sys {sol['t_sys_k']:.0f} K · Sun {sol['sun_peak_k']:.0f} K<img src='{png(out + '/reduced/solar/solar_transit.png')}' alt=solar></div>")
    for k, r in hi.items(): h.append(f"<div class=card><b>M31 H I — {k}</b> ({', '.join(r['blocks'])}): v_sys {r['vsys_helio_kms']:.1f} ± {r['vsys_err_kms']:.1f} km/s · line S/N {r['line_snr']:.1f} (ideal {r['ideal_line_snr']:.1f}) · noise {r['rms_ratio']:.2f}× the radiometer prediction · {r['on_source_hours']:.2f} h on source<br><img src='{png(out + '/reduced/' + ('hitrk' if 'tracked' in k else 'hidft') + '/hi_m31.png')}' alt='hi {k}'></div>")
    if vs:
        eff = []
        for k, r in hi.items(): hrs = sum(e["t_stop"] - e["t_start"] for e in arc.entries() if e["obtype"].startswith("hi_") and any(e["obid"] == b for b in r["blocks"])) / 3600; eff.append([k, f"{hrs:.2f}", f"{r['on_source_hours']:.2f}", f"{r['ideal_line_snr']:.1f}", f"{r['ideal_line_snr'] / np.sqrt(hrs):.1f}", f"{r['line_snr']:.1f}"])
        h.append("<h3>Tracking versus a fixed dish</h3>" + tbl(["mode", "telescope time (h)", "full-response time on source (h)", "ideal S/N", "ideal S/N per √(telescope hour)", "fitted S/N"], eff) + f"<div class=card><img src='{fig_hi(R, comb)}' alt='v_sys'></div><p>Combining the two independent measurements (inverse-variance weights): v_sys = <b>{comb['v']:.1f} ± {comb['e']:.1f} km/s</b>, line S/N ≈ {comb['snr']:.1f}. The ideal S/N is the white-noise prediction; the fitted S/N is lower because the baseline, the Milky Way blend and the width prior are fitted too.</p>")
    h.append("<h2>6 · Archive catalogue</h2>" + tbl(["block", "type", "target", "origin", "start (UTC)", "minutes", "status", "calibration used", "FITS file"], [[e["obid"], e["obtype"], e["target"], e["origin"], hm(e["t_start"]), f"{(e['t_stop'] - e['t_start']) / 60:.0f}", e["status"], e["cal_ref"], e["path"].split("archive/")[-1]] for e in arc.entries()]))
    h.append("<h2>7 · Plugging in your real recordings</h2><p>The reduction reads only the archive. A real solar-transit scan or hydrogen-line recording enters through an adapter that writes the same kind of FITS block with <code>ORIGIN=REAL</code>; nothing downstream changes.</p>")
    if ing:
        h.append(tbl(["recording (simulated stand-in)", "quantity", "from the live archive", "after export → ingest as REAL → reduce", "difference"], ing) + "<p class=mute>This round trip exports blocks to the plain raw-file format, re-ingests them through the adapter, and reduces them again: identical results show that the real-data path is wired.</p>")
    h.append("<ul><li><code>rtops.ingest.ingest_solar(arc, cfg, npz, meta, obid)</code> — arrays <code>unix, counts, hot, cold</code> and the nominal <code>az_deg, el_deg, el_offset_set_deg</code>.</li><li><code>rtops.ingest.ingest_hi(arc, cfg, npz, meta, obid)</code> — arrays <code>freq_hz, unix, spec, hot</code> (+ <code>is_on</code> for ON/OFF data), and <code>mode</code>, target, az/el or the OFF position.</li></ul>")
    h.append("<h2>Assumptions and limits</h2><div class=card><ul><li>The telescope is a <b>digital twin</b>: signals are computed from where the virtual beam really points (command + hidden pointing errors + wind). Nothing here is an observation.</li><li>Hidden truth (pointing terms, wind, satellite brightness, M31 flux × effective area) is used only by the simulator and for the 'true error' columns; the control system and the reduction never read it.</li><li>The simulator and the pointing model use the same functional form, so recovery of the terms is optimistic. Real structures add flexure and other terms.</li><li>Satellite brightness and the GEO beacon level are assumptions. The catalogue error (2 km) is typical of fresh elements.</li><li>Beams are Gaussian (70λ/D); the Sun is a Gaussian-equivalent disc; refraction uses a standard formula in both the simulator and the controller.</li><li>L-band pointing is insensitive to a 0.4° error; Ku-band is not. Results scale with <code>telescope.dish_m</code> in the configuration.</li></ul></div></main></body></html>")
    open(out + "/report.html", "w").write("\n".join(h)); return out + "/report.html"
