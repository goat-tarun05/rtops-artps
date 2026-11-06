"""Telescope Control System: executes observation blocks on the (virtual) telescope.
Primitive operations - slew, track, scan, acquire - are combined into block types: calibration, pointing run, solar drift scan,
satellite tracking, H I ON/OFF tracking, H I stare. Weather triggers stow-and-resume; every block is archived with full metadata."""
import math, os, json, numpy as np
from scipy.optimize import curve_fit
from rtpipe.common import ang_sep
from rtpipe.satellites import satrec_from_elements
from .sky import Ephemeris, Track, FixedTrack, shell_elements
from .pointing import PointingModel
from .mount import Mount
from .receiver import Receiver, Scene, beam_resp
from .monitor import Monitor, hhmm_to_unix
from .archive import Archive
from .scheduler import Dispatcher, build_obs, OB

def wrap180(x): return (np.asarray(x, float) + 180.0) % 360.0 - 180.0
def gauss(t, B, A, t0, w): return B + A * np.exp(-4 * math.log(2) * (t - t0) ** 2 / w ** 2)

class TCS:
    def __init__(self, cfg, outdir):
        self.cfg, self.out = cfg, outdir; os.makedirs(outdir, exist_ok=True); s = cfg["session"]; self.date = s["date"]
        self.t0 = hhmm_to_unix(self.date, s["start"]); self.t_end = hhmm_to_unix(self.date, s["end"]); self.dt = 1.0 / cfg["mount"]["rate_hz"]
        self.rng = np.random.default_rng(s["seed"]); self.mon = Monitor(cfg); self.arc = Archive(os.path.join(outdir, "archive")); self.eph = Ephemeris(cfg)
        self.model, self.truth = PointingModel(), PointingModel(cfg["pointing_truth"]); st = cfg["mount"]["stow"]; self.mount = Mount(cfg, self.rng, az0=st["az"], el0=st["el"])
        self.rx = Receiver(cfg, np.random.default_rng(s["seed"] + 1), self.t0); self.now = self.t0; self.last_cal = {}; self.cal_ref = {}
        self.results = dict(pointing_runs=[], tracking=[], ob_log=[], stow_events=[]); W0, W1 = self.t0 - 1800, self.t_end + 1800; tg = cfg["targets"]
        self.sun = self.eph.sun(W0, W1); self.m31 = self.eph.radec("M31", tg["m31"]["ra_deg"], tg["m31"]["dec_deg"], W0, W1)
        self.m31_off = self.eph.radec("M31-OFF", tg["m31_off"]["ra_deg"], tg["m31_off"]["dec_deg"], W0, W1); self.geo = {}
        ep = (self.t0 / 86400.0 + 2440587.5); self.true_elems = shell_elements(cfg, ep); rng = np.random.default_rng(s["seed"] + 2); self.pub_elems = []
        for e in self.true_elems:
            d = dict(e); dma = rng.normal(0, cfg["satellites"]["catalogue_err_km"]) / (6378.135 + cfg["satellites"]["alt_km"]) * 180 / math.pi; d["ma"] = (e["ma"] + dma) % 360; self.pub_elems.append(d)
        self._pass_cache = {}; self.az_track = None; self.F_hi = self.rx.F_hi
        self.obs = build_obs(cfg, self); self.disp = Dispatcher(cfg, self.obs, self)
    # ------------------------------------------------------------------ context for the dispatcher
    def m31_transit(self): t = np.arange(self.t0, self.t_end + 7200, 60.0); return float(t[np.argmax(self.m31.at(t)[1])])
    def window_ok(self, ob, t0, t1):
        if not ob.type.startswith("hi_"): return True
        t = np.arange(t0, t1 + 1, 300.0); az, el = self.m31.at(t); azs, els = self.sun.at(t); c = self.cfg["constraints"]
        return bool(el.min() >= c.get("hi_min_el", 30.0) and el.max() <= c["max_el"] and ang_sep(az, el, azs, els).min() >= c["sun_avoid_deg"])
    def select_pass(self, ob, t_from, t_to):
        """pick a satellite pass satisfying el / az-rate / keyhole / weather / Sun constraints (uses the PUBLISHED catalogue)"""
        p = ob.params; lo, hi = p["max_el_range"]
        if ob.id not in self._pass_cache: self._pass_cache[ob.id] = self.eph.find_passes(self.pub_elems, t_from, t_to + 600, 20.0, p["min_el"])      # searched once per block, for ITS window
        good = []
        for q in self._pass_cache[ob.id]:
            if q["t_aos"] < t_from or q["t_los"] > t_to or not (lo <= q["max_el"] <= hi) or q["max_az_rate"] > p["max_az_rate"] or q["t_los"] - q["t_aos"] < 240 or q["az_span"] > (self.cfg["mount"]["az_range"][1] - self.cfg["mount"]["az_range"][0]) - 20: continue
            tm = np.array([q["t_aos"], 0.5 * (q["t_aos"] + q["t_los"]), q["t_los"]])
            if self.mon.wind(tm).max() > self.cfg["weather"]["stow_limit_ms"] - 1.5: continue
            good.append(q)
        if not good: return False
        q = min(good, key=lambda q: (abs(q["max_el"] - 60.0), q["t_aos"])); ob.info["pass"] = q; ob.fixed_start = q["t_aos"] - 150.0; ob.duration = q["t_los"] - q["t_aos"] + 300.0; ob.min_duration = ob.duration; return True
    # ------------------------------------------------------------------ motion primitives
    def _mount_az(self, a_mod, t=None):
        """unwrapped mount azimuth for the commanded angles. For a tracked pass (self.az_track) the branch is chosen so the WHOLE pass stays inside the
        cable-wrap range; otherwise the branch nearest the current position is used"""
        a = np.degrees(np.unwrap(np.radians(np.atleast_1d(a_mod)))); lo, hi = self.cfg["mount"]["az_range"]
        if self.az_track is not None and t is not None:
            tr = self.az_track; pmin, pmax = float(tr.azu.min()), float(tr.azu.max()); ks = [k for k in range(-3, 4) if pmin + 360 * k >= lo and pmax + 360 * k <= hi]
            if ks:
                p0 = float(np.interp(t[0], tr.t, tr.azu)); k = min(ks, key=lambda k: abs(p0 + 360 * k - self.mount.az)); return a + 360.0 * round((p0 + 360 * k - a[0]) / 360.0)
        a = a + 360.0 * round((self.mount.az - a[0]) / 360.0)
        if a.min() < lo: a = a + 360.0
        if a.max() > hi: a = a - 360.0
        return a
    def _drive(self, t, az_sky, el_sky):
        enc_az_mod, enc_el = self.model.sky_to_enc(az_sky, el_sky); caz = self._mount_az(enc_az_mod, t); n = len(caz)
        ffa = np.gradient(caz, self.dt) if n > 1 else np.zeros(n); ffe = np.gradient(enc_el, self.dt) if n > 1 else np.zeros(n); wind = self.mon.wind(t)
        az_e, el_e, d_az, d_el, az_t, el_t = self.mount.run(caz, enc_el, ffa, ffe, wind); da, de = self.truth.delta(az_t % 360, el_t)
        az_b = (az_t + da + d_az / np.cos(np.radians(el_t))) % 360.0; el_b = el_t + de + d_el
        ferr = np.hypot((az_e - caz) * np.cos(np.radians(el_e)), el_e - enc_el)
        self.mon.housekeeping(t.mean(), wind.mean(), float(np.std(az_e - caz)), float(np.std(el_e - enc_el)))
        return dict(t=t, az_sky=az_sky, el_sky=el_sky, az_cmd=caz, el_cmd=enc_el, az_enc=az_e, el_enc=el_e, az_b=az_b, el_b=el_b, ferr=ferr, wind=wind)
    def slew(self, target_fn, label, max_s=400.0):
        m = self.cfg["mount"]; tol = m["settle_tol_deg"]; n = int(2.0 / self.dt); t_start = self.now
        while self.now - t_start < max_s:
            t = self.now + self.dt * np.arange(1, n + 1); az_s, el_s = target_fn(t); d = self._drive(t, az_s, el_s); self.now = float(t[-1])
            if (d["ferr"] < tol).all(): break
        self.mon.log(self.now, "INFO", "TCS", f"slew {label} done in {self.now - t_start:.0f} s")
        return self.now - t_start
    def stow_and_wait(self, reason):
        st = self.cfg["mount"]["stow"]; w = self.cfg["weather"]; t_a = self.now; self.mon.log(self.now, "WARN", "SAFETY", f"STOW: {reason} (wind {float(self.mon.wind(self.now)):.1f} m/s)")
        self.slew(lambda t: (np.full(len(t), st["az"]), np.full(len(t), st["el"])), "to stow"); calm = 0.0
        while calm < w["resume_hold_min"] * 60 or self.mon.wind(self.now) > w["resume_limit_ms"]:
            calm = calm + 60.0 if self.mon.wind(self.now) <= w["resume_limit_ms"] else 0.0; self.now += 60.0
            if self.now > self.t_end: break
        self.mon.log(self.now, "INFO", "SAFETY", f"wind below {w['resume_limit_ms']} m/s for {w['resume_hold_min']} min: resuming"); self.results["stow_events"].append(dict(t0=t_a, t1=self.now, reason=reason)); return self.now - t_a
    def unsafe(self): return float(self.mon.wind(self.now)) > self.cfg["weather"]["stow_limit_ms"]
    def idle_until(self, t):
        if t > self.now: self.mount.vaz = self.mount.vel = 0.0; self.now = float(t)
    # ------------------------------------------------------------------ data acquisition
    def _pack(self, chunks):
        return {k: np.concatenate([c[k] for c in chunks]) for k in chunks[0]} if chunks else {}
    def acquire_ku(self, duration, cmd_fn, scene, cadence, qa=None, wind_abort=True, chunk=60.0):
        k = int(round(cadence / self.dt)); t_stop = self.now + duration; pts, cnt = [], []; aborted = False
        while self.now < t_stop - 1e-6:
            if wind_abort and self.unsafe(): aborted = True; break
            nc = min(int(chunk / self.dt), int(round((t_stop - self.now) / self.dt))) // k * k
            if nc < k: break
            t = self.now + self.dt * np.arange(1, nc + 1); az_s, el_s = cmd_fn(t); d = self._drive(t, az_s, el_s); tc, c, elm = self.rx.ku_counts(t, d["az_b"], d["el_b"], scene, cadence, self.dt)
            if qa is not None: azq, elq = qa.at(t); d["qa_err"] = ang_sep(d["az_b"], d["el_b"], azq, elq)
            pts.append(d); cnt.append(dict(t=tc, counts=c)); self.now = float(t[-1])
        return dict(pt=self._pack(pts), t=np.concatenate([c["t"] for c in cnt]) if cnt else np.array([]), counts=np.concatenate([c["counts"] for c in cnt]) if cnt else np.array([]), aborted=aborted)
    def acquire_hi(self, n_dumps, dump_s, cmd_fn, is_on_fn=None):
        tcs, resp, theta, pts, flags = [], [], [], [], []; nd = int(round(dump_s / self.dt)); abort = False
        for i in range(n_dumps):
            if self.unsafe(): abort = True; break
            t = self.now + self.dt * np.arange(1, nd + 1); az_s, el_s = cmd_fn(t); d = self._drive(t, az_s, el_s); azm, elm = self.m31.at(t)
            resp.append(float(beam_resp(d["az_b"], d["el_b"], azm, elm, self.F_hi).mean())); theta.append(float(ang_sep(d["az_b"], d["el_b"], azm, elm).mean())); tcs.append(float(t.mean())); pts.append(d); self.now = float(t[-1])
            flags.append(True if is_on_fn is None else is_on_fn(self.now))
        return dict(t=np.array(tcs), resp=np.array(resp), theta=np.array(theta), is_on=np.array(flags, bool), pt=self._pack(pts), aborted=abort)
    def _pointing_tables(self, d, step=10):
        p = d["pt"]; s = slice(None, None, step); tab = dict(t=p["t"][s], az_sky_cmd=p["az_sky"][s], el_sky_cmd=p["el_sky"][s], az_enc=p["az_enc"][s], el_enc=p["el_enc"][s], following_err=p["ferr"][s], wind=p["wind"][s])
        return {"pointing": tab, "truth": dict(t=p["t"][s], az_bore=p["az_b"][s], el_bore=p["el_b"][s])}
    def _hdr(self, ob, t0, t1, **kw):
        h = dict(TSTART=float(t0), TSTOP=float(t1), TELESCOP=self.cfg["telescope"]["name"], SITELAT=self.cfg["site"]["lat_deg"], SITELONG=self.cfg["site"]["lon_deg"], SITEELEV=self.cfg["site"]["height_m"], DISH_M=self.cfg["telescope"]["dish_m"],
                 PTGMODEL=json.dumps({k: round(float(v), 4) for k, v in self.model.p.items()}), BASEID=ob.base_id or ob.id); h.update(kw); return h
    def _qa(self, ob, d, label):
        p = d["pt"]
        if "qa_err" in p and len(p["qa_err"]): self.results["tracking"].append(dict(ob=ob.id, label=label, rms_true_deg=float(np.sqrt(np.mean(p["qa_err"] ** 2))), p95_true_deg=float(np.percentile(p["qa_err"], 95)), rms_following_deg=float(np.sqrt(np.mean(p["ferr"] ** 2))), n=int(len(p["qa_err"]))))
    # ------------------------------------------------------------------ block types
    def run_cal(self, kind):
        t_a = self.now; obid = f"CAL_{kind.upper()}_{int((t_a - self.t0) // 60):04d}"; m = self.cfg["mount"]
        class _O: id = obid; base_id = obid
        ob = _O()
        if kind == "ku":
            cs = m["cold_sky"]; self.slew(lambda t: (np.full(len(t), cs["az"]), np.full(len(t), cs["el"])), "cold sky"); d = self.acquire_ku(30.0, lambda t: (np.full(len(t), cs["az"]), np.full(len(t), cs["el"])), Scene(), 1.0, wind_abort=False)
            th = self.now + np.arange(1, 31.0); hot = self.rx.ku_load(th, 1.0, self.cfg["ku"]["t_hot_k"]); self.now += 30.0
            tabs = dict(cold=dict(t=d["t"], counts=d["counts"]), hot=dict(t=th, counts=hot)); tabs.update(self._pointing_tables(d))
        else:
            if self.rx.freq is None: self.rx.hi_tune(self.now)
            th = self.now + 10.0 * np.arange(1, 11); hot = self.rx.hi_hot(th, 10.0); self.now += 100.0
            tabs = dict(hot=dict(t=th, spec=hot), freq=dict(freq_hz=self.rx.freq))
        self.last_cal[kind] = self.now; self.cal_ref[kind] = obid; self.arc.write(obid, f"cal_{kind}", self._hdr(ob, t_a, self.now), tabs, target="calibration"); self.mon.log(self.now, "INFO", "CAL", f"{kind} calibration {obid} stored")
        return obid
    def _scan(self, trk, scene, axis, L, rate, ident):
        def cmd_start(t):
            az, el = trk.at(t); return (az - L / np.cos(np.radians(el)), el) if axis == "az" else (az, el - L)
        self.slew(cmd_start, f"{ident} {axis}-scan start"); n = int(2 * L / rate / self.dt); t = self.now + self.dt * np.arange(1, n + 1); off = -L + rate * (t - t[0] + self.dt)
        az_t, el_t = trk.at(t); az_s = az_t + (off / np.cos(np.radians(el_t)) if axis == "az" else 0.0); el_s = el_t + (off if axis == "el" else 0.0)
        d = self._drive(t, az_s, el_s); self.now = float(t[-1]); tc, c, _ = self.rx.ku_counts(t, d["az_b"], d["el_b"], scene, 0.1, self.dt)
        tr = tc - tc[0]                                                              # fit in scan-relative time (Unix seconds would ruin the finite-difference steps)
        try:
            p, cov = curve_fit(gauss, tr, c, p0=[c.min(), c.max() - c.min(), tr[np.argmax(c)], 5.0], bounds=([0, 0, tr[0], 1.0], [np.inf, np.inf, tr[-1], 60.0])); B, A, t0, w = p; t0 = t0 + tc[0]
        except Exception: return None, dict(t=tc, counts=c)
        snr = A / max(1e-9, np.std(c[: max(10, n // 8)] - B)); azs, els = trk.at(t0); ae, ee = float(np.interp(t0, t, d["az_enc"])), float(np.interp(t0, t, d["el_enc"]))
        da, de = self.model.delta(ae, ee)
        res = float(wrap180(azs - ae) - da) if axis == "az" else float((els - ee) - de)
        return dict(target=ident, axis=axis, az_enc=ae, el_enc=ee, resid=res, snr=float(snr), width_s=float(w), t=float(t0)), dict(t=tc, counts=c)
    def run_pointing_run(self, ob):
        t_a = self.now; p = ob.params; scans = []; meas = []; amp_geo = self.cfg["receiver"]["geo_beacon_k"]; F = self.cfg["ku"]["fwhm_theory_deg"]; L = 2.4 * F; rate = 0.25
        geos = [(lon, self.geo.setdefault(lon, self.eph.geo(lon))) for lon in sorted(p["geo"], key=lambda l: self.eph.geo(l).az)]; items = [("geo", g) for g in geos]; ns = p.get("sun", 0)
        for i in range(ns): items.insert(int(round(i * len(items) / max(1, ns - 1 if ns > 1 else 1))) if ns > 1 else 0, ("sun", None))
        for kind, g in items:
            if kind == "sun": trk, scene, name = self.sun, Scene().add("sun", self.sun, self.rx.T_sun), "Sun"
            else:
                trk, name = g[1], g[1].name; az0, el0 = trk.at(np.array([self.now])); azs, els = self.sun.at(np.array([self.now]))
                if ang_sep(az0, el0, azs, els)[0] < 12.0 and els[0] > 0: self.mon.log(self.now, "INFO", "PTG", f"skip {name}: {ang_sep(az0, el0, azs, els)[0]:.1f} deg from the Sun"); continue
                if el0[0] < self.cfg["constraints"]["min_el"]: continue
                scene = Scene().add("geo", trk, amp_geo)
            for axis in ("az", "el"):
                m, raw = self._scan(trk, scene, axis, L, rate, name); scans.append(dict(target=name, axis=axis, **raw))
                if m is None or m["snr"] < 20: self.mon.log(self.now, "WARN", "PTG", f"{name} {axis}-scan: no usable peak"); continue
                meas.append(m)
        a = [m for m in meas if m["axis"] == "az"]; e = [m for m in meas if m["axis"] == "el"]; base = PointingModel(self.model.vector())
        rms_before = float(np.sqrt(np.mean([m["resid"] ** 2 for m in meas]))) if meas else float("nan"); new, err, rms_after = PointingModel.fit([m["az_enc"] for m in a], [m["el_enc"] for m in a], [m["resid"] for m in a], [m["az_enc"] for m in e], [m["el_enc"] for m in e], [m["resid"] for m in e], base)
        g_az, g_el = np.meshgrid(np.arange(60, 301, 20.0), np.arange(20, 71, 10.0)); g_az, g_el = g_az.ravel(), g_el.ravel()
        def sky_rms(model):
            da, de = (self.truth.delta(g_az, g_el)[0] - model.delta(g_az, g_el)[0]) * np.cos(np.radians(g_el)), self.truth.delta(g_az, g_el)[1] - model.delta(g_az, g_el)[1]; return float(np.sqrt(np.mean(da ** 2 + de ** 2)))
        res = dict(ob=ob.id, t0=t_a, t1=self.now, n_meas=len(meas), rms_resid_before_deg=rms_before, rms_resid_after_fit_deg=rms_after, params_before=base.vector().tolist(), params_fit=new.vector().tolist(), params_err=err.tolist(),
                   map_rms_true_before_deg=sky_rms(base), map_rms_true_after_deg=sky_rms(new), applied=bool(p.get("update_model", False)), meas=meas)
        if p.get("update_model", False): self.model = new
        self.mon.log(self.now, "INFO", "PTG", f"{ob.id}: residual rms {rms_before:.3f} deg -> {rms_after:.3f} deg; model {'APPLIED' if res['applied'] else 'not applied (verification)'}"); self.results["pointing_runs"].append(res)
        flat = {k: np.concatenate([s[k] for s in scans]) for k in ("t", "counts")}; sid = np.concatenate([np.full(len(s["t"]), i) for i, s in enumerate(scans)])
        tabs = dict(scans=dict(t=flat["t"], counts=flat["counts"], scan_id=sid), measurements=dict(target=np.array([m["target"] for m in meas]), axis=np.array([m["axis"] for m in meas]), az_enc=np.array([m["az_enc"] for m in meas]), el_enc=np.array([m["el_enc"] for m in meas]), resid=np.array([m["resid"] for m in meas]), snr=np.array([m["snr"] for m in meas])))
        self.arc.write(ob.id, "pointing_run", self._hdr(ob, t_a, self.now, NMEAS=len(meas), RMSB=rms_before, RMSA=rms_after), tabs, target="Sun+GEO beacons"); return "DONE", {}
    def run_solar_drift(self, ob, start):
        tc = ob.info.get("t_center", start + ob.duration / 2); az, el = self.sun.at(np.array([tc])); off = ob.params["el_offset_deg"]; az0, el0 = float(az[0]), float(el[0]) + off
        self.slew(lambda t: (np.full(len(t), az0), np.full(len(t), el0)), f"{ob.id} stare position"); self.idle_until(max(self.now, start)); t_a = self.now
        d = self.acquire_ku(ob.duration, lambda t: (np.full(len(t), az0), np.full(len(t), el0)), Scene().add("sun", self.sun, self.rx.T_sun), 1.0)
        tabs = dict(data=dict(t=d["t"], counts=d["counts"])); tabs.update(self._pointing_tables(d)); st = "ABORTED_WIND" if d["aborted"] else "DONE"; self._qa(ob, d, "solar stare")
        self.arc.write(ob.id, "solar_drift", self._hdr(ob, t_a, self.now, AZ_CMD=az0, EL_CMD=el0, EL_OFFS=off, FREQ_GHZ=self.cfg["ku"]["freq_ghz"]), tabs, target="Sun", status=st, cal_ref=self.cal_ref.get("ku", "")); return st, {}
    def run_sat_track(self, ob, start):
        q = ob.info["pass"]; name = q["name"]; tp = [e for e in self.pub_elems if e["name"] == name][0]; tt = [e for e in self.true_elems if e["name"] == name][0]
        t0, t1 = q["t_aos"] - 90.0, q["t_los"] + 90.0; pub = self.eph.satellite(name, satrec_from_elements(tp), t0, t1, 1.0); tru = self.eph.satellite(name, satrec_from_elements(tt), t0, t1, 1.0); sat = self.cfg["satellites"]
        self.az_track = pub                                                          # keep the whole pass inside the cable-wrap range (it may cross north)
        self.slew(lambda t: pub.at(np.maximum(t, q["t_aos"])), f"{name} acquisition"); self.idle_until(max(self.now, q["t_aos"] - 5)); t_a = self.now
        try: d = self.acquire_ku(q["t_los"] - self.now, lambda t: pub.at(t), Scene().add("sat", tru, sat["t_sat_k"], sat["ref_range_km"]), 0.5, qa=tru, chunk=30.0)
        finally: self.az_track = None
        self._qa(ob, d, f"track {name}")
        tg = np.arange(q["t_aos"], q["t_los"], 1.0); az_p, el_p = pub.at(tg); tabs = dict(data=dict(t=d["t"], counts=d["counts"]), ephem=dict(t=tg, az=az_p, el=el_p, range_km=pub.range_at(tg))); tabs.update(self._pointing_tables(d))
        st = "ABORTED_WIND" if d["aborted"] else "DONE"; self.arc.write(ob.id, "sat_track", self._hdr(ob, t_a, self.now, SATNAME=name, MAX_EL=q["max_el"]), tabs, target=name, status=st, cal_ref=self.cal_ref.get("ku", ""))
        ob.info["sat"] = dict(q); return st, {}
    def _hi_common(self, ob, d, t_a, mode, extra, tabs_extra=None):
        spec = self.rx.hi_spectra(d["t"], d["resp"], d["theta"], self._dump_s) if len(d["t"]) else np.zeros((0, len(self.rx.freq)), np.float32)
        tabs = dict(data=dict(t=d["t"], spec=spec, resp=d["resp"], theta=d["theta"], is_on=d["is_on"]), freq=dict(freq_hz=self.rx.freq)); tabs.update(self._pointing_tables(d) if len(d["t"]) else {})
        st = "ABORTED_WIND" if d["aborted"] else "DONE"; seg = ob.attempts; oid = ob.id if seg == 0 else f"{ob.base_id}_r{seg}"
        self.arc.write(oid, f"hi_{mode}", self._hdr(ob, t_a, self.now, MODE=mode, FC_HZ=float(self.rx.lo_log[-1][1]), DUMP_S=self._dump_s, **extra), tabs, target="M31", status=st, cal_ref=self.cal_ref.get("hi", "")); return st
    def run_hi_track_nod(self, ob, start, duration):
        p = ob.params; self._dump_s = p["dump_s"]; cyc = p["cycle_min"] * 60.0
        if self.rx.freq is None: self.rx.hi_tune(self.now)
        t_end = max(self.now, start) + duration
        self.slew(lambda t: self.m31.at(t), f"{ob.id} M31 ON"); self.idle_until(max(self.now, start)); t_a = self.now; parts = []; k = 0; abort = False
        while self.now < t_end - self._dump_s and not abort:
            on = (k % 2 == 0); trk = self.m31 if on else self.m31_off
            self.slew(lambda t: trk.at(t), f"{ob.id} {'ON' if on else 'OFF'}"); n = int(min(cyc, t_end - self.now) // self._dump_s)
            if n < 1: break
            d = self.acquire_hi(n, self._dump_s, lambda t, trk=trk: trk.at(t), lambda tt, on=on: on); parts.append(d); abort = d["aborted"]
            for a_ in (d,): 
                if len(a_["t"]): qa = ang_sep(a_["pt"]["az_b"], a_["pt"]["el_b"], *(trk.at(a_["pt"]["t"]))); a_["pt"]["qa_err"] = qa
            k += 1
        d = {key: (np.concatenate([x[key] for x in parts]) if key != "pt" else self._pack([x["pt"] for x in parts if len(x["t"])])) for key in ("t", "resp", "theta", "is_on", "pt")}; d["aborted"] = abort
        st = self._hi_common(ob, d, t_a, "track_nod", dict(ON_ID="M31", OFF_RA=self.cfg["targets"]["m31_off"]["ra_deg"], OFF_DEC=self.cfg["targets"]["m31_off"]["dec_deg"])); self._qa(ob, d, "M31 nod")
        got = float(len(d["t"]) * self._dump_s); return st, dict(got_s=got, remaining_s=max(0.0, duration - (self.now - t_a)))
    def run_hi_stare(self, ob, start, duration):
        self._dump_s = ob.params["dump_s"]; tc = ob.info.get("t_center", start + duration / 2); az, el = self.m31.at(np.array([tc])); az0, el0 = float(az[0]), float(el[0])
        if self.rx.freq is None: self.rx.hi_tune(self.now)
        self.slew(lambda t: (np.full(len(t), az0), np.full(len(t), el0)), f"{ob.id} stare position"); self.idle_until(max(self.now, start)); t_a = self.now
        d = self.acquire_hi(int(duration // self._dump_s), self._dump_s, lambda t: (np.full(len(t), az0), np.full(len(t), el0)))
        st = self._hi_common(ob, d, t_a, "stare", dict(AZ_CMD=az0, EL_CMD=el0)); return st, dict(got_s=float(len(d["t"]) * self._dump_s), remaining_s=0.0)
    # ------------------------------------------------------------------ the session loop
    def execute(self, ob, start, dur):
        t_a = self.now; self.mon.log(self.now, "INFO", "TCS", f"BEGIN {ob.id} ({ob.type}) planned {start - self.t0:.0f}..{start + dur - self.t0:.0f} s")
        if ob.type == "pointing_run": st, info = self.run_pointing_run(ob)
        elif ob.type == "solar_drift": st, info = self.run_solar_drift(ob, start)
        elif ob.type == "sat_track": st, info = self.run_sat_track(ob, start)
        elif ob.type == "hi_track_nod": st, info = self.run_hi_track_nod(ob, start, dur)
        elif ob.type == "hi_stare": st, info = self.run_hi_stare(ob, start, dur)
        else: raise ValueError(ob.type)
        self.results["ob_log"].append(dict(id=ob.id, type=ob.type, t0=t_a, t1=self.now, status=st, **{k: v for k, v in info.items()})); self.mon.log(self.now, "INFO", "TCS", f"END {ob.id}: {st}")
        return st, info
    def run_session(self):
        self.mon.log(self.now, "INFO", "TCS", f"session start; {len(self.obs)} blocks queued")
        self.results["planned"] = self.disp.dry_run(self.t0, self.t_end)
        while True:
            act = self.disp.next(self.now, self.last_cal, self.t_end)
            if act[0] == "done": break
            if act[0] == "idle": self.idle_until(act[1]); continue
            _, ob, start, dur, cal = act
            if self.unsafe(): self.stow_and_wait("wind above limit before block start"); continue
            if cal: self.run_cal(cal)
            st, info = self.execute(ob, start, dur)
            if st == "ABORTED_WIND":
                self.stow_and_wait(f"wind during {ob.id}"); rem = info.get("remaining_s", 0.0)
                if ob.type == "hi_track_nod" and rem >= ob.min_duration * 0.5: ob.attempts += 1; ob.duration = rem; ob.min_duration = min(ob.min_duration, rem * 0.5); self.mon.log(self.now, "INFO", "SCHED", f"{ob.id} re-queued with {rem / 60:.0f} min remaining")
                else: ob.status = "DONE"
            else: ob.status = "DONE"
        self.mon.log(self.now, "INFO", "TCS", "session complete"); self.mon.save(os.path.join(self.out, "events.jsonl")); return self.results
