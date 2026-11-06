"""Observation blocks, constraints and the dispatcher (priority + gap filling + fixed-time blocks + dynamic re-planning)."""
import math, numpy as np
from dataclasses import dataclass, field
from rtpipe.common import ang_sep
from .monitor import hhmm_to_unix

@dataclass
class OB:
    id: str; type: str; priority: int = 5; earliest: float = 0.0; deadline: float = 1e12; fixed_start: float = None
    duration: float = 0.0; min_duration: float = 0.0; params: dict = field(default_factory=dict); status: str = "PENDING"; attempts: int = 0
    info: dict = field(default_factory=dict); base_id: str = None; est_overhead: float = 0.0

def build_obs(cfg, ctx):
    date = cfg["session"]["date"]; out = []
    for o in cfg["observations"]:
        ob = OB(id=o["id"], type=o["type"], priority=o.get("priority", 5), params=dict(o.get("params", {})), base_id=o["id"])
        ob.earliest = hhmm_to_unix(date, o["earliest"]) if "earliest" in o else 0.0; ob.deadline = hhmm_to_unix(date, o["deadline"]) if "deadline" in o else 1e12
        ob.duration = o.get("duration_min", 0) * 60.0; ob.min_duration = o.get("min_duration_min", o.get("duration_min", 0)) * 60.0
        if "fixed_center" in o:
            c = o["fixed_center"]; tc = ctx.m31_transit() if c == "transit" else hhmm_to_unix(date, c); ob.fixed_start = tc - ob.duration / 2; ob.info["t_center"] = tc
        if ob.type == "pointing_run": n = len(ob.params.get("geo", [])) + ob.params.get("sun", 0); ob.duration = ob.min_duration = n * 110.0 + 60.0
        out.append(ob)
    return out

class Dispatcher:
    def __init__(self, cfg, obs, ctx):
        self.cfg, self.obs, self.ctx = cfg, obs, ctx; c = cfg["constraints"]; self.margin = c["slew_margin_s"]; self.cal_interval = cfg["session"]["cal_interval_min"] * 60.0
        self.cal_ku, self.cal_hi = 200.0, 170.0                                    # estimated durations of the calibration blocks
    def pending(self): return [o for o in self.obs if o.status == "PENDING"]
    def cal_needed(self, ob, now, last_cal):
        kind = "hi" if ob.type.startswith("hi_") else ("ku" if ob.type in ("solar_drift", "sat_track") else None)
        if kind is None: return None
        t = last_cal.get(kind); return kind if (t is None or now - t > self.cal_interval) else None
    def next(self, now, last_cal, t_end):
        """-> ('run', ob, start, duration, cal_kind) | ('idle', until) | ('done',)"""
        while True:
            pend = self.pending()
            if not pend or now >= t_end - 60: return ("done",)
            for ob in pend:                                                         # target-of-opportunity blocks: choose the pass when their window opens
                if ob.type == "sat_track" and ob.fixed_start is None and now >= ob.earliest - 900:
                    if not self.ctx.select_pass(ob, max(now + self.margin + 300, ob.earliest), ob.deadline): ob.status = "SKIPPED"; ob.info["reason"] = "no pass satisfies the constraints"
            pend = self.pending(); fixed = sorted([o for o in pend if o.fixed_start is not None], key=lambda o: o.fixed_start); nf = fixed[0] if fixed else None
            lead = self.margin + max(self.cal_ku, self.cal_hi)
            if nf is not None and nf.fixed_start - now <= lead:
                if now > nf.fixed_start + 600: nf.status = "SKIPPED"; nf.info["reason"] = "start window missed"; continue
                cal = self.cal_needed(nf, now, last_cal); return ("run", nf, max(now, nf.fixed_start), nf.duration, cal)
            cands = [o for o in pend if o.fixed_start is None and o.type != "sat_track" and o.earliest - 60 <= now and now + o.min_duration <= o.deadline]
            for ob in sorted(cands, key=lambda o: (-o.priority, o.deadline)):
                cal = self.cal_needed(ob, now, last_cal); ovh = self.margin + (0 if cal is None else (self.cal_hi if cal == "hi" else self.cal_ku))
                avail = min(ob.deadline, (nf.fixed_start - self.margin) if nf is not None else 1e12, t_end) - now - ovh
                dur = min(ob.duration, avail)
                if dur >= ob.min_duration and self.ctx.window_ok(ob, now + ovh, now + ovh + dur): return ("run", ob, now + ovh, dur, cal)
            ev = [o.earliest for o in pend if o.earliest > now and o.fixed_start is None]
            if nf is not None: ev.append(nf.fixed_start - lead)
            ev += [now + 900.0]; return ("idle", max(now + 60.0, min(ev)))
    def dry_run(self, t0, t_end):
        """planned timeline using estimated durations (no physics, no weather)"""
        saved = [(o.status, dict(o.info), o.fixed_start) for o in self.obs]; now, last, plan = t0, {}, []
        while True:
            a = self.next(now, last, t_end)
            if a[0] == "done": break
            if a[0] == "idle": now = a[1]; continue
            _, ob, start, dur, cal = a
            if cal: plan.append(dict(id=f"CAL_{cal}", type="cal", t0=start - (self.cal_hi if cal == "hi" else self.cal_ku), t1=start)); last[cal] = start
            plan.append(dict(id=ob.id, type=ob.type, t0=start, t1=start + dur)); ob.status = "DONE"; now = start + dur
        for o, (s, i, f) in zip(self.obs, saved): o.status, o.info, o.fixed_start = s, i, f
        return plan
