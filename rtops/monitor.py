"""Weather model, event log and housekeeping - the 'operator console' of the digital observatory."""
import json, math, numpy as np
from datetime import datetime, timezone

def hhmm_to_unix(date, hhmm):
    h, m = map(int, hhmm.split(":")); return datetime(*map(int, date.split("-")), h, m, tzinfo=timezone.utc).timestamp()

class Monitor:
    def __init__(self, cfg):
        self.cfg = cfg; self.date = cfg["session"]["date"]; self.events = []; self.hk = []
        w = cfg["weather"]; self.base = w["base_ms"]; self.gusts = [(hhmm_to_unix(self.date, g["start"]), hhmm_to_unix(self.date, g["end"]), g["speed"]) for g in w["gusts"]]
        self.t0 = hhmm_to_unix(self.date, cfg["session"]["start"])
    def wind(self, t):
        t = np.asarray(t, float); w = self.base + 1.5 * np.sin(2 * np.pi * (t - self.t0) / (6 * 3600.0)) + 0.4 * np.sin(2 * np.pi * t / 311.0)
        for a, b, s in self.gusts:
            ramp = np.clip((t - a) / 300.0, 0, 1) * np.clip((b - t) / 300.0, 0, 1); w = w + (s - self.base) * ramp
        return np.clip(w, 0.2, None)
    def log(self, t, level, src, msg, **kw):
        self.events.append(dict(t=float(t), utc=datetime.fromtimestamp(t, timezone.utc).strftime("%H:%M:%S"), level=level, src=src, msg=msg, **kw))
    def housekeeping(self, t, wind, az_err, el_err): self.hk.append((float(t), float(wind), float(az_err), float(el_err)))
    def save(self, path):
        with open(path, "w") as f:
            for e in self.events: f.write(json.dumps(e, default=float) + "\n")
