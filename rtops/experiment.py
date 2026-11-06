"""Run a variant of the session with a few settings changed - the workhorse of every experiment in the assignment.
    from rtops.experiment import run_variant
    r = run_variant({"pointing_truth": [0.8, -0.5, 0.1, 0.2, -0.3]}, out="runs/big_errors")
    r["results"]["pointing_runs"][0]["map_rms_true_after_deg"]
Override keys are dotted paths; list positions are numbers, e.g. "observations.0.params.geo" or "weather.gusts.0.speed"."""
import copy, os
from . import config
from .tcs import TCS
from .reduce import reduce_all

def set_path(cfg, path, value):
    d = cfg; keys = path.split(".")
    for k in keys[:-1]: d = d[int(k)] if isinstance(d, list) else d[k]
    last = keys[-1]
    if isinstance(d, list): d[int(last)] = value
    else: d[last] = value

def run_variant(overrides=None, out="runs/variant", config_path=None, reduce=True):
    cfg = config.load(config_path)
    for k, v in (overrides or {}).items(): set_path(cfg, k, copy.deepcopy(v))
    cfg = config.derive(cfg); tcs = TCS(cfg, out); res = tcs.run_session(); res["t0"] = tcs.t0; res["events"] = tcs.mon.events
    R = reduce_all(tcs.arc, cfg, os.path.join(out, "reduced")) if reduce else None
    return dict(cfg=cfg, results=res, reduced=R, tcs=tcs)
