import math, os, yaml, copy
def derive(cfg):
    """derived quantities (effective area, beam widths) - call again after changing telescope.dish_m or the frequencies"""
    t = cfg["telescope"]; aeff = t["efficiency"] * math.pi * (t["dish_m"] / 2) ** 2; cfg["telescope"]["aeff_m2"] = aeff
    for k in ("ku", "hi"): cfg[k]["dish_m"] = t["dish_m"]; cfg[k]["aeff_m2"] = aeff           # rtpipe modules read these
    cfg["ku"]["fwhm_theory_deg"] = 70.0 * (2.998e8 / (cfg["ku"]["freq_ghz"] * 1e9)) / t["dish_m"]; return cfg
def load(path=None):
    path = path or os.environ.get("RTOPS_CONFIG", "config/session.yaml"); return derive(yaml.safe_load(open(path)))
