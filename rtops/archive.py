"""Data archive: one FITS file per observation block (primary header + binary tables) and an SQLite catalogue."""
import os, json, sqlite3, numpy as np
from astropy.io import fits
from datetime import datetime, timezone

def _table(name, cols):
    c = []
    for k, v in cols.items():
        v = np.asarray(v)
        if v.ndim == 1: c.append(fits.Column(name=k.upper(), format={"f": "D", "i": "K", "b": "L", "U": "16A", "S": "16A"}.get(v.dtype.kind, "D") if v.dtype.kind not in "US" else "16A", array=v))
        else: c.append(fits.Column(name=k.upper(), format=f"{v.shape[1]}E", dim=f"({v.shape[1]})", array=v.astype(np.float32)))
    h = fits.BinTableHDU.from_columns(c); h.name = name.upper(); return h

class Archive:
    def __init__(self, root):
        self.root = root; os.makedirs(root, exist_ok=True); self.db = sqlite3.connect(os.path.join(root, "catalogue.sqlite"))
        self.db.execute("CREATE TABLE IF NOT EXISTS obs(obid TEXT PRIMARY KEY, obtype TEXT, target TEXT, origin TEXT, t_start REAL, t_stop REAL, path TEXT, status TEXT, qc TEXT, cal_ref TEXT, notes TEXT)"); self.db.commit()
    def write(self, obid, obtype, header, tables, origin="SIM", target="", status="DONE", qc="OK", cal_ref="", notes=""):
        day = datetime.fromtimestamp(header["TSTART"], timezone.utc).strftime("%Y-%m-%d"); d = os.path.join(self.root, day); os.makedirs(d, exist_ok=True); path = os.path.join(d, f"{obid}.fits")
        ph = fits.PrimaryHDU(); ph.header["OBID"], ph.header["OBTYPE"], ph.header["ORIGIN"], ph.header["TARGET"] = obid, obtype, origin, target
        for k, v in header.items():
            ph.header[k[:8] if len(k) > 8 else k] = v if not isinstance(v, (dict, list)) else json.dumps(v)[:60]
        ph.header["COMMENT"] = "rtops archive; SIM = produced by the digital twin, REAL = ingested recording"
        fits.HDUList([ph] + [_table(n, c) for n, c in tables.items()]).writeto(path, overwrite=True)
        self.db.execute("INSERT OR REPLACE INTO obs VALUES (?,?,?,?,?,?,?,?,?,?,?)", (obid, obtype, target, origin, header["TSTART"], header["TSTOP"], path, status, qc, cal_ref, notes)); self.db.commit()
        return path
    def entries(self, **where):
        q = "SELECT obid,obtype,target,origin,t_start,t_stop,path,status,qc,cal_ref,notes FROM obs"; vals = []
        if where: q += " WHERE " + " AND ".join(f"{k}=?" for k in where); vals = list(where.values())
        cols = ["obid", "obtype", "target", "origin", "t_start", "t_stop", "path", "status", "qc", "cal_ref", "notes"]
        return [dict(zip(cols, r)) for r in self.db.execute(q + " ORDER BY t_start", vals)]
    @staticmethod
    def read(path):
        with fits.open(path) as h:
            hdr = dict(h[0].header); tabs = {x.name.lower(): {c.name.lower(): np.array(x.data[c.name]) for c in x.columns} for x in h[1:]}
        return hdr, tabs
