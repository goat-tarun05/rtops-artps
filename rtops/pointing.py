"""Pointing model (TPOINT-style alt-az terms) and its least-squares fit from a pointing run."""
import numpy as np
PARAMS = ["IA", "IE", "NPAE", "AN", "AW"]      # az zero point, el zero point, az-el non-perpendicularity, az-axis tilt north / west (degrees)

class PointingModel:
    def __init__(self, values=None):
        v = [0.0] * len(PARAMS) if values is None else list(values); self.p = dict(zip(PARAMS, v))
    def vector(self): return np.array([self.p[k] for k in PARAMS])
    def delta(self, az, el):
        """boresight minus encoder, degrees: (azimuth-angle offset, elevation offset)"""
        A, E = np.radians(az), np.radians(el); tE = np.tan(E); p = self.p
        return (p["IA"] + p["NPAE"] * tE + p["AN"] * np.sin(A) * tE - p["AW"] * np.cos(A) * tE, p["IE"] + p["AN"] * np.cos(A) + p["AW"] * np.sin(A))
    def sky_to_enc(self, az, el):
        """encoder angles that put the boresight on the sky position (az, el)"""
        az, el = np.asarray(az, float), np.asarray(el, float); ea, ee = az.copy(), el.copy()
        for _ in range(4):
            da, de = self.delta(ea, ee); ea, ee = az - da, el - de
        return ea % 360.0, ee
    @staticmethod
    def design(az_a, el_a, az_e, el_e):
        """rows for the azimuth measurements (taken at az_a, el_a) and for the elevation measurements (taken at az_e, el_e)"""
        Aa, Ea = np.radians(az_a), np.radians(el_a); Ae = np.radians(az_e); na, ne = len(Aa), len(Ae); tE = np.tan(Ea)
        ra = np.column_stack([np.ones(na), np.zeros(na), tE, np.sin(Aa) * tE, -np.cos(Aa) * tE]); re = np.column_stack([np.zeros(ne), np.ones(ne), np.zeros(ne), np.cos(Ae), np.sin(Ae)])
        return np.vstack([ra, re])
    @classmethod
    def fit(cls, az_a, el_a, d_az, az_e, el_e, d_el, base=None):
        """fit the INCREMENT to `base` from residual offsets (source - encoder - base model); returns model, 1-sigma errors, rms residual"""
        X = cls.design(np.asarray(az_a), np.asarray(el_a), np.asarray(az_e), np.asarray(el_e)); y = np.r_[d_az, d_el]; coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        res = y - X @ coef; s2 = float(res @ res) / max(1, len(y) - len(coef)); cov = np.linalg.pinv(X.T @ X) * s2
        new = (np.zeros(len(PARAMS)) if base is None else base.vector()) + coef
        return cls(new), np.sqrt(np.diag(cov)), float(np.sqrt(np.mean(res ** 2)))
