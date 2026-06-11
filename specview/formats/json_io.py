"""Reader/writer for spectra stored as JSON.

Handles the common shapes:
  * {"x": [...], "y": [...]}                       (any aliased key names)
  * {"name":..,"x_unit":"nm", "wavelength":[...], "absorbance":[...]}
  * {"data": [[x, y], ...]}                         (list of pairs)
  * {"spectra": [ {..}, {..} ]}                     (multiple spectra)
  * [ [x, y], ... ]   or   [ {..}, {..} ]           (top-level list)
"""
from __future__ import annotations

import json
import os

import numpy as np

from ..spectrum import Spectrum
from . import _hints as H


def _as_array(v):
    """Return a 1-D float array if v is a numeric list, else None."""
    if isinstance(v, (list, tuple)) and v and isinstance(v[0], (int, float)):
        return np.asarray(v, dtype=float)
    return None


def _spectrum_from_obj(obj, path, default_name):
    if not isinstance(obj, dict):
        return None
    keys = list(obj.keys())
    name = obj.get(H.match_key(keys, H.NAME_KEYS) or "", default_name) or default_name
    x_unit = H.norm_xunit(obj.get(H.match_key(keys, H.XUNIT_KEYS) or "", None))
    y_unit = H.norm_yunit(obj.get(H.match_key(keys, H.YUNIT_KEYS) or "", None))
    laser = obj.get(H.match_key(keys, H.LASER_KEYS) or "", None)

    xk = H.match_key(keys, H.X_KEYS)
    yk = H.match_key(keys, H.Y_KEYS)
    x = _as_array(obj.get(xk)) if xk else None
    y = _as_array(obj.get(yk)) if yk else None
    if x is None or y is None:
        # fall back to a "data" array of [x, y] pairs
        d = obj.get("data") or obj.get("points")
        if isinstance(d, list) and d and isinstance(d[0], (list, tuple)):
            arr = np.asarray(d, dtype=float)
            if arr.ndim == 2 and arr.shape[1] >= 2:
                x, y = arr[:, 0], arr[:, 1]
    if x is None or y is None:
        return None
    meta = {"source": path}
    if laser:
        try:
            meta["laser_nm"] = float(laser)
        except (TypeError, ValueError):
            pass
    return Spectrum(x=x, y=y, name=str(name), x_unit=x_unit or "pixel",
                    y_unit=y_unit or "intensity", meta=meta)


def load_json(path: str) -> list[Spectrum]:
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    base = os.path.splitext(os.path.basename(path))[0]
    spectra: list[Spectrum] = []

    # {"spectra": [...]}
    if isinstance(data, dict) and isinstance(data.get("spectra"), list):
        for i, item in enumerate(data["spectra"]):
            s = _spectrum_from_obj(item, path, f"{base} [{i + 1}]")
            if s:
                spectra.append(s)
    # top-level list
    elif isinstance(data, list) and data:
        if isinstance(data[0], (list, tuple)):           # list of [x, y] pairs
            arr = np.asarray(data, dtype=float)
            if arr.ndim == 2 and arr.shape[1] >= 2:
                spectra.append(Spectrum(arr[:, 0], arr[:, 1], name=base,
                                        meta={"source": path}))
        elif isinstance(data[0], dict):                  # list of spectrum objects
            for i, item in enumerate(data):
                s = _spectrum_from_obj(item, path, f"{base} [{i + 1}]")
                if s:
                    spectra.append(s)
    # single object
    elif isinstance(data, dict):
        s = _spectrum_from_obj(data, path, base)
        if s:
            spectra.append(s)

    if not spectra:
        raise ValueError(f"Could not find x/y data in {os.path.basename(path)}.")
    return spectra


def save_json(spec: Spectrum, path: str) -> None:
    obj = {"name": spec.name, "x_unit": spec.x_unit, "y_unit": spec.y_unit,
           "x": spec.x.tolist(), "y": spec.y.tolist()}
    if spec.meta.get("laser_nm"):
        obj["laser_nm"] = spec.meta["laser_nm"]
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=2)
