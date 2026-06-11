"""Spectral reference library: build/save/load, and similarity search.

A library is just a named collection of reference spectra, persisted as JSON.
Searching an unknown spectrum returns a ranked hit list with several similarity
scores. Spectra are put on a common x-grid (overlap) before comparison.
"""
from __future__ import annotations

import json
import os

import numpy as np

from .spectrum import Spectrum

_EPS = 1e-12


# ----------------------------------------------------------- similarity
def _common(a: Spectrum, b: Spectrum):
    """Resample both onto the overlapping x-grid of ``a``; return (ya, yb) or None."""
    lo = max(a.x[0], b.x[0])
    hi = min(a.x[-1], b.x[-1])
    if hi <= lo:
        return None
    gx = a.x[(a.x >= lo) & (a.x <= hi)]
    if gx.size < 3:
        return None
    return np.interp(gx, a.x, a.y), np.interp(gx, b.x, b.y)


def similarity_scores(a: Spectrum, b: Spectrum) -> dict:
    """Several similarity measures between two spectra (higher = more similar,
    except ``sam`` and ``euclid`` where lower = more similar)."""
    pair = _common(a, b)
    if pair is None:
        return {"correlation": 0.0, "cosine": 0.0, "sam": float("inf"),
                "euclid": float("inf")}
    ya, yb = pair
    # Pearson correlation (offset- and scale-invariant)
    da, db = ya - ya.mean(), yb - yb.mean()
    denom = np.sqrt((da @ da) * (db @ db))
    corr = float(da @ db / denom) if denom > _EPS else 0.0
    # cosine similarity & spectral angle (scale-invariant)
    na, nb = np.linalg.norm(ya), np.linalg.norm(yb)
    cos = float(ya @ yb / (na * nb)) if na > _EPS and nb > _EPS else 0.0
    sam = float(np.arccos(np.clip(cos, -1.0, 1.0)))      # radians
    # Euclidean distance of unit-vector-normalised spectra
    ua = ya / na if na > _EPS else ya
    ub = yb / nb if nb > _EPS else yb
    euclid = float(np.linalg.norm(ua - ub))
    return {"correlation": corr, "cosine": cos, "sam": sam, "euclid": euclid}


# ----------------------------------------------------------- the library
class SpectralLibrary:
    """A named collection of reference spectra."""

    def __init__(self, name: str = "library"):
        self.name = name
        self.entries: list[Spectrum] = []

    def __len__(self):
        return len(self.entries)

    def add(self, spec: Spectrum) -> None:
        self.entries.append(spec.copy())

    def remove_at(self, i: int) -> None:
        if 0 <= i < len(self.entries):
            del self.entries[i]

    def clear(self) -> None:
        self.entries.clear()

    # ---- persistence (JSON) ---------------------------------------------
    def save(self, path: str) -> None:
        obj = {"name": self.name, "library": [
            {"name": s.name, "x_unit": s.x_unit, "y_unit": s.y_unit,
             "x": s.x.tolist(), "y": s.y.tolist(),
             **({"laser_nm": s.meta["laser_nm"]} if s.meta.get("laser_nm") else {})}
            for s in self.entries]}
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(obj, fh, ensure_ascii=False, indent=1)

    @classmethod
    def load(cls, path: str) -> "SpectralLibrary":
        with open(path, "r", encoding="utf-8") as fh:
            obj = json.load(fh)
        lib = cls(name=obj.get("name") or os.path.splitext(os.path.basename(path))[0])
        for e in obj.get("library", []):
            meta = {"source": path}
            if e.get("laser_nm"):
                meta["laser_nm"] = e["laser_nm"]
            lib.entries.append(Spectrum(
                x=np.asarray(e["x"], float), y=np.asarray(e["y"], float),
                name=e.get("name", "ref"), x_unit=e.get("x_unit", "pixel"),
                y_unit=e.get("y_unit", "intensity"), meta=meta))
        return lib

    # ---- search ----------------------------------------------------------
    def search(self, unknown: Spectrum, top_n: int = 10,
               rank_by: str = "correlation") -> list[dict]:
        """Rank library entries by similarity to ``unknown``.

        Returns a list of dicts {entry, name, scores}, best first. ``rank_by`` is
        one of correlation/cosine (higher better) or sam/euclid (lower better).
        """
        results = []
        for entry in self.entries:
            scores = similarity_scores(unknown, entry)
            results.append({"entry": entry, "name": entry.name, "scores": scores})
        reverse = rank_by in ("correlation", "cosine")
        results.sort(key=lambda r: r["scores"].get(rank_by, 0.0), reverse=reverse)
        return results[:top_n]
