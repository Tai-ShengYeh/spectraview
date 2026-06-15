"""Fluorescence Excitation-Emission Matrix (EEM): model, readers, scatter removal.

An EEM holds excitation and emission wavelength axes and an intensity matrix
``Z`` of shape (n_ex, n_em): Z[i, j] is the intensity at ex[i], em[j].
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

import numpy as np

_EPS = 1e-12


@dataclass(eq=False)
class EEM:
    ex: np.ndarray            # excitation wavelengths (ascending)
    em: np.ndarray            # emission wavelengths (ascending)
    Z: np.ndarray             # intensity, shape (n_ex, n_em)
    name: str = "EEM"
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.ex = np.asarray(self.ex, dtype=float).ravel()
        self.em = np.asarray(self.em, dtype=float).ravel()
        self.Z = np.asarray(self.Z, dtype=float)
        if self.Z.shape != (self.ex.size, self.em.size):
            if self.Z.shape == (self.em.size, self.ex.size):
                self.Z = self.Z.T
            else:
                raise ValueError(
                    f"Z shape {self.Z.shape} doesn't match ex×em "
                    f"({self.ex.size}×{self.em.size}).")
        # sort both axes ascending
        oe, om = np.argsort(self.ex), np.argsort(self.em)
        self.ex, self.em = self.ex[oe], self.em[om]
        self.Z = self.Z[np.ix_(oe, om)]

    def copy(self) -> "EEM":
        return EEM(self.ex.copy(), self.em.copy(), self.Z.copy(), self.name,
                  dict(self.meta))

    @property
    def zrange(self):
        z = self.Z[np.isfinite(self.Z)]
        return (float(z.min()), float(z.max())) if z.size else (0.0, 1.0)


# ----------------------------------------------------------- readers
def _sniff(line: str):
    for d in (",", "\t", ";"):
        if d in line:
            return d
    return None


def read_eem_matrix(path: str, ex_in_columns: bool = True) -> EEM:
    """Read an EEM stored as a matrix CSV/TSV.

    The first row holds one wavelength axis and the first column the other; the
    top-left cell is a blank/label. ``ex_in_columns=True`` means the first row is
    the excitation axis (columns) and the first column is emission (rows).
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = [ln for ln in fh.read().splitlines() if ln.strip()]
    if len(lines) < 2:
        raise ValueError("EEM file has too few rows.")
    delim = _sniff(lines[0])

    def toks(ln):
        return ln.split(delim) if delim else ln.split()

    header = toks(lines[0])
    col_axis = np.array([float(t) for t in header[1:] if t.strip() != ""], dtype=float)
    row_axis, body = [], []
    for ln in lines[1:]:
        parts = toks(ln)
        if not parts or parts[0].strip() == "":
            continue
        try:
            row_axis.append(float(parts[0]))
        except ValueError:
            continue
        vals = []
        for p in parts[1:]:
            p = p.strip()
            try:
                vals.append(float(p))
            except ValueError:
                vals.append(np.nan)
        body.append(vals)
    row_axis = np.asarray(row_axis, dtype=float)
    ncol = col_axis.size
    body = np.array([r[:ncol] + [np.nan] * (ncol - len(r)) for r in body], dtype=float)

    name = os.path.splitext(os.path.basename(path))[0]
    if ex_in_columns:           # rows = emission, columns = excitation
        return EEM(ex=col_axis, em=row_axis, Z=body.T, name=name, meta={"source": path})
    return EEM(ex=row_axis, em=col_axis, Z=body, name=name, meta={"source": path})


def _excitation_of(spec, fallback: float) -> float:
    """Best-effort excitation wavelength for a 1-D emission spectrum."""
    if spec.meta.get("excitation_nm"):
        return float(spec.meta["excitation_nm"])
    m = re.search(r"(?:ex|EX|Ex)[ _]?(\d{2,4}(?:\.\d+)?)", spec.name)
    if m:
        return float(m.group(1))
    return float(fallback)


def eem_from_spectra(spectra, ex_values=None) -> EEM:
    """Assemble an EEM from several 1-D emission spectra (one per excitation)."""
    specs = [s for s in spectra if s.npoints]
    if len(specs) < 2:
        raise ValueError("Need at least two emission spectra to build an EEM.")
    if ex_values is not None and len(ex_values) == len(specs):
        exs = [float(v) for v in ex_values]
    else:
        exs = [_excitation_of(s, i) for i, s in enumerate(specs)]
    lo = max(s.x[0] for s in specs)
    hi = min(s.x[-1] for s in specs)
    if hi <= lo:
        raise ValueError("Emission spectra do not share a common range.")
    n = max(s.npoints for s in specs)
    em = np.linspace(lo, hi, n)
    Z = np.vstack([np.interp(em, s.x, s.y) for s in specs])   # (n_ex, n_em)
    ex = np.asarray(exs, dtype=float)
    return EEM(ex=ex, em=em, Z=Z, name="EEM (from spectra)")


# ----------------------------------------------------------- scatter removal
def remove_scatter(eem: EEM, first: bool = True, second: bool = True,
                   raman: bool = False, width: float = 15.0,
                   raman_shift: float = 3400.0, fill: str = "nan") -> EEM:
    """Mask Rayleigh (em≈ex, em≈2·ex) and optional Raman scatter bands."""
    out = eem.copy()
    EX = out.ex[:, None]          # (n_ex, 1)
    EM = out.em[None, :]          # (1, n_em)
    mask = np.zeros(out.Z.shape, dtype=bool)
    if first:
        mask |= np.abs(EM - EX) <= width
    if second:
        mask |= np.abs(EM - 2.0 * EX) <= width
    if raman:
        # water Raman emission wavelength for each excitation (nm)
        em_raman = 1.0 / (1.0 / np.clip(EX, _EPS, None) - raman_shift * 1e-7)
        mask |= np.abs(EM - em_raman) <= width
    fillval = 0.0 if fill == "zero" else np.nan
    out.Z[mask] = fillval
    out.meta["scatter_removed"] = True
    return out
