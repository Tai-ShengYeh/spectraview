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


@dataclass
class ParafacResult:
    """PARAFAC decomposition of a stack of EEMs into ``rank`` components."""
    ex: np.ndarray            # excitation axis
    em: np.ndarray            # emission axis
    ex_load: np.ndarray       # (n_ex, rank) excitation loadings
    em_load: np.ndarray       # (n_em, rank) emission loadings
    scores: np.ndarray        # (n_samples, rank) sample scores
    fit: float                # 1 - ||X - X̂|| / ||X||
    names: list               # sample names

    @property
    def rank(self) -> int:
        return self.scores.shape[1]

    def component_eem(self, f: int) -> "EEM":
        """The rank-1 EEM of component f (excitation ⊗ emission loading)."""
        Z = np.outer(self.ex_load[:, f], self.em_load[:, f])
        return EEM(self.ex, self.em, Z, name=f"PARAFAC comp {f + 1}")


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
# ----------------------------------------------------------- PARAFAC
def _khatri_rao(B: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Column-wise Khatri-Rao product: KR[j*K+k, f] = B[j,f]·C[k,f]."""
    return (B[:, None, :] * C[None, :, :]).reshape(B.shape[0] * C.shape[0], B.shape[1])


def parafac(X: np.ndarray, rank: int, n_iter: int = 400, tol: float = 1e-9,
            nonneg: bool = True, seed: int = 0):
    """PARAFAC (CANDECOMP) of a 3-way array via alternating least squares.

    X has shape (n_samples, n_ex, n_em). Returns (scores, ex_load, em_load, fit)
    with X[i,j,k] ≈ Σ_f scores[i,f]·ex_load[j,f]·em_load[k,f]. ``nonneg`` projects
    every factor to be non-negative (appropriate for fluorescence).
    """
    X = np.asarray(X, dtype=float)
    I, J, K = X.shape
    rng = np.random.default_rng(seed)
    A = rng.random((I, rank))
    B = rng.random((J, rank))
    C = rng.random((K, rank))
    X0 = X.reshape(I, J * K)
    X1 = np.moveaxis(X, 1, 0).reshape(J, I * K)
    X2 = np.moveaxis(X, 2, 0).reshape(K, I * J)
    normX = np.linalg.norm(X) or 1.0
    prev = None
    err = 1.0
    for _ in range(n_iter):
        A = X0 @ _khatri_rao(B, C) @ np.linalg.pinv((B.T @ B) * (C.T @ C))
        if nonneg:
            A = np.clip(A, 0.0, None)
        B = X1 @ _khatri_rao(A, C) @ np.linalg.pinv((A.T @ A) * (C.T @ C))
        if nonneg:
            B = np.clip(B, 0.0, None)
        C = X2 @ _khatri_rao(A, B) @ np.linalg.pinv((A.T @ A) * (B.T @ B))
        if nonneg:
            C = np.clip(C, 0.0, None)
        Xhat = (A @ _khatri_rao(B, C).T).reshape(I, J, K)
        err = float(np.linalg.norm(X - Xhat) / normX)
        if prev is not None and abs(prev - err) < tol:
            break
        prev = err
    # resolve scaling: unit-norm excitation & emission loadings, push size to scores
    for f in range(rank):
        nb = np.linalg.norm(B[:, f]) or 1.0
        nc = np.linalg.norm(C[:, f]) or 1.0
        B[:, f] /= nb
        C[:, f] /= nc
        A[:, f] *= nb * nc
    return A, B, C, 1.0 - err


def parafac_from_eems(eems, rank: int, nonneg: bool = True) -> ParafacResult:
    """Run PARAFAC on a list of EEMs that share the same ex/em grid."""
    eems = list(eems)
    if len(eems) < 2:
        raise ValueError("PARAFAC needs at least two EEMs (samples).")
    ref = eems[0]
    for e in eems:
        if (e.ex.shape != ref.ex.shape or e.em.shape != ref.em.shape
                or not np.allclose(e.ex, ref.ex) or not np.allclose(e.em, ref.em)):
            raise ValueError("All EEMs must share the same excitation/emission grid.")
    X = np.stack([np.nan_to_num(e.Z, nan=0.0) for e in eems])   # (n_samples, n_ex, n_em)
    if rank < 1 or rank > min(ref.ex.size, ref.em.size, len(eems)):
        raise ValueError("Number of components is out of range for this data.")
    A, B, C, fit = parafac(X, rank, nonneg=nonneg)
    # order components by descending mean score (largest first)
    order = np.argsort(A.mean(axis=0))[::-1]
    return ParafacResult(ex=ref.ex, em=ref.em, ex_load=B[:, order], em_load=C[:, order],
                         scores=A[:, order], fit=float(fit),
                         names=[e.name for e in eems])


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
