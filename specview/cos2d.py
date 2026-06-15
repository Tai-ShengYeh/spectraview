"""Generalized two-dimensional correlation spectroscopy (2D-COS), Noda.

Given a perturbation series of spectra (rows = observations, columns =
wavenumbers), computes the synchronous (Φ) and asynchronous (Ψ) correlation
maps. Also implements two-trace 2D (2T2D) correlation for exactly two spectra.

Conventions follow shigemorita's 2Dpy:
    Φ = Dᵀ·D / (m-1)
    Ψ = Dᵀ·(H·D) / (m-1)
where D is the dynamic spectra (mean-subtracted) and H is the Hilbert-Noda
matrix  H[i,j] = 1/(π·(j-i))  for i≠j, 0 on the diagonal.
"""
from __future__ import annotations

import numpy as np

from .spectrum import Spectrum

_EPS = 1e-12


def hilbert_noda_matrix(m: int) -> np.ndarray:
    """The m×m Hilbert-Noda transformation matrix."""
    H = np.zeros((m, m))
    idx = np.arange(m)
    for i in range(m):
        d = idx - i
        nz = d != 0
        H[i, nz] = 1.0 / (np.pi * d[nz])
    return H


def correlation(M: np.ndarray, ref: str = "mean"):
    """Synchronous & asynchronous maps from a data matrix.

    M: shape (m_observations, n_wavenumbers). Returns (sync, asyn), each n×n.
    ref: 'mean' subtracts the average spectrum (dynamic spectra); 'none' keeps M.
    """
    M = np.asarray(M, dtype=float)
    if M.ndim != 2 or M.shape[0] < 2:
        raise ValueError("Need at least two spectra (rows) for 2D-COS.")
    D = M - M.mean(axis=0) if ref == "mean" else M
    m = D.shape[0]
    sync = (D.T @ D) / (m - 1)
    asyn = (D.T @ (hilbert_noda_matrix(m) @ D)) / (m - 1)
    return sync, asyn


def hetero_correlation(M1: np.ndarray, M2: np.ndarray, ref: str = "mean"):
    """Hetero-spectral 2D correlation between two datasets sharing a perturbation.

    M1 (m × n1) and M2 (m × n2) are two techniques measured on the SAME m
    perturbation points. Returns (sync, asyn), each n1 × n2.
    """
    M1 = np.asarray(M1, dtype=float)
    M2 = np.asarray(M2, dtype=float)
    if M1.shape[0] != M2.shape[0]:
        raise ValueError("Both datasets need the same number of perturbation rows.")
    m = M1.shape[0]
    if m < 2:
        raise ValueError("Need at least two perturbation points.")
    D1 = M1 - M1.mean(axis=0) if ref == "mean" else M1
    D2 = M2 - M2.mean(axis=0) if ref == "mean" else M2
    sync = (D1.T @ D2) / (m - 1)
    asyn = (D1.T @ (hilbert_noda_matrix(m) @ D2)) / (m - 1)
    return sync, asyn


def hetero_from_spectra(group1: list[Spectrum], group2: list[Spectrum],
                        ref: str = "mean"):
    """Hetero-correlation from two equal-length groups of spectra.

    Returns (x1, x2, sync, asyn); sync/asyn have shape (n1, n2).
    """
    if len(group1) != len(group2):
        raise ValueError("The two sets must contain the same number of spectra.")
    if len(group1) < 2:
        raise ValueError("Each set needs at least two spectra.")
    g1, M1 = _common_grid(group1)
    g2, M2 = _common_grid(group2)
    sync, asyn = hetero_correlation(M1, M2, ref)
    return g1, g2, sync, asyn


def two_trace(ya: np.ndarray, yb: np.ndarray):
    """Two-trace 2D (2T2D) correlation, Noda 2018.

    Works from exactly two spectra (no mean-centering, so the asynchronous map
    is non-degenerate):
        Φ = ½(yₐ⊗yₐ + y_b⊗y_b)
        Ψ = ½(yₐ⊗y_b − y_b⊗yₐ)
    """
    ya = np.asarray(ya, dtype=float)
    yb = np.asarray(yb, dtype=float)
    sync = 0.5 * (np.outer(ya, ya) + np.outer(yb, yb))
    asyn = 0.5 * (np.outer(ya, yb) - np.outer(yb, ya))
    return sync, asyn


# ---------------------------------------------------- from Spectrum objects
def _common_grid(spectra: list[Spectrum]):
    lo = max(s.x[0] for s in spectra)
    hi = min(s.x[-1] for s in spectra)
    if hi <= lo:
        raise ValueError("Spectra do not share a common x-range.")
    base = spectra[0]
    gx = base.x[(base.x >= lo) & (base.x <= hi)]
    if gx.size < 3:
        raise ValueError("Too little overlap between spectra.")
    M = np.vstack([np.interp(gx, s.x, s.y) for s in spectra])
    return gx, M


def correlation_from_spectra(spectra: list[Spectrum], ref: str = "mean"):
    """Resample a list of spectra to a common grid and compute 2D-COS.

    Returns (x, sync, asyn) where x is the shared wavenumber/axis grid.
    """
    if len(spectra) < 2:
        raise ValueError("2D-COS needs at least two spectra.")
    gx, M = _common_grid(spectra)
    sync, asyn = correlation(M, ref)
    return gx, sync, asyn


def two_trace_from_spectra(spec_a: Spectrum, spec_b: Spectrum):
    """2T2D correlation from two spectra (resampled to their overlap)."""
    gx, M = _common_grid([spec_a, spec_b])
    sync, asyn = two_trace(M[0], M[1])
    return gx, sync, asyn
