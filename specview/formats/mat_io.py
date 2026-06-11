"""Reader for MATLAB .mat files (via scipy.io — no extra dependency).

Resolves x / y by, in order:
  1. variables whose names match x/y aliases (x, wavelength, intensity, ...)
  2. a single 2-D matrix: N×2 / 2×N -> (x, y); N×M (M>2) -> col0 = x, rest = spectra
  3. exactly two 1-D numeric arrays of equal length -> (first = x, second = y)

MATLAB v7.3 files are HDF5; scipy can't read them, so a clear hint is raised.
"""
from __future__ import annotations

import os

import numpy as np

from ..spectrum import Spectrum
from . import _hints as H
from .binary_io import MissingDependency


def _vector(v):
    """Return a 1-D float array if v squeezes to one, else None."""
    a = np.asarray(v, dtype=float).squeeze()
    return a if a.ndim == 1 and a.size > 1 else None


def load_mat(path: str) -> list[Spectrum]:
    from scipy.io import loadmat
    try:
        raw = loadmat(path, squeeze_me=True)
    except NotImplementedError as exc:
        raise MissingDependency(
            "This is a MATLAB v7.3 (HDF5) .mat file, which SciPy can't read.\n"
            "Re-save it in MATLAB as v7:  save('file.mat', 'var', '-v7')\n"
            "(or install h5py and export to another format)."
        ) from exc
    except Exception as exc:  # noqa: BLE001
        raise ValueError(f"Could not read {os.path.basename(path)}: {exc}") from exc

    data = {k: v for k, v in raw.items() if not k.startswith("__")}
    base = os.path.splitext(os.path.basename(path))[0]
    keys = list(data.keys())

    # 1) named x / y variables
    xk = H.match_key(keys, H.X_KEYS)
    yk = H.match_key(keys, H.Y_KEYS)
    if xk and yk:
        x, y = _vector(data[xk]), _vector(data[yk])
        if x is not None and y is not None:
            return [Spectrum(x, y, name=base, x_unit=H.norm_xunit(xk) or "pixel",
                             y_unit=H.norm_yunit(yk) or "intensity",
                             meta={"source": path})]

    # 2) a single 2-D matrix
    matrices = [(k, np.asarray(v, dtype=float)) for k, v in data.items()
                if np.asarray(v).ndim == 2 and np.asarray(v).size >= 4]
    if len(matrices) == 1:
        k, arr = matrices[0]
        if arr.shape[0] < arr.shape[1]:      # wider than tall -> rows are series
            arr = arr.T
        if arr.shape[1] == 2:
            return [Spectrum(arr[:, 0], arr[:, 1], name=base, meta={"source": path})]
        if arr.shape[1] > 2:                  # col 0 = x, remaining = spectra
            x = arr[:, 0]
            return [Spectrum(x, arr[:, j], name=f"{base} [{j}]", meta={"source": path})
                    for j in range(1, arr.shape[1])]

    # 3) exactly two equal-length 1-D arrays
    vectors = [(k, _vector(v)) for k, v in data.items()]
    vectors = [(k, a) for k, a in vectors if a is not None]
    if len(vectors) == 2 and vectors[0][1].size == vectors[1][1].size:
        (xk, x), (yk, y) = vectors
        return [Spectrum(x, y, name=base, x_unit=H.norm_xunit(xk) or "pixel",
                         y_unit=H.norm_yunit(yk) or "intensity",
                         meta={"source": path})]
    # ...or one lone vector -> use sample index as x
    if len(vectors) == 1:
        y = vectors[0][1]
        return [Spectrum(np.arange(y.size, dtype=float), y, name=base,
                         x_unit="pixel", meta={"source": path})]

    raise ValueError(
        f"Couldn't identify x/y arrays in {os.path.basename(path)}. "
        f"Variables found: {', '.join(keys) or '(none)'}.")
