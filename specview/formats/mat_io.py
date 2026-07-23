"""Reader for MATLAB .mat files.

* v4 / v5 / v7 ("-v7"):  via scipy.io (no extra dependency)
* v7.3 ("-v7.3", HDF5):  via h5py, when installed (``pip install h5py``)

Resolves x / y by, in order:
  1. variables whose names match x/y aliases (x, wavelength, intensity, ...)
  2. a single 2-D matrix: N x 2 / 2 x N -> (x, y); N x M (M>2) -> col0 = x, rest = spectra
  3. exactly two 1-D numeric arrays of equal length -> (first = x, second = y)
     ...or one lone vector -> plotted against sample index

MATLAB *objects* saved from PLS_Toolbox / Eigenvector "dataset" class (e.g.
the classic corn.mat benchmark) are also handled: each object's ``data``
matrix is unpacked row-by-row, with the wavelength axis pulled from the
matching ``axisscale`` cell.
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


def _unwrap(a):
    """Peel 0-d object arrays (MatlabObject nesting) down to the payload."""
    while isinstance(a, np.ndarray) and a.shape == ():
        a = a.item()
    return a


def _guess_xunit(x: np.ndarray) -> str:
    """Crude axis-unit guess for dataset objects lacking a labelled axis."""
    lo, hi = float(np.min(x)), float(np.max(x))
    if 100.0 <= lo and hi <= 3200.0:
        return "nm"          # UV/Vis/NIR wavelengths
    if 3200.0 < lo and hi <= 13000.0:
        return "cm-1"        # NIR/IR wavenumbers
    return "pixel"


def _dataset_spectra(key: str, obj, base: str, path: str):
    """Unpack one MATLAB object with Eigenvector dataset-like fields.

    Returns a list of Spectrum, or None if the object doesn't look spectral
    (e.g. a property-value table without a matching axisscale).
    """
    names = getattr(getattr(obj, "dtype", None), "names", None)
    if not names or "data" not in names:
        return None
    try:
        d = np.asarray(_unwrap(obj["data"]), dtype=float)
    except (TypeError, ValueError):
        return None
    if d.ndim == 1:
        d = d[None, :]
    if d.ndim != 2 or d.size == 0:
        return None

    # find a wavelength axis among the axisscale cells
    x = None
    if "axisscale" in names:
        try:
            cells = np.asarray(_unwrap(obj["axisscale"]), dtype=object).ravel()
        except (TypeError, ValueError):
            cells = []
        for want, transpose in ((d.shape[1], False), (d.shape[0], True)):
            for c in cells:
                try:
                    v = np.asarray(_unwrap(c), dtype=float).squeeze()
                except (TypeError, ValueError):
                    continue
                if v.ndim == 1 and v.size == want and v.size > 1:
                    x = v
                    if transpose:
                        d = d.T
                    break
            if x is not None:
                break
    if x is None:
        return None          # not spectral (no axis that matches the data)
    if x.size > 2 and np.array_equal(x, np.arange(x[0], x[0] + x.size)):
        return None          # axis is a plain 1..N index -> a table, not spectra

    label = key
    if "name" in names:
        nm = _unwrap(obj["name"])
        if isinstance(nm, str) and nm.strip():
            label = nm.strip()
    unit = _guess_xunit(x)
    many = d.shape[0] > 1
    return [Spectrum(x, d[i],
                     name=f"{base}:{label}" + (f" [{i}]" if many else ""),
                     x_unit=unit, meta={"source": path, "variable": key})
            for i in range(d.shape[0])]


def _spectra_from_vars(data: dict, base: str, path: str) -> list[Spectrum]:
    """Shared x/y resolution for both the scipy and h5py loaders."""
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
        # named y is a matrix (rows or columns of spectra) with a named x vector
        if x is not None:
            ym = np.asarray(data[yk], dtype=float)
            if ym.ndim == 2:
                if ym.shape[0] == x.size:      # spectra in columns
                    ym = ym.T
                if ym.shape[1] == x.size:
                    return [Spectrum(x, ym[i], name=f"{base} [{i}]",
                                     x_unit=H.norm_xunit(xk) or "pixel",
                                     meta={"source": path})
                            for i in range(ym.shape[0])]

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
        "Couldn't identify x/y arrays in %s. Variables found: %s."
        % (os.path.basename(path), ", ".join(keys) or "(none)"))


def _load_mat73(path: str) -> dict:
    """Extract numeric variables from a MATLAB v7.3 (HDF5) file via h5py."""
    try:
        import h5py
    except ImportError as exc:
        raise MissingDependency(
            "This is a MATLAB v7.3 (HDF5) .mat file. Install h5py to read it:\n"
            "    pip install h5py\n"
            "(or re-save it in MATLAB as v7:  save('file.mat','var','-v7') )"
        ) from exc

    data: dict = {}
    with h5py.File(path, "r") as f:
        def visit(name, obj):
            if not isinstance(obj, h5py.Dataset):
                return
            if name.startswith("#refs#"):          # cell/struct internals
                return
            mcls = obj.attrs.get("MATLAB_class", b"")
            if isinstance(mcls, bytes):
                mcls = mcls.decode("ascii", "ignore")
            if mcls in ("char", "cell", "struct", "function_handle", "logical"):
                return
            try:
                arr = np.asarray(obj[()], dtype=float)
            except (TypeError, ValueError):
                return
            # MATLAB v7.3 stores arrays column-major -> transpose back
            if arr.ndim >= 2:
                arr = arr.T
            key = name.split("/")[-1]
            data[key] = np.squeeze(arr)
        f.visititems(visit)
    return data


def load_mat(path: str) -> list[Spectrum]:
    from scipy.io import loadmat
    base = os.path.splitext(os.path.basename(path))[0]
    try:
        raw = loadmat(path, squeeze_me=True)
        data = {k: v for k, v in raw.items() if not k.startswith("__")}
    except NotImplementedError:
        # scipy can't read v7.3 (HDF5) -> h5py fallback
        data = _load_mat73(path)
    except Exception as exc:  # noqa: BLE001
        raise ValueError("Could not read %s: %s" % (os.path.basename(path), exc)) from exc

    # PLS_Toolbox / Eigenvector dataset objects (corn.mat and friends)
    from_objects = []
    plain = {}
    for k, v in data.items():
        specs = _dataset_spectra(k, v, base, path)
        if specs:
            from_objects.extend(specs)
        elif getattr(getattr(v, "dtype", None), "names", None):
            continue          # non-spectral struct/object (e.g. property table)
        else:
            plain[k] = v
    if from_objects:
        return from_objects
    return _spectra_from_vars(plain or data, base, path)
