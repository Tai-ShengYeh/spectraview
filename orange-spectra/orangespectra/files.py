"""Bulk spectra loading: files, folders, .zip archives, and NetCDF (.cdf).

Qt-free, so it's testable headlessly. Text formats reuse the URL parsers in
:mod:`orangespectra.core` (JCAMP-DX AFFN, two-column CSV) and add the two
matrix layouts SpectraView's "combined export" writes:

  * wide:  header = numeric wavelengths, each data row = one spectrum
           (optionally a leading name column)
  * long:  first column = wavelength, each further column = one spectrum
           (header names the spectra)

NetCDF (.cdf, classic/NetCDF3 — e.g. chemometrics datasets like applewine, or
ANDI chromatography exports) is read with ``scipy.io.netcdf_file`` (no extra
dependency): the largest 2-D numeric variable becomes rows = spectra, with a
matching 1-D variable as the x-axis when one exists.
"""
from __future__ import annotations

import io
import os
import re
import zipfile

import numpy as np

from .core import make_spectrum, parse_csv, parse_jcamp, split_fields

TEXT_EXTS = {".dx", ".jdx", ".jcm", ".csv", ".txt", ".tsv", ".dat", ".asc",
             ".xy", ".prn"}
NETCDF_EXTS = {".cdf", ".nc"}
SPECTRA_EXTS = TEXT_EXTS | NETCDF_EXTS

_NUMRE = re.compile(r"^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$")


def _is_num(tok: str) -> bool:
    return bool(_NUMRE.match(tok.strip()))


# ============================================================ text formats
def parse_matrix_text(text: str, source: str = ""):
    """Parse the two matrix CSV layouts into a list of spectra (or None)."""
    rows = []
    for line in text.splitlines():
        parts = split_fields(line)
        if len(parts) >= 3:
            rows.append(parts)
    if len(rows) < 2:
        return None
    header, data = rows[0], rows[1:]
    stem = os.path.splitext(os.path.basename(source))[0] or "spectrum"

    # wide layout: header (except maybe first cell) is all numeric = x axis
    hdr_num_from = 1 if not _is_num(header[0]) else 0
    if all(_is_num(t) for t in header[hdr_num_from:]) and \
            len(header) - hdr_num_from >= 3:
        x = np.array([float(t) for t in header[hdr_num_from:]])
        specs = []
        for i, r in enumerate(data):
            has_name = not _is_num(r[0])
            name = r[0] if has_name else f"{stem} {i + 1}"
            vals = r[1:] if has_name else r
            if len(vals) != x.size or not all(_is_num(v) for v in vals):
                continue
            specs.append(make_spectrum(x, [float(v) for v in vals],
                                       name=name, source=source))
        if specs:
            return specs

    # long layout: first column numeric x; other columns = spectra
    if all(_is_num(r[0]) for r in data):
        ncol = min(len(r) for r in data)
        if ncol >= 2 and all(all(_is_num(v) for v in r[:ncol]) for r in data):
            x = np.array([float(r[0]) for r in data])
            names = (header[1:ncol] if not _is_num(header[1])
                     else [f"{stem} {i}" for i in range(1, ncol)])
            return [make_spectrum(x, [float(r[c]) for r in data],
                                  name=str(names[c - 1]), source=source)
                    for c in range(1, ncol)]
    return None


def parse_text_spectra(text: str, source: str = ""):
    """Text content -> list of spectra: JCAMP, matrix CSV, or two-column CSV."""
    if "##TITLE" in text.upper():
        spec = parse_jcamp(text, source=source)
        return [spec] if spec is not None else []
    specs = parse_matrix_text(text, source=source)
    if specs:
        return specs
    spec = parse_csv(text, source=source)
    return [spec] if spec is not None else []


# ================================================================ NetCDF
_X_NAMES = ("wavelength", "wavenumber", "wave", "lambda", "x", "axis",
            "scan_acquisition_time", "retention_time", "time", "mz")


def load_netcdf(path_or_bytes, source: str = ""):
    """Read spectra from a classic NetCDF (.cdf/.nc) file.

    Strategy: the largest 2-D numeric variable is the data matrix (rows =
    spectra/samples; transposed if that matches an x-candidate). The x-axis is
    a 1-D variable whose length equals the number of columns (preferring
    wavelength-ish names); sample names come from a matching char matrix, else
    are numbered. With no 2-D variable, equal-length 1-D pairs are used
    (monotonic one = x).
    """
    from scipy.io import netcdf_file

    if isinstance(path_or_bytes, bytes):
        fh = netcdf_file(io.BytesIO(path_or_bytes), "r", mmap=False)
        stem = os.path.splitext(os.path.basename(source))[0] or "netcdf"
    else:
        fh = netcdf_file(path_or_bytes, "r", mmap=False)
        stem = os.path.splitext(os.path.basename(path_or_bytes))[0]
        source = source or path_or_bytes
    try:
        num_vars, char_vars = {}, {}
        for name, var in fh.variables.items():
            arr = np.asarray(var[:])
            if arr.dtype.kind in "fiu":
                num_vars[name] = arr.astype(float)
            elif arr.dtype.kind in "SU":
                char_vars[name] = arr

        two_d = {k: v for k, v in num_vars.items()
                 if v.ndim == 2 and min(v.shape) >= 1 and max(v.shape) >= 3}
        one_d = {k: v for k, v in num_vars.items() if v.ndim == 1 and v.size >= 3}

        if two_d:
            key = max(two_d, key=lambda k: two_d[k].size)
            X = two_d[key]
            xvec, xname = _pick_x(one_d, X.shape[1])
            if xvec is None:
                alt, altname = _pick_x(one_d, X.shape[0])
                if alt is not None:
                    X, xvec, xname = X.T, alt, altname
            if xvec is None:
                xvec = np.arange(X.shape[1], dtype=float)
                xname = "index"
            names = _sample_names(char_vars, X.shape[0], stem)
            return [make_spectrum(xvec, X[i], name=names[i],
                                  x_label=xname, source=source)
                    for i in range(X.shape[0])]

        # 1-D fallback: (x, y) variable pairs of equal length
        if len(one_d) >= 2:
            items = sorted(one_d.items(), key=lambda kv: -kv[1].size)
            for xk, xv in items:
                if np.all(np.diff(xv) > 0) or np.all(np.diff(xv) < 0):
                    ys = [(k, v) for k, v in items
                          if k != xk and v.size == xv.size]
                    if ys:
                        return [make_spectrum(xv, yv, name=f"{stem}:{k}",
                                              x_label=xk, source=source)
                                for k, yv in ys]
        raise ValueError(
            f"No spectra-like variables found in {source or 'NetCDF file'} "
            f"(variables: {', '.join(fh.variables) or 'none'}).")
    finally:
        fh.close()


def _pick_x(one_d: dict, n: int):
    cands = {k: v for k, v in one_d.items() if v.size == n}
    if not cands:
        return None, None
    for pref in _X_NAMES:
        for k in cands:
            if pref in k.lower():
                return cands[k], k
    # otherwise prefer a monotonic one
    for k, v in cands.items():
        if np.all(np.diff(v) > 0) or np.all(np.diff(v) < 0):
            return v, k
    k = next(iter(cands))
    return cands[k], k


def _sample_names(char_vars: dict, n: int, stem: str):
    for arr in char_vars.values():
        if arr.ndim == 2 and arr.shape[0] == n:
            try:
                return [bytes(row).decode("latin-1").strip("\x00 ").strip()
                        or f"{stem} {i + 1}" for i, row in enumerate(arr)]
            except Exception:                          # noqa: BLE001
                break
    return [f"{stem} {i + 1}" for i in range(n)]


# ======================================================= files / folders / zip
def load_spectra_path(path: str) -> list:
    """One file path -> list of spectra (dispatch on extension)."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".zip":
        return load_spectra_zip(path)
    if ext in NETCDF_EXTS:
        return load_netcdf(path)
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return parse_text_spectra(fh.read(), source=path)


def load_spectra_zip(path: str) -> list:
    """All spectra files inside a .zip archive -> list of spectra."""
    specs = []
    with zipfile.ZipFile(path) as zf:
        for info in zf.infolist():
            if info.is_dir():
                continue
            ext = os.path.splitext(info.filename)[1].lower()
            if ext not in SPECTRA_EXTS:
                continue
            raw = zf.read(info)
            src = f"{path}!{info.filename}"
            try:
                if ext in NETCDF_EXTS:
                    specs.extend(load_netcdf(raw, source=src))
                else:
                    specs.extend(parse_text_spectra(
                        raw.decode("utf-8", errors="replace"), source=src))
            except (ValueError, OSError):
                continue
    return specs


def load_spectra_folder(path: str, recursive: bool = True) -> list:
    """All spectra files in a folder (optionally recursive) -> list of spectra."""
    specs = []
    if recursive:
        walker = os.walk(path)
    else:
        walker = [(path, [], sorted(os.listdir(path)))]
    for root, _dirs, files in walker:
        for fn in sorted(files):
            ext = os.path.splitext(fn)[1].lower()
            if ext not in SPECTRA_EXTS and ext != ".zip":
                continue
            try:
                specs.extend(load_spectra_path(os.path.join(root, fn)))
            except (ValueError, OSError):
                continue
    return specs
