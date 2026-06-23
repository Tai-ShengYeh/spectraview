"""File IO: load by extension, save to CSV / JCAMP-DX.

``load_any`` is the single entry point used by the UI. It dispatches on the
file extension and returns a list of :class:`Spectrum`.
"""
from __future__ import annotations

import os

import numpy as np

from ..spectrum import Spectrum
from .ascii_io import load_ascii
from .binary_io import MissingDependency, load_opus, load_spc, load_sp
from .jcamp import load_jcamp
from .json_io import load_json, save_json
from .mat_io import load_mat

# Map file extensions to reader functions.
_READERS = {
    ".csv": load_ascii, ".txt": load_ascii, ".dat": load_ascii,
    ".tsv": load_ascii, ".asc": load_ascii, ".prn": load_ascii, ".xy": load_ascii,
    ".spectrum": load_ascii,   # NeoSpectra export: tab-separated x/y with a unit header
    ".dx": load_jcamp, ".jdx": load_jcamp, ".jcm": load_jcamp, ".jcamp": load_jcamp,
    ".spc": load_spc,
    ".sp": load_sp,   # PerkinElmer Spectrum (FTIR/UV-Vis)
    ".json": load_json,
    ".mat": load_mat,
}

# OPUS files often have numeric extensions (.0, .1, ...).
OPEN_FILTER = (
    "Spectra (*.csv *.txt *.dat *.tsv *.asc *.prn *.xy *.Spectrum *.dx *.jdx "
    "*.jcamp *.spc *.sp *.json *.mat);;"
    "ASCII (*.csv *.txt *.dat *.tsv *.asc *.prn *.xy);;"
    "NeoSpectra (*.Spectrum);;"
    "JCAMP-DX (*.dx *.jdx *.jcamp);;"
    "JSON (*.json);;"
    "MATLAB (*.mat);;"
    "GRAMS SPC (*.spc);;"
    "PerkinElmer (*.sp);;"
    "Bruker OPUS (*.0 *.1 *.2);;"
    "All files (*.*)"
)


def load_any(path: str) -> list[Spectrum]:
    """Load one file into a list of spectra, dispatching on its extension."""
    ext = os.path.splitext(path)[1].lower()
    reader = _READERS.get(ext)
    if reader is None:
        # OPUS files use numeric extensions like .0 / .1
        if ext[1:].isdigit():
            reader = load_opus
        else:
            reader = load_ascii  # last-resort: try to parse as ASCII
    return reader(path)


def save_csv(spec: Spectrum, path: str) -> None:
    """Write a single spectrum as a two-column CSV with a unit header."""
    header = f"{spec.x_label},{spec.y_label}"
    # Explicit UTF-8 handle: unit labels contain non-ASCII (e.g. cm⁻¹, µm),
    # which would crash np.savetxt under the Windows cp950 default codec.
    with open(path, "w", encoding="utf-8", newline="") as fh:
        np.savetxt(fh, np.column_stack([spec.x, spec.y]), delimiter=",",
                   header=header, comments="", fmt="%.6g")


def save_jcamp(spec: Spectrum, path: str) -> None:
    """Write a single spectrum as a JCAMP-DX 4.24 (XYPOINTS, AFFN) file."""
    xunit = {"cm-1": "1/CM", "nm": "NANOMETERS", "um": "MICROMETERS",
             "raman_cm-1": "1/CM", "eV": "EV"}.get(spec.x_unit, "ARBITRARY UNITS")
    yunit = {"absorbance": "ABSORBANCE", "transmittance": "TRANSMITTANCE",
             "%T": "TRANSMITTANCE", "reflectance": "REFLECTANCE",
             "KM": "KUBELKA-MUNK", "counts": "COUNTS"}.get(spec.y_unit, "ARBITRARY UNITS")
    lines = [
        f"##TITLE={spec.name}",
        "##JCAMP-DX=4.24",
        "##DATA TYPE=INFRARED SPECTRUM",
        "##ORIGIN=SpectraView",
        "##OWNER=public domain",
        f"##XUNITS={xunit}",
        f"##YUNITS={yunit}",
        f"##FIRSTX={spec.x[0]:.6g}",
        f"##LASTX={spec.x[-1]:.6g}",
        f"##NPOINTS={spec.npoints}",
        "##XFACTOR=1.0",
        "##YFACTOR=1.0",
        "##XYPOINTS=(XY..XY)",
    ]
    for xi, yi in zip(spec.x, spec.y):
        lines.append(f"{xi:.6g}, {yi:.6g}")
    lines.append("##END=")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def merge_spectra(spectra):
    """Put a set of spectra onto one common x-grid.

    Returns ``(common_x, [y_on_common, ...], x_label)``. If every spectrum
    already shares the same grid it is used verbatim; otherwise all spectra are
    interpolated onto a uniform grid spanning their overlapping x-range.
    """
    spectra = [s for s in spectra if s.npoints]
    if not spectra:
        raise ValueError("No spectra to export.")
    x0 = spectra[0].x
    same_grid = all(s.npoints == x0.size and np.allclose(s.x, x0) for s in spectra)
    if same_grid:
        common_x = x0
    else:
        lo = max(s.x[0] for s in spectra)
        hi = min(s.x[-1] for s in spectra)
        if hi <= lo:
            raise ValueError("Spectra do not share a common x-range; cannot merge.")
        n = max(s.npoints for s in spectra)
        common_x = np.linspace(lo, hi, n)
    ys = [np.interp(common_x, s.x, s.y) for s in spectra]
    return common_x, ys, spectra[0].x_label


def save_combined_csv(spectra, path: str, layout: str = "columns") -> dict:
    """Export several spectra merged into one CSV on a shared x-grid.

    layout='columns': first column = x, then one column per spectrum (re-imports
                      into SpectraView as multiple spectra).
    layout='rows':    one row per spectrum (an X-matrix for PLS/PCA: the first
                      row is the x-axis, the first column is the spectrum name).
    Returns a small info dict (n_spectra, n_points, resampled).
    """
    common_x, ys, x_label = merge_spectra(spectra)
    names = [s.name for s in spectra if s.npoints]
    resampled = not all(s.npoints == common_x.size and np.allclose(s.x, common_x)
                        for s in spectra if s.npoints)
    import csv
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        if layout == "rows":
            w.writerow([f"x:{x_label}"] + [f"{v:.8g}" for v in common_x])
            for name, y in zip(names, ys):
                w.writerow([name] + [f"{v:.8g}" for v in y])
        else:
            w.writerow([x_label] + names)
            ymat = np.column_stack(ys) if ys else np.empty((common_x.size, 0))
            for i in range(common_x.size):
                w.writerow([f"{common_x[i]:.8g}"]
                           + [f"{ymat[i, j]:.8g}" for j in range(ymat.shape[1])])
    return {"n_spectra": len(names), "n_points": int(common_x.size),
            "resampled": resampled}


__all__ = ["load_any", "save_csv", "save_jcamp", "save_json", "save_combined_csv",
           "merge_spectra", "OPEN_FILTER", "MissingDependency"]
