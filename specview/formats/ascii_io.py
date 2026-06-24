"""Reader for ASCII spectra: CSV / TXT / DAT / TSV / PRN.

Handles comma / tab / semicolon / whitespace delimiters, optional text headers,
European decimal commas, and multiple y-columns sharing a single x-column.
"""
from __future__ import annotations

import os
import re

import numpy as np

from ..spectrum import Spectrum

_NUMBER = re.compile(r"[-+]?\d")


def _looks_like_data(line: str) -> bool:
    """A data row starts with something numeric and has >=2 numeric-ish tokens."""
    s = line.strip()
    if not s or s.startswith(("#", "//", ";")):
        return False
    if not _NUMBER.match(s):
        return False
    # Tokens must split on the delimiter, so the comma is NOT part of a number.
    nums = re.findall(r"[-+]?\d[\d.eE+\-]*", s)
    return len(nums) >= 2


def _sniff_delimiter(line: str) -> str | None:
    for delim in (",", "\t", ";"):
        if delim in line:
            return delim
    if re.search(r"\s", line.strip()):
        return None  # whitespace split
    return None


def _guess_units(header_text: str):
    """Guess (x_unit, y_unit) from header keywords."""
    t = header_text.lower()
    x_unit = "pixel"
    if "raman" in t and "shift" in t:
        x_unit = "raman_cm-1"
    elif "cm-1" in t or "cm^-1" in t or "wavenumber" in t or "wavenumbers" in t:
        x_unit = "cm-1"
    elif "raman" in t:
        x_unit = "raman_cm-1"
    elif "nm" in t or "wavelength" in t:
        x_unit = "nm"
    elif "µm" in t or "micron" in t or "micromet" in t or re.search(r"\bum\b", t):
        x_unit = "um"   # \bum\b so "spectrum"/"aluminum" don't read as micrometres
    elif "ev" in t:
        x_unit = "eV"

    y_unit = "intensity"
    if "absorb" in t or re.search(r"\babs\b", t):
        y_unit = "absorbance"
    elif "%reflect" in t or "% reflect" in t or "%r " in t:
        y_unit = "%R"
    elif "reflect" in t:
        y_unit = "reflectance"
    elif "%transmit" in t or "% transmit" in t or "%t " in t or "%t" in t:
        y_unit = "%T"
    elif "transmit" in t:
        y_unit = "transmittance"
    elif "kubelka" in t or "k-m" in t or "k/m" in t:
        y_unit = "KM"
    elif "log(1/r)" in t or "log1r" in t:
        y_unit = "log1R"
    elif "count" in t:
        y_unit = "counts"
    return x_unit, y_unit


# PerkinElmer ASCII export ("PEDS") states its axes explicitly in a #GR block,
# so we read them instead of guessing from free-form header text.
_PE_ASC_XUNITS = {"CM-1": "cm-1", "NM": "nm", "NANOMETERS": "nm", "UM": "um",
                  "MICROMETERS": "um", "RAMAN": "raman_cm-1"}
_PE_ASC_YUNITS = {"A": "absorbance", "ABS": "absorbance", "ABSORBANCE": "absorbance",
                  "%T": "%T", "T": "transmittance", "TRANSMITTANCE": "transmittance",
                  "%R": "%R", "R": "reflectance", "KM": "KM", "LOG(1/R)": "log1R"}


def _is_perkinelmer_asc(first_line: str) -> bool:
    """PerkinElmer ASCII ('PEDS') files start e.g. 'PE IR ... ASCII PEDS 1.60'."""
    s = first_line.upper()
    return s.startswith("PE ") and "PEDS" in s


def _load_perkinelmer_asc(path: str, lines: list[str]) -> list[Spectrum]:
    """PerkinElmer ASCII export: units in the #GR block, x/y pairs after #DATA.

    More reliable than free-form unit guessing — the file declares its axes
    (e.g. CM-1 / A), so we read them rather than infer from surrounding text.
    """
    def _find(marker: str):
        return next((i for i, ln in enumerate(lines) if ln.strip() == marker), None)

    data_at = _find("#DATA")
    if data_at is None:
        raise ValueError(f"No #DATA block in {os.path.basename(path)}.")
    x_unit, y_unit = "cm-1", "absorbance"
    gr = _find("#GR")
    if gr is not None and gr + 2 < len(lines):
        x_unit = _PE_ASC_XUNITS.get(lines[gr + 1].strip().upper(), x_unit)
        y_unit = _PE_ASC_YUNITS.get(lines[gr + 2].strip().upper(), y_unit)
    xs, ys = [], []
    for ln in lines[data_at + 1:]:
        parts = ln.replace(",", " ").split()
        if len(parts) >= 2:
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
            except ValueError:
                continue
    if not xs:
        raise ValueError(f"No numeric data after #DATA in {os.path.basename(path)}.")
    base = os.path.splitext(os.path.basename(path))[0]
    return [Spectrum(x=np.asarray(xs, float), y=np.asarray(ys, float), name=base,
                     x_unit=x_unit, y_unit=y_unit, meta={"source": path})]


def load_ascii(path: str) -> list[Spectrum]:
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        raw_lines = fh.read().splitlines()
    if raw_lines and _is_perkinelmer_asc(raw_lines[0]):
        return _load_perkinelmer_asc(path, raw_lines)

    header_lines, data_lines = [], []
    started = False
    for ln in raw_lines:
        if not started and not _looks_like_data(ln):
            header_lines.append(ln)
        else:
            started = True
            if ln.strip():
                data_lines.append(ln)
    if not data_lines:
        raise ValueError(f"No numeric data found in {os.path.basename(path)}.")

    delim = _sniff_delimiter(data_lines[0])

    # European format heuristic: ';' delimiter usually means ',' is the decimal mark.
    euro_decimal = delim == ";"

    def parse_row(line: str):
        if euro_decimal:
            line = line.replace(",", ".")
        parts = line.split(delim) if delim else line.split()
        out = []
        for p in parts:
            p = p.strip()
            if p == "":
                continue
            try:
                out.append(float(p))
            except ValueError:
                out.append(np.nan)
        return out

    rows = [parse_row(ln) for ln in data_lines]
    ncol = max(len(r) for r in rows)
    rows = [r + [np.nan] * (ncol - len(r)) for r in rows if len(r) >= 2 or ncol == 1]
    arr = np.array(rows, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)

    header_text = "\n".join(header_lines)
    x_unit, y_unit = _guess_units(header_text)
    base = os.path.splitext(os.path.basename(path))[0].strip()

    # Try to use a header row of column names (the last header line) as legends.
    col_names = None
    if header_lines:
        last = header_lines[-1]
        d = _sniff_delimiter(last)
        toks = [t.strip() for t in (last.split(d) if d else last.split())]
        if len(toks) >= arr.shape[1]:
            col_names = toks

    x = arr[:, 0]
    spectra: list[Spectrum] = []
    ncols = arr.shape[1]
    for j in range(1, max(2, ncols)):
        if j >= ncols:
            break
        y = arr[:, j]
        if np.all(np.isnan(y)):
            continue
        # For a 2-column file the single column header is just the y-axis/unit
        # label, so use the filename. Only use column names to tell apart the
        # series of a genuinely multi-column file.
        name = base
        if ncols > 2:
            name = f"{base} [{j}]"
            if col_names and j < len(col_names) and col_names[j]:
                name = col_names[j]
        spectra.append(Spectrum(x=x, y=y, name=name, x_unit=x_unit, y_unit=y_unit,
                                meta={"source": path}))
    if not spectra:  # single column of y only -> use sample index as x
        spectra.append(Spectrum(x=np.arange(arr.shape[0]), y=arr[:, 0], name=base,
                                x_unit="pixel", y_unit=y_unit, meta={"source": path}))
    return spectra
