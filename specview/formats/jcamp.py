"""A self-contained JCAMP-DX reader (no third-party dependency).

Supports the common cases used by FTIR / UV-Vis / Raman instruments:
  * ##XYDATA= (X++(Y..Y))           with AFFN, PAC, SQZ, DIF and DUP encodings
  * ##XYPOINTS= / ##PEAK TABLE=     as (X, Y) pairs (AFFN)
  * multi-block (LINK) files        -> one Spectrum per data block

ASDF reference: McDonald & Wilks, JCAMP-DX 4.24 / 5.00.
"""
from __future__ import annotations

import os
import re

import numpy as np

from ..spectrum import Spectrum

# --- ASDF pseudo-digit tables ------------------------------------------------
_SQZ = {"@": 0, "A": 1, "B": 2, "C": 3, "D": 4, "E": 5, "F": 6, "G": 7, "H": 8, "I": 9,
        "a": -1, "b": -2, "c": -3, "d": -4, "e": -5, "f": -6, "g": -7, "h": -8, "i": -9}
_DIF = {"%": 0, "J": 1, "K": 2, "L": 3, "M": 4, "N": 5, "O": 6, "P": 7, "Q": 8, "R": 9,
        "j": -1, "k": -2, "l": -3, "m": -4, "n": -5, "o": -6, "p": -7, "q": -8, "r": -9}
_DUP = {"S": 1, "T": 2, "U": 3, "V": 4, "W": 5, "X": 6, "Y": 7, "Z": 8, "s": 9}

# Chars that *unambiguously* mark a compressed (ASDF) data section.
# 'E'/'e' are excluded because they double as exponent markers in plain AFFN.
_COMPRESSED_CHARS = set("@%ABCDFGHIJKLMNOPQRSTUVWXYZabcdfghijklmnopqrs")

_XUNIT_MAP = {
    "1/CM": "cm-1", "1/ CM": "cm-1", "CM-1": "cm-1", "CM^-1": "cm-1",
    "NANOMETERS": "nm", "NANOMETER": "nm", "NM": "nm",
    "MICROMETERS": "um", "MICROMETER": "um", "MICRONS": "um",
    "ELECTRON VOLT": "eV", "EV": "eV", "RAMANSHIFT": "raman_cm-1",
}
_YUNIT_MAP = {
    "ABSORBANCE": "absorbance", "TRANSMITTANCE": "transmittance",
    "%TRANSMITTANCE": "%T", "% TRANSMITTANCE": "%T", "%T": "%T",
    "REFLECTANCE": "reflectance", "%REFLECTANCE": "%R",
    "KUBELKA-MUNK": "KM", "KUBELKAMUNK": "KM",
    "ARBITRARY UNITS": "a.u.", "ARBITRARY": "a.u.", "COUNTS": "counts",
    "INTENSITY": "intensity",
}


def _decode_line(line: str, compressed: bool):
    """Decode one data line into a list of (mode, value).

    mode is one of 'A' (absolute), 'D' (difference) or 'P' (dup-count).
    """
    tokens = []
    i, n = 0, len(line)
    while i < n:
        c = line[i]
        if c in " ,\t":
            i += 1
            continue
        mode = "A"
        digits = ""
        if compressed and c in _SQZ:
            v = _SQZ[c]
            digits = str(abs(v))
            sign = -1 if v < 0 else 1
            i += 1
        elif compressed and c in _DIF:
            v = _DIF[c]
            digits = str(abs(v))
            sign = -1 if v < 0 else 1
            mode = "D"
            i += 1
        elif compressed and c in _DUP:
            tokens.append(("P", _DUP[c]))
            i += 1
            continue
        elif c in "+-":
            sign = -1 if c == "-" else 1
            i += 1
        elif c.isdigit() or c == ".":
            sign = 1
        else:
            i += 1
            continue
        # consume the rest of the number
        while i < n and (line[i].isdigit() or line[i] == "."):
            digits += line[i]
            i += 1
        # AFFN exponent (only meaningful when not compressed)
        if not compressed and i < n and line[i] in "eE":
            exp = line[i]
            j = i + 1
            if j < n and line[j] in "+-":
                exp += line[j]
                j += 1
            while j < n and line[j].isdigit():
                exp += line[j]
                j += 1
            digits += exp
            i = j
        if digits in ("", ".", "+", "-"):
            continue
        try:
            tokens.append((mode, sign * float(digits)))
        except ValueError:
            continue
    return tokens


def _decode_xydata(data_lines: list[str], yfactor: float):
    """Decode an X++(Y..Y) block into a flat array of ordinate values."""
    compressed = any(ch in _COMPRESSED_CHARS for ln in data_lines for ch in ln)
    all_y: list[float] = []
    prev_line_was_dif = False
    for li, line in enumerate(data_lines):
        toks = _decode_line(line, compressed)
        if not toks:
            continue
        # first token is the abscissa checkpoint -> drop it
        y_toks = toks[1:]
        line_ys: list[float] = []
        y = 0.0
        prev_mode, prev_val = None, None
        for mode, val in y_toks:
            if mode == "P":
                for _ in range(int(val) - 1):
                    if prev_mode == "D":
                        y += prev_val
                    else:
                        y = prev_val
                    line_ys.append(y)
            elif mode == "D":
                y += val
                line_ys.append(y)
                prev_mode, prev_val = "D", val
            else:  # absolute
                y = val
                line_ys.append(y)
                prev_mode, prev_val = "A", val
        has_dif = any(m == "D" for m, _ in y_toks)
        # DIF line-boundary check: first ordinate repeats previous line's last.
        if li > 0 and prev_line_was_dif and all_y and line_ys:
            if abs(line_ys[0] - all_y[-1]) <= max(1.0, abs(all_y[-1]) * 1e-6):
                line_ys = line_ys[1:]
        all_y.extend(line_ys)
        prev_line_was_dif = has_dif
    return np.asarray(all_y, dtype=float) * yfactor


def _decode_pairs(data_lines: list[str], xfactor: float, yfactor: float):
    """Decode (XY..XY) / peak-table style pairs into x, y arrays."""
    nums: list[float] = []
    for line in data_lines:
        for tok in re.split(r"[\s,;]+", line.strip()):
            if not tok:
                continue
            try:
                nums.append(float(tok))
            except ValueError:
                pass
    arr = np.asarray(nums[: len(nums) // 2 * 2], dtype=float).reshape(-1, 2)
    return arr[:, 0] * xfactor, arr[:, 1] * yfactor


def _parse_block(block: str, path: str, idx: int) -> Spectrum | None:
    ldr = {}
    data_kind = None
    data_lines: list[str] = []
    capturing = None
    for raw in block.splitlines():
        if raw.strip().startswith("##"):
            m = re.match(r"##\s*([^=]+)=(.*)", raw.strip())
            if not m:
                continue
            key = m.group(1).strip().upper()
            val = m.group(2).strip()
            if key in ("XYDATA", "XYPOINTS", "PEAK TABLE", "PEAKTABLE"):
                data_kind = "xydata" if key == "XYDATA" else "pairs"
                capturing = data_kind
                continue
            capturing = None
            ldr[key] = val
        elif capturing:
            data_lines.append(raw)

    if not data_lines:
        return None

    def num(key, default):
        try:
            return float(ldr.get(key, default))
        except (TypeError, ValueError):
            return default

    xfactor = num("XFACTOR", 1.0)
    yfactor = num("YFACTOR", 1.0)
    title = ldr.get("TITLE", "") or f"{os.path.splitext(os.path.basename(path))[0]}"
    if idx > 0 and ldr.get("TITLE"):
        title = ldr["TITLE"]

    if data_kind == "xydata":
        y = _decode_xydata(data_lines, yfactor)
        firstx = num("FIRSTX", 0.0)
        lastx = num("LASTX", float(len(y) - 1))
        npts = int(num("NPOINTS", len(y)))
        npts = npts if npts > 0 else len(y)
        if len(y) != npts:
            npts = len(y)
        x = np.linspace(firstx, lastx, npts)
        if len(x) != len(y):
            x = np.linspace(firstx, lastx, len(y))
    else:
        x, y = _decode_pairs(data_lines, xfactor, yfactor)

    if y.size == 0:
        return None

    xunit = _XUNIT_MAP.get(ldr.get("XUNITS", "").upper().strip(), "pixel")
    yunit = _YUNIT_MAP.get(ldr.get("YUNITS", "").upper().strip(), "intensity")
    meta = {"source": path, "jcamp": {k: ldr.get(k) for k in
            ("DATA TYPE", "ORIGIN", "OWNER", "DATE", "INSTRUMENT") if k in ldr}}
    return Spectrum(x=x, y=y, name=title.strip() or f"block {idx + 1}",
                    x_unit=xunit, y_unit=yunit, meta=meta)


def load_jcamp(path: str) -> list[Spectrum]:
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    # Split a LINK file into blocks; each spectrum block starts at ##TITLE.
    parts = re.split(r"(?m)^##TITLE\s*=", text)
    blocks = []
    if len(parts) <= 1:
        blocks = [text]
    else:
        blocks = ["##TITLE=" + p for p in parts[1:]]
    spectra = []
    for i, blk in enumerate(blocks):
        try:
            s = _parse_block(blk, path, i)
        except Exception:
            s = None
        if s is not None:
            spectra.append(s)
    if not spectra:
        raise ValueError(f"No spectral data blocks found in {os.path.basename(path)}.")
    return spectra
