"""Readers for binary instrument formats (optional dependencies).

These degrade gracefully: if the needed package is not installed, a clear
ImportError with a ``pip install`` hint is raised instead of crashing.
"""
from __future__ import annotations

import os

import numpy as np

from ..spectrum import Spectrum


class MissingDependency(ImportError):
    """Raised when an optional reader's backend package is not installed."""


def load_spc(path: str) -> list[Spectrum]:
    """GRAMS/Galactic .spc — needs ``pip install spc-spectra``."""
    try:
        import spc_spectra as spc  # type: ignore
    except ImportError:
        try:
            import spc  # type: ignore  # older package name
        except ImportError as exc:
            raise MissingDependency(
                "Reading .spc files needs the 'spc-spectra' package.\n"
                "Install it with:  pip install spc-spectra"
            ) from exc

    f = spc.File(path)
    base = os.path.splitext(os.path.basename(path))[0]
    # Map SPC unit codes to our axis units (best-effort).
    xlabel = (getattr(f, "xlabel", "") or "").lower()
    ylabel = (getattr(f, "ylabel", "") or "").lower()
    x_unit = ("cm-1" if "cm-1" in xlabel or "wavenumber" in xlabel else
              "nm" if "nm" in xlabel or "nanomet" in xlabel else
              "raman_cm-1" if "raman" in xlabel else "pixel")
    y_unit = ("absorbance" if "absorb" in ylabel else
              "transmittance" if "trans" in ylabel else
              "counts" if "count" in ylabel else "intensity")

    spectra = []
    for i, sub in enumerate(getattr(f, "sub", [])):
        x = np.asarray(getattr(f, "x", getattr(sub, "x", [])), dtype=float)
        y = np.asarray(sub.y, dtype=float)
        if x.size != y.size:
            x = np.arange(y.size, dtype=float)
        name = base if len(f.sub) == 1 else f"{base} [{i + 1}]"
        spectra.append(Spectrum(x=x, y=y, name=name, x_unit=x_unit, y_unit=y_unit,
                                meta={"source": path}))
    if not spectra:
        raise ValueError(f"No subfiles found in {os.path.basename(path)}.")
    return spectra


def load_opus(path: str) -> list[Spectrum]:
    """Bruker OPUS (.0, .1, ...) — needs ``pip install brukeropusreader``."""
    try:
        from brukeropusreader import read_file  # type: ignore
    except ImportError as exc:
        raise MissingDependency(
            "Reading Bruker OPUS files needs the 'brukeropusreader' package.\n"
            "Install it with:  pip install brukeropusreader"
        ) from exc

    data = read_file(path)
    base = os.path.basename(path)
    spectra = []
    # OPUS files contain several blocks; pick the spectral ones we recognise.
    for block in ("AB", "Absorbance", "Transmittance", "ScSm", "ScRf", "Raman"):
        if block not in data:
            continue
        y = np.asarray(data[block], dtype=float)
        try:
            x = np.asarray(data.get_range(block), dtype=float)
        except Exception:
            params = data.get(block + " Data Parameter", {})
            fx, lx = params.get("FXV"), params.get("LXV")
            x = (np.linspace(fx, lx, y.size) if fx is not None and lx is not None
                 else np.arange(y.size, dtype=float))
        n = min(x.size, y.size)
        y_unit = ("absorbance" if block in ("AB", "Absorbance") else
                  "transmittance" if block == "Transmittance" else "intensity")
        spectra.append(Spectrum(x=x[:n], y=y[:n], name=f"{base} ({block})",
                                x_unit="cm-1", y_unit=y_unit, meta={"source": path}))
    if not spectra:
        raise ValueError(f"No recognised spectral block in {base}.")
    return spectra


# ---------------------------------------------------------------------------
# PerkinElmer Spectrum (.sp)  --  pure-Python, no external dependency.
# ---------------------------------------------------------------------------
# The ".sp" file is a tagged-block binary ("PEPE2D constant interval DataSet
# file"). After a NUL-terminated text signature comes a sequence of blocks,
# each:  <block-id : uint16 LE>  <length : int32 LE>  <content[length]>.
# Every real block-id has high byte 0x8b. The blocks we need:
#   0x8b72  X range   -> member tag (2 B) + 2 float64  (firstX, lastX)
#   0x8b74  delta X   -> member tag (2 B) + 1 float64
#   0x8b75  N points  -> member tag (2 B) + 1 int32
#   0x8b77  X units   -> member tag (2 B) + uint16 len + ASCII
#   0x8b78  Y units   -> member tag (2 B) + uint16 len + ASCII
#   0x8b7b  data type -> member tag (2 B) + uint16 len + ASCII ("Spectrum")
#   0x8b7c  Y data    -> member tag (2 B) + int32 byte-count + float64[N]
import struct

_PE_SIG = b"PEPE"

# PerkinElmer unit strings -> SpectraView canonical units.
_PE_XUNITS = {
    "cm-1": "cm-1", "1/cm": "cm-1", "cm^-1": "cm-1", "wavenumber": "cm-1",
    "nm": "nm", "nanometers": "nm", "wavelength": "nm",
    "um": "um", "µm": "um", "micrometers": "um", "micron": "um",
    "ev": "eV", "raman": "raman_cm-1",
}
_PE_YUNITS = {
    "%t": "%T", "%transmittance": "%T", "transmittance": "transmittance",
    "a": "absorbance", "abs": "absorbance", "absorbance": "absorbance",
    "%r": "%R", "reflectance": "reflectance", "log(1/r)": "log1R",
    "kubelka-munk": "KM", "km": "KM", "counts": "counts",
}


def _pe_blocks(raw: bytes) -> dict:
    """Return {block_id: (content_offset, length)} for a PerkinElmer .sp file."""
    start = raw.find(b"\x00")
    i = (start + 1) if start != -1 else 4
    blocks: dict[int, tuple[int, int]] = {}
    n = len(raw)
    while i < n - 6:
        if raw[i + 1] != 0x8B:          # second byte of every block-id is 0x8b
            i += 1
            continue
        bid = struct.unpack_from("<H", raw, i)[0]
        length = struct.unpack_from("<i", raw, i + 2)[0]
        if length < 0 or i + 6 + length > n:
            i += 1                       # spurious 0x8b inside data: resync
            continue
        blocks.setdefault(bid, (i + 6, length))
        i += 6 + length
    return blocks


def load_sp(path: str) -> list[Spectrum]:
    """PerkinElmer Spectrum ``.sp`` (FTIR/UV-Vis). No extra package required."""
    with open(path, "rb") as fh:
        raw = fh.read()
    if not raw.startswith(_PE_SIG):
        raise ValueError(
            f"{os.path.basename(path)} is not a PerkinElmer .sp file "
            "(missing 'PEPE' signature)."
        )
    blocks = _pe_blocks(raw)

    def _doubles(bid, count):
        if bid not in blocks:
            return None
        off, ln = blocks[bid]
        off += 2                         # skip the 2-byte member tag
        if off + 8 * count > len(raw):
            return None
        return struct.unpack_from("<%dd" % count, raw, off)

    def _string(bid):
        if bid not in blocks:
            return ""
        off, ln = blocks[bid]
        slen = struct.unpack_from("<H", raw, off + 2)[0]   # tag(2) + len(2)
        return raw[off + 4:off + 4 + slen].decode("latin1").strip()

    # --- Y data (block 0x8b7c): member tag (2 B) + int32 byte-count + doubles
    if 0x8B7C not in blocks:
        raise ValueError(f"No data block found in {os.path.basename(path)}.")
    doff, _ = blocks[0x8B7C]
    nbytes = struct.unpack_from("<i", raw, doff + 2)[0]
    npts = nbytes // 8
    y = np.frombuffer(raw, "<f8", npts, doff + 6).astype(float)

    # --- X axis from (firstX, lastX, npoints); fall back to delta if needed
    xr = _doubles(0x8B72, 2)
    np_block = _doubles(0x8B75, 1)
    n_meta = int(struct.unpack_from("<i", raw, blocks[0x8B75][0] + 2)[0]) \
        if 0x8B75 in blocks else npts
    if xr is not None:
        first_x, last_x = xr
        x = np.linspace(first_x, last_x, npts)
    else:                                # reconstruct from firstX + deltaX
        dx = _doubles(0x8B74, 1)
        first_x = (xr[0] if xr else 0.0)
        step = dx[0] if dx else 1.0
        x = first_x + step * np.arange(npts)

    x_unit = _PE_XUNITS.get(_string(0x8B77).lower(), "cm-1")
    y_unit = _PE_YUNITS.get(_string(0x8B78).lower(), "transmittance")

    base = os.path.splitext(os.path.basename(path))[0]
    return [Spectrum(x=x, y=y, name=base, x_unit=x_unit, y_unit=y_unit,
                     meta={"source": path, "format": "PerkinElmer .sp",
                           "npoints_meta": n_meta})]
