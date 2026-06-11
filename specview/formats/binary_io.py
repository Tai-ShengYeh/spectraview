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
