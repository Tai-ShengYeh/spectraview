"""Shared heuristics for structured readers (JSON, MATLAB): key aliases that
name the x / y arrays and units, plus unit-string normalisation.
"""
from __future__ import annotations

X_KEYS = ["x", "wavelength", "wavelengths", "wavenumber", "wavenumbers",
          "raman_shift", "ramanshift", "shift", "energy", "freq", "frequency",
          "x_axis", "xaxis", "xdata", "wn", "nm"]
Y_KEYS = ["y", "intensity", "intensities", "absorbance", "abs", "transmittance",
          "reflectance", "counts", "signal", "spectrum", "y_axis", "yaxis",
          "ydata", "value", "values", "data_y"]
NAME_KEYS = ["name", "title", "label", "id", "sample"]
XUNIT_KEYS = ["x_unit", "xunit", "x_units", "xunits", "xlabel", "x_label"]
YUNIT_KEYS = ["y_unit", "yunit", "y_units", "yunits", "ylabel", "y_label"]
LASER_KEYS = ["laser", "laser_nm", "excitation", "excitation_nm",
              "laser_wavelength", "excitation_wavelength"]


def match_key(keys_present, aliases):
    """Return the actual key among ``keys_present`` matching any alias (case-insensitive)."""
    lower = {str(k).lower(): k for k in keys_present}
    for a in aliases:
        if a.lower() in lower:
            return lower[a.lower()]
    return None


def norm_xunit(s):
    if not s:
        return None
    t = str(s).lower().strip()
    if "raman" in t:
        return "raman_cm-1"
    if "cm-1" in t or "cm^-1" in t or "1/cm" in t or "wavenumber" in t:
        return "cm-1"
    if t in ("um", "µm", "micrometer", "micrometers", "micron", "microns"):
        return "um"
    if "nm" in t or "nanometer" in t or "wavelength" in t:
        return "nm"
    if t in ("ev", "electronvolt", "electron volt"):
        return "eV"
    if "thz" in t:
        return "THz"
    if t in ("pixel", "px", "point", "points", "index"):
        return "pixel"
    return None


def norm_yunit(s):
    if not s:
        return None
    t = str(s).lower().strip()
    if "absorb" in t or t == "abs":
        return "absorbance"
    if "%t" in t or t == "%transmittance":
        return "%T"
    if "transmit" in t:
        return "transmittance"
    if "%r" in t or "%reflect" in t:
        return "%R"
    if "reflect" in t:
        return "reflectance"
    if "kubelka" in t:
        return "KM"
    if "log(1/r)" in t or "log1r" in t:
        return "log1R"
    if "count" in t:
        return "counts"
    if "intensity" in t or t == "a.u." or "arb" in t:
        return "intensity"
    return None
