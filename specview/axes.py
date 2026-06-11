"""Axis-type conversions for the x (spectral) and y (intensity) axes.

X conversions route through vacuum wavelength in nanometres as the canonical
representation. Y conversions route through linear transmittance/reflectance
fraction (0..1). Conversions return a NEW Spectrum (data re-sorted ascending).
"""
from __future__ import annotations

import numpy as np

from .spectrum import Spectrum

# Physical constants for the optical conversions.
_HC_EV_NM = 1239.841984  # h*c in eV·nm  -> E[eV] = 1239.84 / λ[nm]
_C_THZ_NM = 299792.458    # c in nm·THz   -> f[THz] = 299792.458 / λ[nm]

# Units that describe an optical x-axis (convertible to wavelength).
X_CONVERTIBLE = {"nm", "um", "cm-1", "raman_cm-1", "eV", "THz"}
# Units that describe an optical y-axis (convertible via T/R fraction).
Y_CONVERTIBLE = {"transmittance", "%T", "absorbance",
                 "reflectance", "%R", "KM", "log1R"}

_EPS = 1e-12


# ======================================================================== X
def x_to_nm(x: np.ndarray, unit: str, laser_nm: float | None = None) -> np.ndarray:
    """Convert an x array in ``unit`` to wavelength in nm."""
    x = np.asarray(x, dtype=float)
    if unit == "nm":
        return x.copy()
    if unit == "um":
        return x * 1000.0
    if unit == "cm-1":
        return 1.0e7 / np.where(np.abs(x) < _EPS, np.nan, x)
    if unit == "eV":
        return _HC_EV_NM / np.where(np.abs(x) < _EPS, np.nan, x)
    if unit == "THz":
        return _C_THZ_NM / np.where(np.abs(x) < _EPS, np.nan, x)
    if unit == "raman_cm-1":
        if not laser_nm:
            raise ValueError("Raman shift conversion needs the laser wavelength (nm).")
        wn_laser = 1.0e7 / laser_nm           # cm-1
        wn_abs = wn_laser - x                  # absolute wavenumber of scattered light
        return 1.0e7 / np.where(np.abs(wn_abs) < _EPS, np.nan, wn_abs)
    raise ValueError(f"Cannot convert x-unit '{unit}' to wavelength.")


def nm_to_x(nm: np.ndarray, unit: str, laser_nm: float | None = None) -> np.ndarray:
    """Convert wavelength in nm to an x array in ``unit``."""
    nm = np.asarray(nm, dtype=float)
    safe = np.where(np.abs(nm) < _EPS, np.nan, nm)
    if unit == "nm":
        return nm.copy()
    if unit == "um":
        return nm / 1000.0
    if unit == "cm-1":
        return 1.0e7 / safe
    if unit == "eV":
        return _HC_EV_NM / safe
    if unit == "THz":
        return _C_THZ_NM / safe
    if unit == "raman_cm-1":
        if not laser_nm:
            raise ValueError("Raman shift conversion needs the laser wavelength (nm).")
        wn_laser = 1.0e7 / laser_nm
        return wn_laser - (1.0e7 / safe)
    raise ValueError(f"Cannot convert wavelength to x-unit '{unit}'.")


def convert_x(spec: Spectrum, target: str, laser_nm: float | None = None) -> Spectrum:
    """Return a copy of ``spec`` with the x-axis expressed in ``target`` unit."""
    if spec.x_unit == target:
        return spec.copy()
    if spec.x_unit not in X_CONVERTIBLE or target not in X_CONVERTIBLE:
        raise ValueError(
            f"x-axis conversion {spec.x_unit!r} → {target!r} is not supported "
            "(pixel axes need a calibration first)."
        )
    laser = laser_nm or spec.meta.get("laser_nm")
    nm = x_to_nm(spec.x, spec.x_unit, laser)
    new_x = nm_to_x(nm, target, laser)
    out = spec.replace_data(new_x, spec.y)  # replace_data re-sorts ascending
    out.x_unit = target
    if laser:
        out.meta["laser_nm"] = laser
    return out


# ======================================================================== Y
def y_to_fraction(y: np.ndarray, unit: str) -> np.ndarray:
    """Convert a y array in ``unit`` to linear transmittance/reflectance (0..1)."""
    y = np.asarray(y, dtype=float)
    if unit in ("transmittance", "reflectance"):
        return y.copy()
    if unit in ("%T", "%R"):
        return y / 100.0
    if unit == "absorbance":
        return np.power(10.0, -y)
    if unit == "log1R":
        return np.power(10.0, -y)
    if unit == "KM":
        # invert KM = (1-R)^2 / (2R)  ->  R = (1+KM) - sqrt((1+KM)^2 - 1)
        k = np.maximum(y, 0.0)
        return (1.0 + k) - np.sqrt(np.maximum((1.0 + k) ** 2 - 1.0, 0.0))
    raise ValueError(f"Cannot convert y-unit '{unit}' to a T/R fraction.")


def fraction_to_y(frac: np.ndarray, unit: str) -> np.ndarray:
    """Convert a linear T/R fraction (0..1) to a y array in ``unit``."""
    frac = np.asarray(frac, dtype=float)
    if unit in ("transmittance", "reflectance"):
        return frac.copy()
    if unit in ("%T", "%R"):
        return frac * 100.0
    if unit == "absorbance":
        return -np.log10(np.clip(frac, _EPS, None))
    if unit == "log1R":
        return -np.log10(np.clip(frac, _EPS, None))
    if unit == "KM":
        r = np.clip(frac, _EPS, None)
        return (1.0 - r) ** 2 / (2.0 * r)
    raise ValueError(f"Cannot convert a T/R fraction to y-unit '{unit}'.")


def convert_y(spec: Spectrum, target: str) -> Spectrum:
    """Return a copy of ``spec`` with the y-axis expressed in ``target`` unit."""
    if spec.y_unit == target:
        return spec.copy()
    if spec.y_unit not in Y_CONVERTIBLE or target not in Y_CONVERTIBLE:
        raise ValueError(
            f"y-axis conversion {spec.y_unit!r} → {target!r} is not supported "
            "(raw intensity/counts cannot be turned into absorbance)."
        )
    frac = y_to_fraction(spec.y, spec.y_unit)
    new_y = fraction_to_y(frac, target)
    out = spec.copy()
    out.y = new_y
    out.y_unit = target
    return out
