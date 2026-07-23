"""Live instrument backends.

Currently:
    * ``oceanoptics``  — Ocean Optics / Ocean Insight via python-seabreeze
                         (USB2000/USB4000/Flame/QE Pro/... )
    * ``godirect_dev`` — Vernier Go Direct (best-effort; see module docstring
                         for the SpectroVis Plus limitation + CSV watch-folder
                         workflow)
    * ``innospectra``  — InnoSpectra NIRScan via the official Python SDK,
                         driven through a 32-bit bridge subprocess
                         (``isc_bridge.py``)

All backends import their hardware library lazily, so this package is safe to
import when no instrument SDK is installed.
"""
from __future__ import annotations

from .base import DeviceError, DeviceInfo, DeviceNotFound, SpectrometerDevice

__all__ = ["DeviceError", "DeviceInfo", "DeviceNotFound", "SpectrometerDevice"]
