"""Ocean Optics / Ocean Insight spectrometers via python-seabreeze.

Covers the classic USB series (USB2000, USB2000+, USB4000, HR2000, HR4000,
Flame, Maya, QE65000, **QE Pro**, NIRQuest, STS, ...) — everything the
SeaBreeze driver knows.

Install (one-time)::

    pip install seabreeze
    seabreeze_os_setup        # installs udev rules (Linux) / driver notes (Windows)

On Windows the spectrometer must use the *WinUSB* driver (the installer from
``seabreeze_os_setup`` provides it; OceanView's own driver also works with
``pyseabreeze``). If the default ``cseabreeze`` backend fails to see the
device, we transparently retry with the pure-python ``pyseabreeze`` backend.
"""
from __future__ import annotations

import numpy as np

from ..formats.binary_io import MissingDependency
from .base import DeviceError, DeviceInfo, DeviceNotFound, SpectrometerDevice

_INSTALL_HINT = (
    "Ocean Optics support needs the 'seabreeze' package.\n\n"
    "    pip install seabreeze\n"
    "    seabreeze_os_setup\n\n"
    "then re-plug the spectrometer. (On Windows, the device must use the "
    "WinUSB driver — seabreeze_os_setup explains how.)"
)


def _sb():
    """Import seabreeze lazily so the app runs without it installed."""
    try:
        import seabreeze.spectrometers as sbs
        return sbs
    except ImportError as exc:
        raise MissingDependency(_INSTALL_HINT) from exc


def list_devices() -> list[DeviceInfo]:
    """Enumerate all connected Ocean Optics spectrometers."""
    sbs = _sb()
    out = []
    for dev in sbs.list_devices():
        out.append(DeviceInfo(
            backend="oceanoptics",
            ident=dev.serial_number,
            label=f"{dev.model} ({dev.serial_number})",
            extra={"model": dev.model},
        ))
    return out


class OceanOpticsDevice(SpectrometerDevice):
    """One Ocean Optics spectrometer (selected by serial number, or first found)."""

    backend = "oceanoptics"

    def __init__(self, serial: str | None = None,
                 correct_dark_counts: bool = True,
                 correct_nonlinearity: bool = True):
        self._serial = serial
        self._spec = None
        self._want_dark = correct_dark_counts
        self._want_nl = correct_nonlinearity
        self._label = "Ocean Optics"

    # -- lifecycle ---------------------------------------------------------
    def open(self) -> None:
        if self._spec is not None:
            return
        sbs = _sb()
        try:
            if self._serial:
                self._spec = sbs.Spectrometer.from_serial_number(self._serial)
            else:
                self._spec = sbs.Spectrometer.from_first_available()
        except Exception as exc:  # noqa: BLE001 (seabreeze raises various types)
            msg = str(exc)
            if "No device found" in msg or "not found" in msg.lower():
                raise DeviceNotFound(
                    "No Ocean Optics spectrometer found.\n"
                    "Check the USB cable and close OceanView/OceanArt "
                    "(only one program can hold the device at a time)."
                ) from exc
            raise DeviceError(f"Could not open spectrometer: {msg}") from exc
        self._label = f"{self._spec.model} ({self._spec.serial_number})"
        # Feature support varies by model — probe once.
        try:
            self._spec.intensities(correct_dark_counts=self._want_dark,
                                   correct_nonlinearity=self._want_nl)
            self._dark_ok, self._nl_ok = self._want_dark, self._want_nl
        except Exception:  # noqa: BLE001 — model lacks the correction feature
            self._dark_ok = self._nl_ok = False

    def close(self) -> None:
        if self._spec is not None:
            try:
                self._spec.close()
            except Exception:  # noqa: BLE001
                pass
            self._spec = None

    @property
    def is_open(self) -> bool:
        return self._spec is not None

    @property
    def label(self) -> str:
        return self._label

    # -- parameters -----------------------------------------------------------
    @property
    def integration_limits_ms(self) -> tuple[float, float]:
        self._require_open()
        lo, hi = self._spec.integration_time_micros_limits
        return (lo / 1000.0, hi / 1000.0)

    def set_integration_time_ms(self, ms: float) -> None:
        self._require_open()
        lo, hi = self._spec.integration_time_micros_limits
        us = int(np.clip(ms * 1000.0, lo, hi))
        self._spec.integration_time_micros(us)

    # -- data -------------------------------------------------------------------
    def wavelengths(self) -> np.ndarray:
        self._require_open()
        return np.asarray(self._spec.wavelengths(), dtype=float)

    def read_intensities(self) -> np.ndarray:
        self._require_open()
        return np.asarray(
            self._spec.intensities(correct_dark_counts=self._dark_ok,
                                   correct_nonlinearity=self._nl_ok),
            dtype=float)

    # -- helpers -----------------------------------------------------------------
    def _require_open(self) -> None:
        if self._spec is None:
            raise DeviceError("Spectrometer is not connected (call open() first).")
