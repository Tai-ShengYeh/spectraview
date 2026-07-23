"""Common interface for live spectrometer acquisition.

Every backend (Ocean Optics / Vernier Go Direct / ...) exposes:

* ``list_devices()``      -> list of ``DeviceInfo``
* a device class implementing :class:`SpectrometerDevice`

The UI (``specview.ui.acquisition``) only talks to this interface, so new
instruments can be added by dropping another module into this package.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..spectrum import Spectrum


class DeviceError(RuntimeError):
    """Any instrument-related failure with a user-readable message."""


class DeviceNotFound(DeviceError):
    """No instrument answered on the bus."""


@dataclass
class DeviceInfo:
    """One entry in a device list (returned by ``list_devices``)."""
    backend: str           # e.g. "oceanoptics"
    ident: str             # backend-specific id (serial number / address)
    label: str             # human-readable, shown in combo boxes
    extra: dict = field(default_factory=dict)


class SpectrometerDevice:
    """Abstract live spectrometer.

    Lifecycle: ``open()`` ... ``read_spectrum()`` ... ``close()``.
    Implementations must be safe to ``close()`` twice.
    """

    backend = "abstract"

    # -- lifecycle ---------------------------------------------------------
    def open(self) -> None:
        raise NotImplementedError

    def close(self) -> None:
        raise NotImplementedError

    @property
    def is_open(self) -> bool:
        raise NotImplementedError

    # -- identity ----------------------------------------------------------
    @property
    def label(self) -> str:
        return "spectrometer"

    # -- acquisition parameters ---------------------------------------------
    @property
    def integration_limits_ms(self) -> tuple[float, float]:
        """(min, max) integration time in milliseconds."""
        return (1.0, 10_000.0)

    def set_integration_time_ms(self, ms: float) -> None:
        raise NotImplementedError

    # -- data ----------------------------------------------------------------
    def wavelengths(self) -> np.ndarray:
        """Wavelength axis in nm (constant per device)."""
        raise NotImplementedError

    def read_intensities(self) -> np.ndarray:
        """One freshly acquired intensity array (blocking)."""
        raise NotImplementedError

    def read_spectrum(self, name: str | None = None,
                      averages: int = 1) -> Spectrum:
        """Acquire ``averages`` scans and return their mean as a Spectrum."""
        if averages < 1:
            averages = 1
        acc = None
        for _ in range(averages):
            y = np.asarray(self.read_intensities(), dtype=float)
            acc = y if acc is None else acc + y
        y = acc / averages
        x = np.asarray(self.wavelengths(), dtype=float)
        return Spectrum(x, y, name=name or self.label,
                        x_unit="nm", y_unit="counts",
                        meta={"source": f"live:{self.backend}",
                              "device": self.label,
                              "averages": averages})

    # -- context manager sugar -------------------------------------------------
    def __enter__(self):
        self.open()
        return self

    def __exit__(self, *exc):
        self.close()
        return False
