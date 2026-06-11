"""Spectrum data model and the document that holds a set of spectra."""
from __future__ import annotations

import copy
import itertools
from dataclasses import dataclass, field

import numpy as np

# A pleasant, high-contrast qualitative palette (Tableau-ish), cycled for new spectra.
PALETTE = [
    "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
    "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f",
    "#393b79", "#b5651d", "#637939", "#8c6d31", "#843c39",
]

# Canonical labels for axis units -> (display text). Conversions live in axes.py.
X_UNIT_LABELS = {
    "pixel": "Pixel",
    "nm": "Wavelength (nm)",
    "um": "Wavelength (µm)",
    "cm-1": "Wavenumber (cm⁻¹)",
    "raman_cm-1": "Raman shift (cm⁻¹)",
    "eV": "Energy (eV)",
    "THz": "Frequency (THz)",
}

Y_UNIT_LABELS = {
    "intensity": "Intensity",
    "counts": "Counts",
    "a.u.": "Intensity (a.u.)",
    "transmittance": "Transmittance",
    "%T": "Transmittance (%)",
    "absorbance": "Absorbance",
    "reflectance": "Reflectance",
    "%R": "Reflectance (%)",
    "KM": "Kubelka-Munk",
    "log1R": "log(1/R)",
}


# eq=False is REQUIRED: a Spectrum holds numpy arrays, and the dataclass-generated
# __eq__ would compare them field-by-field. Membership tests (``spec in set``,
# ``list.index``) used by SpectrumSet would then trigger array comparison between
# spectra of different lengths -> "operands could not be broadcast" ValueError.
# We track spectra by identity, so default object identity equality is exactly right.
@dataclass(eq=False)
class Spectrum:
    """A single spectrum: paired x / y arrays plus presentation metadata.

    Internally ``x`` is always stored ascending (sorted on construction). The
    display can be flipped independently (e.g. IR is shown high→low wavenumber).
    """

    x: np.ndarray
    y: np.ndarray
    name: str = "spectrum"
    x_unit: str = "pixel"
    y_unit: str = "intensity"
    color: str | None = None
    visible: bool = True
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.x = np.asarray(self.x, dtype=float).ravel()
        self.y = np.asarray(self.y, dtype=float).ravel()
        if self.x.size != self.y.size:
            n = min(self.x.size, self.y.size)
            self.x, self.y = self.x[:n], self.y[:n]
        # Drop Na/inf samples, then sort ascending by x and de-duplicate x.
        good = np.isfinite(self.x) & np.isfinite(self.y)
        self.x, self.y = self.x[good], self.y[good]
        if self.x.size:
            order = np.argsort(self.x, kind="stable")
            self.x, self.y = self.x[order], self.y[order]
            uniq = np.concatenate(([True], np.diff(self.x) != 0))
            self.x, self.y = self.x[uniq], self.y[uniq]

    # ---- basic helpers ---------------------------------------------------
    @property
    def npoints(self) -> int:
        return int(self.x.size)

    @property
    def x_label(self) -> str:
        return X_UNIT_LABELS.get(self.x_unit, self.x_unit)

    @property
    def y_label(self) -> str:
        return Y_UNIT_LABELS.get(self.y_unit, self.y_unit)

    @property
    def xrange(self) -> tuple[float, float]:
        return (float(self.x.min()), float(self.x.max())) if self.x.size else (0.0, 1.0)

    @property
    def yrange(self) -> tuple[float, float]:
        return (float(self.y.min()), float(self.y.max())) if self.y.size else (0.0, 1.0)

    def copy(self) -> "Spectrum":
        return copy.deepcopy(self)

    def replace_data(self, x: np.ndarray, y: np.ndarray) -> "Spectrum":
        """Return a copy with new data but identical presentation metadata."""
        s = self.copy()
        s.x = np.asarray(x, dtype=float).ravel()
        s.y = np.asarray(y, dtype=float).ravel()
        s.__post_init__()
        return s

    def value_at(self, x0: float) -> float:
        """Linearly interpolated y at x0 (edge-clamped)."""
        if self.x.size == 0:
            return float("nan")
        return float(np.interp(x0, self.x, self.y))

    def uniform_step(self) -> float:
        """Median sample spacing (used by filters that assume uniform x)."""
        if self.x.size < 2:
            return 1.0
        return float(np.median(np.diff(self.x)))

    def is_uniform(self, rtol: float = 1e-3) -> bool:
        if self.x.size < 3:
            return True
        d = np.diff(self.x)
        return bool(np.allclose(d, d[0], rtol=rtol))


class SpectrumSet:
    """The open document: an ordered list of spectra + colour assignment."""

    def __init__(self) -> None:
        self.spectra: list[Spectrum] = []
        self._color_cycle = itertools.cycle(PALETTE)

    # ---- container protocol ---------------------------------------------
    def __len__(self) -> int:
        return len(self.spectra)

    def __iter__(self):
        return iter(self.spectra)

    def __getitem__(self, i: int) -> Spectrum:
        return self.spectra[i]

    # ---- mutation --------------------------------------------------------
    def add(self, spec: Spectrum, at: int | None = None) -> Spectrum:
        if spec.color is None:
            spec.color = next(self._color_cycle)
        if at is None:
            self.spectra.append(spec)
        else:
            self.spectra.insert(at, spec)
        return spec

    def add_many(self, specs) -> list[Spectrum]:
        return [self.add(s) for s in specs]

    def remove(self, spec: Spectrum) -> None:
        if spec in self.spectra:
            self.spectra.remove(spec)

    def remove_at(self, i: int) -> None:
        if 0 <= i < len(self.spectra):
            del self.spectra[i]

    def clear(self) -> None:
        self.spectra.clear()

    def replace(self, old: Spectrum, new: Spectrum) -> None:
        """Swap a spectrum in place, keeping its colour and list position."""
        if old in self.spectra:
            i = self.spectra.index(old)
            if new.color is None:
                new.color = old.color
            self.spectra[i] = new

    def index(self, spec: Spectrum) -> int:
        return self.spectra.index(spec)

    def visible_spectra(self) -> list[Spectrum]:
        return [s for s in self.spectra if s.visible]

    def snapshot(self) -> list[Spectrum]:
        """Deep copy of all spectra, for the undo stack."""
        return copy.deepcopy(self.spectra)

    def restore(self, snap: list[Spectrum]) -> None:
        self.spectra = copy.deepcopy(snap)
