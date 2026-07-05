"""Webcam/grating spectrometer: a spectrum photo -> calibrated intensity curve.

Implements the core of a Theremino-style spectrometer, working on a captured
image (a photo of a diffraction-grating spectrum) rather than a live camera:

  1. take a horizontal strip (ROI) of the image and average its rows,
  2. read an intensity profile along the columns (luminance or a colour
     channel),
  3. map pixel column -> wavelength with a linear or quadratic calibration
     fitted from known reference lines (e.g. a fluorescent lamp's mercury peaks
     at 436, 546, 611 nm).

Pure numpy (+ Pillow only to read image files), so it is unit-tested headlessly.
"""
from __future__ import annotations

import io
import os

import numpy as np

from .core import make_spectrum

# ITU-R 601 luma weights, matching orange-assay's grayscale convention.
_LUMA = np.array([0.299, 0.587, 0.114])
CHANNELS = ["luminance", "red", "green", "blue", "sum(RGB)"]
CAL_MODELS = [("linear", 1), ("quadratic", 2)]

# A common calibration reference: fluorescent-lamp mercury/terbium lines (nm).
FLUORESCENT_LINES = [405.4, 435.8, 546.1, 611.6]


def load_rgb(source) -> np.ndarray:
    """Return an (H, W, 3) float RGB array in 0..255 from a path/bytes/array."""
    if isinstance(source, (str, bytes)):
        from PIL import Image
        buf = io.BytesIO(source) if isinstance(source, bytes) else source
        arr = np.asarray(Image.open(buf).convert("RGB"), dtype=float)
    else:
        arr = np.asarray(source, dtype=float)
        if arr.ndim == 2:
            arr = np.repeat(arr[:, :, None], 3, axis=2)
    if arr.max() <= 1.0 + 1e-9:
        arr = arr * 255.0
    return arr[..., :3]


def channel_image(rgb: np.ndarray, channel: str) -> np.ndarray:
    """Reduce an RGB image to a single-channel (H, W) intensity image."""
    if channel == "red":
        return rgb[..., 0]
    if channel == "green":
        return rgb[..., 1]
    if channel == "blue":
        return rgb[..., 2]
    if channel == "sum(RGB)":
        return rgb.sum(axis=2)
    return rgb @ _LUMA                       # luminance (default)


def extract_profile(rgb: np.ndarray, channel: str = "luminance",
                    row_center: float | None = None,
                    row_frac: float = 0.15) -> np.ndarray:
    """Average a horizontal ROI band into a 1-D intensity profile along columns.

    ``row_center`` is the ROI centre row (defaults to the image middle);
    ``row_frac`` is the band height as a fraction of image height (min 1 row).
    """
    img = channel_image(np.asarray(rgb, float), channel)
    h = img.shape[0]
    cy = h / 2.0 if row_center is None else float(row_center)
    half = max(0.5, row_frac * h / 2.0)
    lo = max(0, int(round(cy - half)))
    hi = min(h, int(round(cy + half)) + 1)
    if hi <= lo:
        lo, hi = 0, h
    return img[lo:hi, :].mean(axis=0)


def fit_calibration(pixels, wavelengths, model: str = "linear"):
    """Fit pixel -> wavelength as a polynomial. Returns (coeffs, r2, degree).

    ``coeffs`` are numpy.polyfit order (highest power first). Needs degree+1
    points; falls back to a lower degree when too few are supplied.
    """
    px = np.asarray(pixels, float)
    wl = np.asarray(wavelengths, float)
    if px.size != wl.size or px.size < 2:
        raise ValueError("Calibration needs at least 2 (pixel, wavelength) points.")
    degree = dict(CAL_MODELS).get(model, 1)
    degree = min(degree, px.size - 1)
    coeffs = np.polyfit(px, wl, degree)
    pred = np.polyval(coeffs, px)
    ss_res = float(np.sum((wl - pred) ** 2))
    ss_tot = float(np.sum((wl - wl.mean()) ** 2)) or 1.0
    return coeffs, 1.0 - ss_res / ss_tot, degree


def image_to_spectrum(source, channel: str = "luminance",
                      row_center: float | None = None, row_frac: float = 0.15,
                      calibration=None, model: str = "linear",
                      name: str = "", flip: bool = False) -> dict:
    """Convert a spectrum photo into a spectrum dict.

    ``calibration`` is an optional list of (pixel, wavelength) pairs; with it,
    the x-axis is wavelength (nm), otherwise it is the pixel column index.
    ``flip`` reverses the profile (when the grating spreads blue-to-red the
    other way). Returns make_spectrum(...) with an extra "calibration" entry.
    """
    rgb = load_rgb(source)
    profile = extract_profile(rgb, channel=channel, row_center=row_center,
                              row_frac=row_frac)
    if flip:
        profile = profile[::-1]
    px = np.arange(profile.size, dtype=float)

    x, x_label, cal = px, "pixel", None
    if calibration:
        pts = [(float(p), float(w)) for p, w in calibration]
        coeffs, r2, degree = fit_calibration([p for p, _ in pts],
                                             [w for _, w in pts], model)
        x = np.polyval(coeffs, px)
        x_label = "wavelength (nm)"
        cal = {"coeffs": coeffs, "r2": r2, "degree": degree, "points": pts}
        if np.all(np.diff(x) < 0):           # keep ascending wavelength
            x, profile = x[::-1], profile[::-1]

    if not name:
        name = (os.path.splitext(os.path.basename(source))[0]
                if isinstance(source, str) else "spectrum")
    spec = make_spectrum(x, profile, name=name, x_label=x_label,
                         source=source if isinstance(source, str) else "")
    spec["calibration"] = cal
    return spec
