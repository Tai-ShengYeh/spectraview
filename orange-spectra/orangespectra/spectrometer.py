"""Webcam/grating spectrometer: a spectrum photo -> calibrated intensity curve.

Implements the core of a Theremino-style spectrometer, working on a captured
image (a photo of a diffraction-grating spectrum) rather than a live camera:

1. optionally rotate the image so the dispersion axis runs left-to-right,
2. take a horizontal strip (ROI) of the image and average its rows,
3. read an intensity profile along the columns (luminance or a colour
   channel),
4. map pixel column -> wavelength with a linear or quadratic calibration
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

# Rotation applied *before* the ROI strip is taken. Labels are what the widget
# shows; the value is the clockwise rotation in degrees.
ROTATIONS = [("0deg (no rotation)", 0),
             ("90deg clockwise", 90),
             ("180deg", 180),
             ("270deg clockwise (= 90deg counter-clockwise)", 270)]

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


def rotate_rgb(rgb: np.ndarray, rotate: int = 0) -> np.ndarray:
    """Rotate an (H, W, 3) image clockwise by 0/90/180/270 degrees.

    Photos taken through a hand spectroscope often have the dispersion axis
    running vertically; rotating by 90 degrees puts it along the columns so the
    horizontal-strip machinery below applies unchanged.
    """
    arr = np.asarray(rgb)
    deg = int(rotate) % 360
    if deg == 0:
        return arr
    if deg not in (90, 180, 270):
        raise ValueError("rotate must be one of 0, 90, 180, 270 (degrees).")
    # np.rot90 turns counter-clockwise, so clockwise by d == CCW by 360-d.
    return np.ascontiguousarray(np.rot90(arr, k=(360 - deg) // 90))


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
    return rgb @ _LUMA  # luminance (default)


def extract_profile(rgb: np.ndarray, channel: str = "luminance",
                    row_center: float | None = None,
                    row_frac: float = 0.15,
                    rotate: int = 0) -> np.ndarray:
    """Average a horizontal ROI band into a 1-D intensity profile along columns.

    ``row_center`` is the ROI centre row (defaults to the image middle);
    ``row_frac`` is the band height as a fraction of image height (min 1 row).
    ``rotate`` (0/90/180/270, clockwise) is applied before the ROI is taken.
    """
    img = channel_image(rotate_rgb(np.asarray(rgb, float), rotate), channel)
    h = img.shape[0]
    cy = h / 2.0 if row_center is None else float(row_center)
    half = max(0.5, row_frac * h / 2.0)
    lo = max(0, int(round(cy - half)))
    hi = min(h, int(round(cy + half)) + 1)
    if hi <= lo:
        lo, hi = 0, h
    return img[lo:hi, :].mean(axis=0)


def brightest_row(rgb: np.ndarray, channel: str = "sum(RGB)",
                  rotate: int = 0, smooth: int = 5) -> int:
    """Row index (in the rotated image) where the spectrum band lies.

    Each row is scored by its 99th-percentile *total* (summed RGB) intensity,
    i.e. by the brightest features it contains rather than by its mean.
    Scoring the row mean lets a large, dim lens flare beat a thin, bright
    spectrum band -- and a single colour channel makes that worse, because
    the band's red or blue part is short. ``channel`` is accepted for API
    symmetry with :func:`extract_profile` but has no effect on the choice of
    row. ``smooth`` (rows) suppresses single hot rows.
    """
    del channel  # see docstring: the band is located on total intensity
    img = channel_image(rotate_rgb(np.asarray(rgb, float), rotate), "sum(RGB)")
    if img.size == 0:
        return 0
    score = np.percentile(img, 99, axis=1)
    return int(np.argmax(smooth_profile(score, smooth)))


def smooth_profile(profile, window: int = 1) -> np.ndarray:
    """Moving-average smoothing with edge padding; ``window <= 1`` is a no-op."""
    y = np.asarray(profile, float)
    w = int(window)
    if w <= 1 or y.size == 0:
        return y
    w = min(w, y.size)
    if w % 2 == 0:
        w += 1
    pad = w // 2
    padded = np.pad(y, pad, mode="edge")
    kernel = np.ones(w) / w
    return np.convolve(padded, kernel, mode="valid")


def _prominences(y: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """Topographic prominence of each peak index (same idea as scipy's)."""
    out = np.empty(idx.size)
    for k, i in enumerate(idx):
        h = y[i]
        j = i - 1
        left_min = h
        while j >= 0 and y[j] <= h:
            left_min = min(left_min, y[j])
            j -= 1
        j = i + 1
        right_min = h
        while j < y.size and y[j] <= h:
            right_min = min(right_min, y[j])
            j += 1
        out[k] = h - max(left_min, right_min)
    return out


def refine_peak(profile, index: int) -> float:
    """Sub-pixel peak position by fitting a parabola to the three top samples."""
    y = np.asarray(profile, float)
    i = int(index)
    if i <= 0 or i >= y.size - 1:
        return float(i)
    denom = y[i - 1] - 2.0 * y[i] + y[i + 1]
    if denom == 0:
        return float(i)
    shift = 0.5 * (y[i - 1] - y[i + 1]) / denom
    return float(i + np.clip(shift, -1.0, 1.0))


def find_profile_peaks(profile, min_prominence: float = 0.0,
                       min_distance: int = 1, smooth: int = 1,
                       max_peaks: int | None = None, subpixel: bool = True):
    """Find peaks in an intensity profile. Pure numpy, no scipy needed.

    ``min_prominence`` is in profile units; pass e.g. 0.05 * peak-to-peak range
    for a scale-free threshold. ``min_distance`` suppresses peaks closer than
    that many pixels, keeping the more prominent one. Returns a list of dicts
    with keys ``index`` (int), ``position`` (float, sub-pixel) and
    ``prominence``, sorted by position.
    """
    y_raw = np.asarray(profile, float)
    if y_raw.size < 3:
        return []
    y = smooth_profile(y_raw, smooth)

    # Plateau-aware local maxima: rising strictly, then non-rising.
    idx = []
    n = y.size
    i = 1
    while i < n - 1:
        if y[i] > y[i - 1]:
            j = i
            while j < n - 1 and y[j + 1] == y[i]:
                j += 1
            if j < n - 1 and y[j + 1] < y[j]:
                idx.append((i + j) // 2)
            i = j + 1
        else:
            i += 1
    idx = np.asarray(idx, int)
    if idx.size == 0:
        return []

    prom = _prominences(y, idx)
    keep = prom >= float(min_prominence)
    idx, prom = idx[keep], prom[keep]
    if idx.size == 0:
        return []

    order = np.argsort(prom)[::-1]  # strongest first
    chosen = []
    for k in order:
        if all(abs(int(idx[k]) - int(idx[c])) >= int(min_distance)
               for c in chosen):
            chosen.append(k)
        if max_peaks is not None and len(chosen) >= int(max_peaks):
            break
    chosen = sorted(chosen, key=lambda k: idx[k])

    peaks = []
    for k in chosen:
        i0 = int(idx[k])
        pos = refine_peak(y, i0) if subpixel else float(i0)
        peaks.append({"index": i0, "position": pos,
                      "prominence": float(prom[k]),
                      "intensity": float(y_raw[i0])})
    return peaks


def fit_calibration(pixels, wavelengths, model: str = "linear"):
    """Fit pixel -> wavelength as a polynomial. Returns (coeffs, r2, degree).

    ``coeffs`` are numpy.polyfit order (highest power first). Needs degree+1
    points; falls back to a lower degree when too few are supplied.
    """
    px = np.asarray(pixels, float)
    wl = np.asarray(wavelengths, float)
    if px.size != wl.size or px.size < 2:
        raise ValueError("Calibration needs at least 2 (pixel, wavelength) points.")
    if np.unique(px).size < 2:
        raise ValueError("Calibration needs at least 2 different pixel positions.")
    if np.unique(wl).size < 2:
        # Typically every peak was written with the same lambda still selected
        # in the combo box; the fit would be flat and the axis collapses.
        raise ValueError(
            f"every calibration line has the same wavelength ({wl[0]:g} nm); "
            "give each line its own wavelength.")
    degree = dict(CAL_MODELS).get(model, 1)
    degree = min(degree, px.size - 1)
    coeffs = np.polyfit(px, wl, degree)
    pred = np.polyval(coeffs, px)
    ss_res = float(np.sum((wl - pred) ** 2))
    ss_tot = float(np.sum((wl - wl.mean()) ** 2)) or 1.0
    return coeffs, 1.0 - ss_res / ss_tot, degree


def pixel_to_wavelength(pixel, calibration=None, model: str = "linear",
                        coeffs=None):
    """Convert a (possibly fractional) pixel position to nm, or None.

    Pass either ready-made ``coeffs`` (numpy.polyfit order) or a list of
    ``(pixel, wavelength)`` calibration pairs. Returns None when uncalibrated,
    so the widget can fall back to showing pixels.
    """
    if coeffs is None:
        if not calibration:
            return None
        coeffs, _, _ = fit_calibration([p for p, _ in calibration],
                                       [w for _, w in calibration], model)
    return float(np.polyval(coeffs, float(pixel)))


def monotonic_segment(x, keep=None):
    """Index range ``(lo, hi)`` of the monotonic branch of ``x`` to keep.

    A quadratic calibration fitted from a handful of lines often has its
    turning point inside the sensor, so pixels past it map back onto
    wavelengths already covered and the spectrum folds onto itself. This picks
    the monotonic run that carries the calibration pixels (``keep``), or the
    longest one when they are ambiguous. Returns ``(0, len(x))`` when ``x`` is
    already monotonic.
    """
    x = np.asarray(x, float)
    if x.size < 3:
        return 0, x.size
    d = np.diff(x)
    if np.all(d > 0) or np.all(d < 0):
        return 0, x.size
    sign = np.sign(d)
    breaks = [0] + [int(i) + 1 for i in np.nonzero(np.diff(sign))[0]] + [x.size - 1]
    segments = [(breaks[i], breaks[i + 1] + 1) for i in range(len(breaks) - 1)]
    keep = [] if keep is None else [float(k) for k in keep]

    def score(seg):
        lo, hi = seg
        inside = sum(1 for k in keep if lo <= k <= hi - 1)
        return (inside, hi - lo)

    return max(segments, key=score)


def image_to_spectrum(source, channel: str = "luminance",
                      row_center: float | None = None, row_frac: float = 0.15,
                      calibration=None, model: str = "linear",
                      name: str = "", flip: bool = False,
                      rotate: int = 0) -> dict:
    """Convert a spectrum photo into a spectrum dict.

    ``calibration`` is an optional list of (pixel, wavelength) pairs; with it,
    the x-axis is wavelength (nm), otherwise it is the pixel column index.
    ``rotate`` (0/90/180/270 clockwise) is applied first, for photos whose
    dispersion axis is not already horizontal. ``flip`` reverses the profile
    (when the grating spreads blue-to-red the other way). Returns
    make_spectrum(...) with an extra "calibration" entry.
    """
    rgb = load_rgb(source)
    profile = extract_profile(rgb, channel=channel, row_center=row_center,
                              row_frac=row_frac, rotate=rotate)
    if flip:
        profile = profile[::-1]
    px = np.arange(profile.size, dtype=float)

    x, x_label, cal, clipped = px, "pixel", None, None
    if calibration:
        pts = [(float(p), float(w)) for p, w in calibration]
        coeffs, r2, degree = fit_calibration([p for p, _ in pts],
                                             [w for _, w in pts], model)
        x = np.polyval(coeffs, px)
        x_label = "wavelength (nm)"
        cal = {"coeffs": coeffs, "r2": r2, "degree": degree, "points": pts}
        # A non-monotonic fit folds the spectrum onto itself; keep one branch.
        lo, hi = monotonic_segment(x, keep=[p for p, _ in pts])
        if (lo, hi) != (0, x.size):
            clipped = (int(lo), int(hi))
            x, profile, px = x[lo:hi], profile[lo:hi], px[lo:hi]
        if x.size and np.all(np.diff(x) < 0):  # keep ascending wavelength
            x, profile = x[::-1], profile[::-1]

    if not name:
        name = (os.path.splitext(os.path.basename(source))[0]
                if isinstance(source, str) else "spectrum")
    spec = make_spectrum(x, profile, name=name, x_label=x_label,
                         source=source if isinstance(source, str) else "")
    spec["calibration"] = cal
    spec["rotate"] = int(rotate) % 360
    spec["clipped"] = clipped
    return spec
