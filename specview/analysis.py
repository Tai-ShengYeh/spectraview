"""Peak analysis: detection (with FWHM), curve fitting / deconvolution, integration.

All functions take a :class:`Spectrum` and return plain result objects; they do
not mutate the input.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks as _sp_find_peaks
from scipy.signal import peak_widths, savgol_filter

from .spectrum import Spectrum

_EPS = 1e-12
_trapz = getattr(np, "trapezoid", getattr(np, "trapz"))
_GAUSS_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))   # FWHM = 2.3548 * sigma
_GAUSS_AREA = np.sqrt(np.pi / (4.0 * np.log(2.0)))  # area ≈ 1.0645 * height * FWHM


# ============================================================ data objects
@dataclass
class Peak:
    center: float
    height: float
    fwhm: float
    prominence: float
    area: float          # Gaussian-equivalent estimate: height * FWHM * 1.0645
    index: int


@dataclass
class FitComponent:
    center: float
    amplitude: float     # peak height
    fwhm: float
    area: float          # numerically integrated component area
    eta: float | None = None   # pseudo-Voigt mixing (None for pure G/L)


@dataclass
class FitResult:
    model: str
    components: list          # list[FitComponent]
    x: np.ndarray
    total: np.ndarray        # fitted sum (peaks + baseline)
    comp_curves: list        # list[np.ndarray], one per component (no baseline)
    baseline: np.ndarray | None
    r_squared: float


# ============================================================ peak shapes
def gaussian(x, amp, center, sigma):
    return amp * np.exp(-((x - center) ** 2) / (2.0 * sigma ** 2))


def lorentzian(x, amp, center, gamma):
    return amp / (1.0 + ((x - center) / gamma) ** 2)


def pseudo_voigt(x, amp, center, fwhm, eta):
    g = np.exp(-4.0 * np.log(2.0) * (x - center) ** 2 / fwhm ** 2)
    lo = 1.0 / (1.0 + 4.0 * (x - center) ** 2 / fwhm ** 2)
    return amp * (eta * lo + (1.0 - eta) * g)


_PER_PARAMS = {"gaussian": 3, "lorentzian": 3, "pseudovoigt": 4}


def _eval_peak(model: str, x, ps):
    if model == "gaussian":
        return gaussian(x, ps[0], ps[1], ps[2])
    if model == "lorentzian":
        return lorentzian(x, ps[0], ps[1], ps[2])
    return pseudo_voigt(x, ps[0], ps[1], ps[2], ps[3])


# ============================================================ peak finding
def find_peaks(spec: Spectrum, min_height_frac: float = 0.05,
               min_prominence_frac: float = 0.03, min_distance: float = 0.0,
               valleys: bool = False, smooth_window: int = 0) -> list[Peak]:
    """Detect peaks (or valleys) and measure their FWHM.

    Fractions are relative to the signal's full range. ``min_distance`` is in
    x-units. Returns peaks sorted by position.
    """
    x = spec.x.astype(float)
    y0 = spec.y.astype(float)
    y = y0
    if smooth_window and smooth_window >= 3 and smooth_window < y0.size:
        w = smooth_window + (smooth_window + 1) % 2  # force odd
        y = savgol_filter(y0, w, min(3, w - 1), mode="interp")
    sig = -y if valleys else y
    rng = float(sig.max() - sig.min()) or 1.0
    height = sig.min() + min_height_frac * rng
    prominence = max(min_prominence_frac * rng, _EPS)
    dx = abs(spec.uniform_step()) or 1.0
    distance = max(1, int(round(min_distance / dx))) if min_distance > 0 else 1

    idx, props = _sp_find_peaks(sig, height=height, prominence=prominence,
                                distance=distance)
    if idx.size == 0:
        return []
    widths, _, lips, rips = peak_widths(sig, idx, rel_height=0.5)
    axis = np.arange(x.size)
    x_left = np.interp(lips, axis, x)
    x_right = np.interp(rips, axis, x)
    fwhm = np.abs(x_right - x_left)

    peaks = []
    for k, i in enumerate(idx):
        h = float(y0[i])
        f = float(fwhm[k])
        peaks.append(Peak(center=float(x[i]), height=h, fwhm=f,
                          prominence=float(props["prominences"][k]),
                          area=abs(h) * f * float(_GAUSS_AREA), index=int(i)))
    peaks.sort(key=lambda p: p.center)
    return peaks


# ============================================================ curve fitting
def _build_model(model: str, n_peaks: int, fit_baseline: bool):
    per = _PER_PARAMS[model]

    def model_fn(x, *params):
        y = np.zeros_like(x, dtype=float)
        for j in range(n_peaks):
            y = y + _eval_peak(model, x, params[j * per:(j + 1) * per])
        if fit_baseline:
            slope, intercept = params[n_peaks * per], params[n_peaks * per + 1]
            y = y + slope * x + intercept
        return y

    return model_fn


def fit_peaks(spec: Spectrum, model: str = "gaussian", max_peaks: int = 6,
              fit_baseline: bool = False, init_peaks: list[Peak] | None = None,
              x_range: tuple[float, float] | None = None) -> FitResult:
    """Fit a sum of peak profiles to the spectrum (auto-initialised by find_peaks).

    model: 'gaussian' | 'lorentzian' | 'pseudovoigt'.
    """
    if model not in _PER_PARAMS:
        raise ValueError(f"Unknown peak model: {model}")
    x = spec.x.astype(float)
    y = spec.y.astype(float)
    if x_range is not None:
        lo, hi = sorted(x_range)
        m = (x >= lo) & (x <= hi)
        x, y = x[m], y[m]
    if x.size < 5:
        raise ValueError("Too few points in the fit range.")

    if init_peaks is None:
        sub = Spectrum(x, y, x_unit=spec.x_unit, y_unit=spec.y_unit)
        init_peaks = find_peaks(sub, min_height_frac=0.02, min_prominence_frac=0.02)
    if not init_peaks:
        raise ValueError("No peaks detected to initialise the fit.")
    init_peaks = sorted(init_peaks, key=lambda p: p.prominence, reverse=True)[:max_peaks]
    init_peaks = sorted(init_peaks, key=lambda p: p.center)

    span = float(x.max() - x.min()) or 1.0
    ymax = float(np.abs(y).max()) or 1.0
    per = _PER_PARAMS[model]
    p0, lo, hi = [], [], []
    for p in init_peaks:
        amp = p.height if abs(p.height) > _EPS else ymax * 0.1
        if model == "gaussian":
            w = max(p.fwhm / _GAUSS_FWHM, span * 1e-3)
        elif model == "lorentzian":
            w = max(p.fwhm / 2.0, span * 1e-3)
        else:
            w = max(p.fwhm, span * 1e-3)
        p0 += [amp, p.center, w]
        lo += [0.0, float(x.min()), span * 1e-4]
        hi += [5.0 * ymax + _EPS, float(x.max()), span]
        if model == "pseudovoigt":
            p0.append(0.5)
            lo.append(0.0)
            hi.append(1.0)
    if fit_baseline:
        p0 += [0.0, float(np.median(y))]
        lo += [-np.inf, -np.inf]
        hi += [np.inf, np.inf]

    model_fn = _build_model(model, len(init_peaks), fit_baseline)
    popt, _ = curve_fit(model_fn, x, y, p0=p0, bounds=(lo, hi), maxfev=20000)

    comps, comp_curves = [], []
    for j in range(len(init_peaks)):
        ps = popt[j * per:(j + 1) * per]
        curve = _eval_peak(model, x, ps)
        comp_curves.append(curve)
        area = float(_trapz(curve, x))
        if model == "gaussian":
            fwhm = float(_GAUSS_FWHM * abs(ps[2]))
            eta = None
        elif model == "lorentzian":
            fwhm = float(2.0 * abs(ps[2]))
            eta = None
        else:
            fwhm = float(abs(ps[2]))
            eta = float(ps[3])
        comps.append(FitComponent(center=float(ps[1]), amplitude=float(ps[0]),
                                  fwhm=fwhm, area=area, eta=eta))
    comps.sort(key=lambda c: c.center)

    baseline = None
    if fit_baseline:
        baseline = popt[-2] * x + popt[-1]
    total = model_fn(x, *popt)
    ss_res = float(np.sum((y - total) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2)) or 1.0
    return FitResult(model=model, components=comps, x=x, total=total,
                     comp_curves=comp_curves, baseline=baseline,
                     r_squared=1.0 - ss_res / ss_tot)


def component_to_spectrum(spec: Spectrum, x: np.ndarray, curve: np.ndarray,
                          label: str) -> Spectrum:
    """Wrap a fitted component curve as a Spectrum (for adding to the document)."""
    out = Spectrum(x=x, y=curve, name=label, x_unit=spec.x_unit, y_unit=spec.y_unit)
    out.meta["fit_component"] = True
    return out


# ============================================================ integration
def integrate(spec: Spectrum, x_lo: float, x_hi: float,
              baseline: str = "none") -> dict:
    """Integrate the spectrum over [x_lo, x_hi].

    baseline='linear' subtracts the straight line joining the two endpoints.
    Returns area, centroid and the endpoints actually used.
    """
    lo, hi = sorted((x_lo, x_hi))
    m = (spec.x >= lo) & (spec.x <= hi)
    xx = spec.x[m].astype(float)
    yy = spec.y[m].astype(float).copy()
    if xx.size < 2:
        raise ValueError("Integration range contains too few points.")
    base = None
    if baseline == "linear":
        base = np.interp(xx, [xx[0], xx[-1]], [yy[0], yy[-1]])
        yy = yy - base
    area = float(_trapz(yy, xx))
    centroid = float(_trapz(xx * yy, xx) / area) if abs(area) > _EPS else float("nan")
    return {"area": area, "centroid": centroid,
            "x_lo": float(xx[0]), "x_hi": float(xx[-1]),
            "baseline": baseline}
