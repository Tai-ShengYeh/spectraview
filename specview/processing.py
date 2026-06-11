"""Signal processing for spectra.

Every public function takes a :class:`Spectrum` and returns a NEW Spectrum,
leaving the input untouched. Array-level helpers are kept private (``_``).
All routines that assume a uniform x-grid resample internally when needed.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.signal import savgol_filter
from scipy.sparse.linalg import spsolve

from .spectrum import Spectrum

_EPS = 1e-12
# np.trapz was deprecated in NumPy 2.0 in favour of np.trapezoid.
_trapz = getattr(np, "trapezoid", getattr(np, "trapz"))


# ----------------------------------------------------------------- helpers
def _odd_window(n: int, npoints: int) -> int:
    """Clamp a Savitzky-Golay window to a valid odd value <= npoints."""
    n = int(n)
    if n < 3:
        n = 3
    if n % 2 == 0:
        n += 1
    if n > npoints:
        n = npoints if npoints % 2 == 1 else npoints - 1
    return max(n, 3)


def _to_uniform(spec: Spectrum):
    """Return (x_uniform, y_on_uniform, was_resampled)."""
    if spec.is_uniform():
        return spec.x, spec.y, False
    n = spec.npoints
    xu = np.linspace(spec.x[0], spec.x[-1], n)
    yu = np.interp(xu, spec.x, spec.y)
    return xu, yu, True


# --------------------------------------------------------------- smoothing
def savitzky_golay(spec: Spectrum, window: int = 11, polyorder: int = 3) -> Spectrum:
    """Savitzky-Golay smoothing (preserves peak shape better than a box filter)."""
    xu, yu, resampled = _to_uniform(spec)
    win = _odd_window(window, yu.size)
    po = min(int(polyorder), win - 1)
    ys = savgol_filter(yu, win, po, mode="interp")
    x_out = spec.x if not resampled else xu
    out = spec.replace_data(x_out, np.interp(spec.x, xu, ys) if resampled else ys)
    out.meta["last_op"] = f"Savitzky-Golay (win={win}, poly={po})"
    return out


def moving_average(spec: Spectrum, window: int = 5) -> Spectrum:
    """Simple boxcar moving average."""
    w = max(1, int(window))
    kernel = np.ones(w) / w
    ys = np.convolve(spec.y, kernel, mode="same")
    out = spec.replace_data(spec.x, ys)
    out.meta["last_op"] = f"Moving average (win={w})"
    return out


# ------------------------------------------------------------- derivatives
def derivative(spec: Spectrum, order: int = 1, window: int = 11,
               polyorder: int = 3) -> Spectrum:
    """n-th derivative via Savitzky-Golay (order 1..4)."""
    order = int(np.clip(order, 1, 4))
    xu, yu, resampled = _to_uniform(spec)
    win = _odd_window(window, yu.size)
    po = min(max(int(polyorder), order), win - 1)
    delta = float(np.median(np.diff(xu))) or 1.0
    yd = savgol_filter(yu, win, po, deriv=order, delta=delta, mode="interp")
    y_out = np.interp(spec.x, xu, yd) if resampled else yd
    out = spec.replace_data(spec.x, y_out)
    out.y_unit = "a.u."
    out.name = f"{spec.name} d{order}"
    out.meta["last_op"] = f"{order}-order derivative (SG win={win}, poly={po})"
    return out


# ---------------------------------------------------------------- baseline
def _rubberband(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Lower convex-hull ('rubber band') baseline."""
    n = x.size
    if n < 3:
        return np.full_like(y, y.min())
    # Monotone chain lower hull on (x, y).
    pts = list(range(n))
    hull = []
    for i in pts:
        while len(hull) >= 2:
            o, a = hull[-2], hull[-1]
            cross = (x[a] - x[o]) * (y[i] - y[o]) - (y[a] - y[o]) * (x[i] - x[o])
            if cross <= 0:
                hull.pop()
            else:
                break
        hull.append(i)
    hx, hy = x[hull], y[hull]
    return np.interp(x, hx, hy)


def baseline_rubberband(spec: Spectrum) -> Spectrum:
    base = _rubberband(spec.x, spec.y)
    out = spec.replace_data(spec.x, spec.y - base)
    out.meta["last_op"] = "Baseline: rubberband"
    return out


def baseline_polynomial(spec: Spectrum, order: int = 3, n_iter: int = 25) -> Spectrum:
    """Iterative polynomial baseline (ModPoly): fit, clip points above fit, repeat."""
    x, y = spec.x, spec.y.copy()
    # Normalise x to [-1, 1] for a well-conditioned Vandermonde fit.
    xs = np.linspace(-1.0, 1.0, x.size)
    work = y.copy()
    base = np.zeros_like(y)
    for _ in range(max(1, int(n_iter))):
        coeffs = np.polyfit(xs, work, int(order))
        base = np.polyval(coeffs, xs)
        new_work = np.minimum(work, base)
        if np.allclose(new_work, work, atol=_EPS):
            work = new_work
            break
        work = new_work
    out = spec.replace_data(x, y - base)
    out.meta["last_op"] = f"Baseline: polynomial (order={order})"
    return out


def baseline_als(spec: Spectrum, lam: float = 1e5, p: float = 0.01,
                 n_iter: int = 10) -> Spectrum:
    """Asymmetric Least Squares baseline (Eilers & Boelens, 2005).

    ``lam`` controls smoothness (1e2..1e9); ``p`` the asymmetry (0.001..0.1).
    """
    y = spec.y
    n = y.size
    if n < 3:
        return spec.copy()
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(n, n - 2))
    D = lam * (D @ D.transpose())
    w = np.ones(n)
    z = y.copy()
    for _ in range(max(1, int(n_iter))):
        W = sparse.spdiags(w, 0, n, n)
        z = spsolve((W + D).tocsc(), w * y)
        w = p * (y > z) + (1.0 - p) * (y < z)
    out = spec.replace_data(spec.x, y - z)
    out.meta["last_op"] = f"Baseline: ALS (lam={lam:g}, p={p:g})"
    return out


# ----------------------------------------------------------- normalisation
def normalize(spec: Spectrum, method: str = "max", x0: float | None = None) -> Spectrum:
    """Normalise a spectrum.

    method: 'max' (peak=1), 'minmax' (0..1), 'area' (∫=1), 'vector' (L2=1),
            'value' (y at x0 = 1).
    """
    y = spec.y.astype(float)
    if method == "max":
        denom = np.max(np.abs(y))
        ys = y / denom if denom > _EPS else y
    elif method == "minmax":
        lo, hi = y.min(), y.max()
        ys = (y - lo) / (hi - lo) if (hi - lo) > _EPS else y - lo
    elif method == "area":
        area = _trapz(np.abs(y), spec.x)
        ys = y / area if abs(area) > _EPS else y
    elif method == "vector":
        norm = np.sqrt(np.sum(y ** 2))
        ys = y / norm if norm > _EPS else y
    elif method == "value":
        if x0 is None:
            raise ValueError("normalize(method='value') needs x0.")
        v = float(np.interp(x0, spec.x, y))
        ys = y / v if abs(v) > _EPS else y
    else:
        raise ValueError(f"Unknown normalize method: {method}")
    out = spec.replace_data(spec.x, ys)
    out.meta["last_op"] = f"Normalize: {method}"
    return out


# --------------------------------------------------- scatter correction
def snv(spec: Spectrum) -> Spectrum:
    """Standard Normal Variate: (y - mean) / std."""
    y = spec.y
    mu, sd = y.mean(), y.std()
    ys = (y - mu) / sd if sd > _EPS else y - mu
    out = spec.replace_data(spec.x, ys)
    out.meta["last_op"] = "SNV"
    return out


def detrend(spec: Spectrum, order: int = 1) -> Spectrum:
    """Remove a polynomial trend of the given order."""
    xs = np.linspace(-1.0, 1.0, spec.x.size)
    coeffs = np.polyfit(xs, spec.y, int(order))
    trend = np.polyval(coeffs, xs)
    out = spec.replace_data(spec.x, spec.y - trend)
    out.meta["last_op"] = f"Detrend (order={order})"
    return out


def msc(spec: Spectrum, reference: Spectrum) -> Spectrum:
    """Multiplicative Scatter Correction against a reference spectrum.

    The reference is interpolated onto ``spec``'s x-grid, then a linear fit
    y ≈ a + b·ref is removed: y_corr = (y - a) / b.
    """
    ref = np.interp(spec.x, reference.x, reference.y)
    b, a = np.polyfit(ref, spec.y, 1)
    ys = (spec.y - a) / b if abs(b) > _EPS else spec.y - a
    out = spec.replace_data(spec.x, ys)
    out.meta["last_op"] = "MSC"
    return out


def msc_set(specs: list[Spectrum], reference: Spectrum | None = None) -> list[Spectrum]:
    """Apply MSC to a list of spectra, using the mean spectrum as reference if none given."""
    if not specs:
        return []
    if reference is None:
        reference = average(specs, method="mean")
    return [msc(s, reference) for s in specs]


# ------------------------------------------------------------ transforms
def cut(spec: Spectrum, x_lo: float, x_hi: float, keep: bool = True) -> Spectrum:
    """Keep (default) or remove the x-range [x_lo, x_hi]."""
    lo, hi = sorted((x_lo, x_hi))
    inside = (spec.x >= lo) & (spec.x <= hi)
    mask = inside if keep else ~inside
    out = spec.replace_data(spec.x[mask], spec.y[mask])
    out.meta["last_op"] = f"Cut {'keep' if keep else 'remove'} [{lo:g}, {hi:g}]"
    return out


def interpolate(spec: Spectrum, step: float | None = None,
                n: int | None = None) -> Spectrum:
    """Resample onto a uniform grid by step width or number of points."""
    x0, x1 = spec.x[0], spec.x[-1]
    if step:
        new_x = np.arange(x0, x1 + step * 0.5, step)
    elif n:
        new_x = np.linspace(x0, x1, int(n))
    else:
        new_x = np.linspace(x0, x1, spec.npoints)
    new_y = np.interp(new_x, spec.x, spec.y)
    out = spec.replace_data(new_x, new_y)
    out.meta["last_op"] = "Interpolate / resample"
    return out


# ------------------------------------------------------------ arithmetic
def _common_grid(a: Spectrum, b: Spectrum):
    """Overlapping x-grid of two spectra with b interpolated onto it."""
    lo = max(a.x[0], b.x[0])
    hi = min(a.x[-1], b.x[-1])
    if hi <= lo:
        raise ValueError("Spectra do not overlap in x.")
    mask = (a.x >= lo) & (a.x <= hi)
    gx = a.x[mask]
    ay = a.y[mask]
    by = np.interp(gx, b.x, b.y)
    return gx, ay, by


def combine(a: Spectrum, b: Spectrum, op: str) -> Spectrum:
    """Arithmetic between two spectra: op in {'add','sub','mul','div'}."""
    gx, ay, by = _common_grid(a, b)
    if op == "add":
        gy = ay + by
    elif op == "sub":
        gy = ay - by
    elif op == "mul":
        gy = ay * by
    elif op == "div":
        gy = ay / np.where(np.abs(by) < _EPS, np.nan, by)
    else:
        raise ValueError(f"Unknown op: {op}")
    out = a.replace_data(gx, gy)
    sym = {"add": "+", "sub": "−", "mul": "×", "div": "÷"}[op]
    out.name = f"({a.name} {sym} {b.name})"
    out.color = None
    out.meta["last_op"] = f"Arithmetic: {a.name} {sym} {b.name}"
    return out


def average(specs: list[Spectrum], method: str = "mean") -> Spectrum:
    """Average a list of spectra (interpolated onto the first's overlapping grid)."""
    if not specs:
        raise ValueError("average() needs at least one spectrum.")
    if len(specs) == 1:
        return specs[0].copy()
    lo = max(s.x[0] for s in specs)
    hi = min(s.x[-1] for s in specs)
    if hi <= lo:
        raise ValueError("Spectra do not share a common x-range.")
    base = specs[0]
    gx = base.x[(base.x >= lo) & (base.x <= hi)]
    stack = np.vstack([np.interp(gx, s.x, s.y) for s in specs])
    gy = np.median(stack, axis=0) if method == "median" else np.mean(stack, axis=0)
    out = base.replace_data(gx, gy)
    out.name = f"{method} of {len(specs)} spectra"
    out.color = None
    out.meta["last_op"] = f"Average ({method}, n={len(specs)})"
    return out
