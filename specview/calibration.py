"""Calibration curve (檢量線): fit a standards series and predict concentrations.

The everyday quantitative workflow in food / instrumental analysis: measure a
set of standards of known concentration, read one number ("signal") from each
spectrum — absorbance at a wavelength, peak height, or integrated band area —
fit signal vs. concentration, then invert the model to estimate the
concentration of unknown samples (Beer–Lambert, A = ε·b·c).

This module is pure NumPy/SciPy (no Qt) so the maths can be unit-tested. It
reports the statistics taught in analytical-chemistry courses: slope and
intercept with their standard errors, R², the residual standard error s(y/x),
the limit of detection / quantification, and a confidence interval on each
predicted concentration. Reference: Miller & Miller, *Statistics and
Chemometrics for Analytical Chemistry*.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from .spectrum import Spectrum

_EPS = 1e-12
_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz  # trapz gone in numpy 2.3

# (label, value) for the signal-reading modes offered in the UI.
SIGNAL_MODES = [
    ("Value at x₀ (e.g. absorbance @ λ)", "value"),
    ("Peak height in range", "height"),
    ("Integrated area in range", "area"),
]

_NUM_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


# ============================================================ signal readers
def signal_at(spec: Spectrum, x0: float) -> float:
    """Signal = linearly-interpolated y at a single x position (absorbance @ λ)."""
    return float(spec.value_at(x0))


def _range(spec: Spectrum, lo: float, hi: float):
    lo, hi = sorted((float(lo), float(hi)))
    m = (spec.x >= lo) & (spec.x <= hi)
    return spec.x[m].astype(float), spec.y[m].astype(float)


def signal_height(spec: Spectrum, lo: float, hi: float,
                  subtract_baseline: bool = False) -> float:
    """Signal = peak height: the maximum y inside [lo, hi].

    With ``subtract_baseline`` the straight line joining the two range endpoints
    is removed first (a local baseline), so the height is measured from it.
    """
    xx, yy = _range(spec, lo, hi)
    if xx.size == 0:
        return float("nan")
    if subtract_baseline and xx.size >= 2:
        yy = yy - np.interp(xx, [xx[0], xx[-1]], [yy[0], yy[-1]])
    return float(np.max(yy))


def signal_area(spec: Spectrum, lo: float, hi: float,
                subtract_baseline: bool = False) -> float:
    """Signal = trapezoidal band area over [lo, hi] (optional endpoint baseline)."""
    xx, yy = _range(spec, lo, hi)
    if xx.size < 2:
        return float("nan")
    if subtract_baseline:
        yy = yy - np.interp(xx, [xx[0], xx[-1]], [yy[0], yy[-1]])
    return float(_trapz(yy, xx))


def read_signal(spec: Spectrum, mode: str, x0: float | None = None,
                lo: float | None = None, hi: float | None = None,
                subtract_baseline: bool = False) -> float:
    """Dispatch to the signal reader named by ``mode`` ('value'|'height'|'area')."""
    if mode == "value":
        return signal_at(spec, x0)
    if mode == "height":
        return signal_height(spec, lo, hi, subtract_baseline)
    if mode == "area":
        return signal_area(spec, lo, hi, subtract_baseline)
    raise ValueError(f"Unknown signal mode: {mode!r}")


def parse_concentration(name: str | None) -> float | None:
    """Best-effort concentration from a spectrum name.

    Returns the first number found ('std_5ppm' → 5.0, '0.25 mg/L' → 0.25) or
    None when the name has no number. The user can always correct it afterwards.
    """
    if not name:
        return None
    m = _NUM_RE.search(str(name))
    return float(m.group()) if m else None


# ============================================================ the model
@dataclass
class CalibrationModel:
    """A fitted calibration: y(signal) = c0 + c1·x + c2·x²  (x = concentration)."""

    degree: int
    coeffs: np.ndarray          # ascending powers [c0, c1, (c2)]; c0=0 if origin-forced
    through_origin: bool
    conc: np.ndarray            # standard concentrations actually used
    signal: np.ndarray          # standard signals actually used
    fitted: np.ndarray
    residuals: np.ndarray
    r_squared: float
    s_yx: float                 # residual standard error  √(SSE/dof)
    dof: int
    se: np.ndarray              # standard errors aligned with ``coeffs``
    conc_unit: str = ""
    signal_label: str = "signal"

    # ---- convenience accessors (the linear terms always exist) ----------
    @property
    def intercept(self) -> float:
        return float(self.coeffs[0])

    @property
    def slope(self) -> float:
        return float(self.coeffs[1])

    @property
    def se_intercept(self) -> float:
        return float(self.se[0])

    @property
    def se_slope(self) -> float:
        return float(self.se[1])

    @property
    def r(self) -> float:
        """Signed correlation coefficient (linear only)."""
        if self.degree != 1 or not np.isfinite(self.r_squared):
            return float("nan")
        return float(np.sign(self.slope) * np.sqrt(max(self.r_squared, 0.0)))

    @property
    def lod(self) -> float | None:
        """Limit of detection = 3.3·s(y/x)/slope (linear models)."""
        if self.degree == 1 and abs(self.slope) > _EPS and np.isfinite(self.s_yx):
            return 3.3 * self.s_yx / abs(self.slope)
        return None

    @property
    def loq(self) -> float | None:
        """Limit of quantification = 10·s(y/x)/slope (linear models)."""
        if self.degree == 1 and abs(self.slope) > _EPS and np.isfinite(self.s_yx):
            return 10.0 * self.s_yx / abs(self.slope)
        return None

    @property
    def conc_range(self) -> tuple[float, float]:
        return (float(self.conc.min()), float(self.conc.max())) if self.conc.size \
            else (0.0, 1.0)

    # ---- evaluation -----------------------------------------------------
    def predict_signal(self, conc):
        """Forward model: expected signal at the given concentration(s)."""
        return np.polyval(self.coeffs[::-1], np.asarray(conc, float))

    def _in_range(self, c: float) -> bool:
        lo, hi = self.conc_range
        return bool(np.isfinite(c) and lo - _EPS <= c <= hi + _EPS)

    def predict_concentration(self, signal_value: float,
                              replicates: int = 1) -> dict:
        """Inverse prediction: estimate concentration from a measured signal.

        Returns ``{conc, ci, in_range}`` where ``ci`` is the 95 % confidence
        half-width (Miller & Miller) for an ordinary linear fit, or None when a
        textbook interval does not apply (origin-forced or quadratic models).
        """
        s = float(signal_value)
        c0, c1 = float(self.coeffs[0]), float(self.coeffs[1])
        if self.degree == 1:
            if abs(c1) < _EPS:
                return {"conc": float("nan"), "ci": None, "in_range": False}
            x0 = (s - c0) / c1
            ci = self._linear_ci(s, x0, replicates)
            return {"conc": x0, "ci": ci, "in_range": self._in_range(x0)}
        # quadratic: solve c2·x² + c1·x + (c0 − s) = 0, pick the sensible root
        c2 = float(self.coeffs[2])
        if abs(c2) < _EPS:
            x0 = (s - c0) / c1 if abs(c1) > _EPS else float("nan")
            return {"conc": x0, "ci": None, "in_range": self._in_range(x0)}
        disc = c1 * c1 - 4.0 * c2 * (c0 - s)
        if disc < 0:
            return {"conc": float("nan"), "ci": None, "in_range": False}
        sq = float(np.sqrt(disc))
        roots = [(-c1 + sq) / (2 * c2), (-c1 - sq) / (2 * c2)]
        lo, hi = self.conc_range
        mid = 0.5 * (lo + hi)
        span = (hi - lo) or 1.0
        inside = [r for r in roots if lo - 0.5 * span <= r <= hi + 0.5 * span]
        x0 = min(inside or roots, key=lambda r: abs(r - mid))
        return {"conc": x0, "ci": None, "in_range": self._in_range(x0)}

    def _linear_ci(self, s: float, x0: float, replicates: int) -> float | None:
        """95 % confidence half-width on x0 for an ordinary (with-intercept) line."""
        if self.through_origin or self.dof <= 0 or not np.isfinite(self.s_yx):
            return None
        from scipy.stats import t as student_t
        n = int(self.conc.size)
        xbar = float(self.conc.mean())
        ybar = float(self.signal.mean())
        sxx = float(((self.conc - xbar) ** 2).sum())
        if sxx <= _EPS:
            return None
        m = max(1, int(replicates))
        sx0 = (self.s_yx / abs(self.slope)) * np.sqrt(
            1.0 / m + 1.0 / n + (s - ybar) ** 2 / (self.slope ** 2 * sxx))
        return float(student_t.ppf(0.975, self.dof) * sx0)

    def equation(self) -> str:
        """Human-readable model equation, e.g. ``y = 0.123·x + 0.004``."""
        c = self.coeffs
        if self.degree == 1:
            if self.through_origin:
                return f"y = {c[1]:.6g}·x"
            sgn = "+" if c[0] >= 0 else "−"
            return f"y = {c[1]:.6g}·x {sgn} {abs(c[0]):.6g}"
        s1 = "+" if c[1] >= 0 else "−"
        s0 = "+" if c[0] >= 0 else "−"
        return f"y = {c[2]:.6g}·x² {s1} {abs(c[1]):.6g}·x {s0} {abs(c[0]):.6g}"


# ============================================================ fitting
def fit_calibration(concentrations, signals, degree: int = 1,
                    through_origin: bool = False, conc_unit: str = "",
                    signal_label: str = "signal") -> CalibrationModel:
    """Least-squares fit of signal vs. concentration.

    degree: 1 (linear) or 2 (quadratic). ``through_origin`` forces the line
    through (0, 0) for the linear case (a no-intercept model). Non-finite pairs
    are dropped. Raises ValueError if too few standards remain.
    """
    if degree not in (1, 2):
        raise ValueError("degree must be 1 (linear) or 2 (quadratic).")
    c = np.asarray(concentrations, float).ravel()
    s = np.asarray(signals, float).ravel()
    n = min(c.size, s.size)
    c, s = c[:n], s[:n]
    good = np.isfinite(c) & np.isfinite(s)
    c, s = c[good], s[good]
    n = c.size

    origin = bool(through_origin and degree == 1)
    p = 1 if origin else degree + 1            # number of free parameters
    if n < p:
        kind = "quadratic" if degree == 2 else "linear"
        raise ValueError(f"Need at least {p} standards to fit a {kind} model; got {n}.")

    if origin:
        A = c[:, None]                          # slope only
    else:
        A = np.vander(c, degree + 1, increasing=True)   # [1, x, x²]
    beta, *_ = np.linalg.lstsq(A, s, rcond=None)
    coeffs = np.array([0.0, float(beta[0])]) if origin else beta.astype(float)

    fitted = A @ beta
    resid = s - fitted
    dof = n - p
    sse = float(resid @ resid)
    s_yx = float(np.sqrt(sse / dof)) if dof > 0 else float("nan")

    se_beta = np.full(p, np.nan)
    if dof > 0:
        try:
            cov = (s_yx ** 2) * np.linalg.inv(A.T @ A)
            se_beta = np.sqrt(np.clip(np.diag(cov), 0.0, None))
        except np.linalg.LinAlgError:
            pass
    se = np.array([0.0, float(se_beta[0])]) if origin else se_beta.astype(float)

    sst = float(((s - s.mean()) ** 2).sum())
    r_squared = (1.0 - sse / sst) if sst > _EPS else float("nan")

    return CalibrationModel(
        degree=degree, coeffs=coeffs, through_origin=origin, conc=c, signal=s,
        fitted=fitted, residuals=resid, r_squared=r_squared, s_yx=s_yx, dof=dof,
        se=se, conc_unit=conc_unit, signal_label=signal_label)
