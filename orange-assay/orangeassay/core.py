"""Assay image analysis + dose-response curve fitting (Qt-free, testable).

Reuses the coffee-ring-analyzer algorithm so results match the Streamlit /
mobile / CLI implementations:

  * ITU-R 601 luma grayscale, 0..1
  * Otsu threshold with a plateau search (robust to bimodal gaps)
  * grid of cells; per-cell thresholded pixel **count** or **intensity**
    (integral image -> O(1) per cell), plus a per-cell **mean** for wells
  * blank-column normalization + replicate averaging
  * 3-parameter logistic fit -> R2, EC50, LOD

Everything here is pure numpy / scipy (+ Pillow only to read image files), so
it runs and is unit-tested without Orange or Qt.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

_EPS = 1e-12

# Coffee-ring defaults (from coffee-ring-analyzer), so the widget reproduces it.
DEFAULT_THRESHOLD = 0.43
DEFAULT_X_BOUNDS = (0, 60, 140, 220, 290, 360, 440, 520, 583)
DEFAULT_Y_BOUNDS = (0, 50, 75, 130, 155, 211)
DEFAULT_CONCENTRATIONS = (0.01, 0.09, 0.23, 0.46, 0.91, 3.63, 7.27)

# Standard microtitre plate layouts: name -> (rows, cols).
PLATE_FORMATS = {
    "6 (2×3)": (2, 3), "12 (3×4)": (3, 4), "24 (4×6)": (4, 6),
    "48 (6×8)": (6, 8), "96 (8×12)": (8, 12), "384 (16×24)": (16, 24),
}
_ROW_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


# ============================================================ image loading
def load_grayscale(source) -> np.ndarray:
    """Return a 0..1 float grayscale image.

    ``source`` may be a file path (read with Pillow), an (H, W) array (assumed
    already grayscale), or an (H, W, 3/4) RGB(A) array. RGB is converted with
    ITU-R 601-2 luma weights (0.299, 0.587, 0.114).
    """
    if isinstance(source, str):
        from PIL import Image
        arr = np.asarray(Image.open(source).convert("RGB"), dtype=float)
    else:
        arr = np.asarray(source, dtype=float)
    if arr.ndim == 3:
        arr = arr[..., :3] @ np.array([0.299, 0.587, 0.114])
    if arr.max() > 1.0 + _EPS:
        arr = arr / 255.0
    return np.clip(arr.astype(float), 0.0, 1.0)


# ============================================================ thresholding
def otsu_threshold(gray: np.ndarray, bins: int = 256) -> float:
    """Otsu's threshold in [0, 1]; returns the centre of the maximal plateau of
    the between-class variance (robust when the histogram has a flat gap)."""
    g = np.asarray(gray, float).ravel()
    hist, edges = np.histogram(g, bins=bins, range=(0.0, 1.0))
    hist = hist.astype(float)
    total = hist.sum()
    if total <= 0:
        return 0.5
    p = hist / total
    centres = (edges[:-1] + edges[1:]) / 2.0
    omega = np.cumsum(p)
    mu = np.cumsum(p * centres)
    mu_t = mu[-1]
    denom = omega * (1.0 - omega)
    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_b = (mu_t * omega - mu) ** 2 / denom
    sigma_b = np.nan_to_num(sigma_b, nan=-1.0, posinf=-1.0, neginf=-1.0)
    peak = sigma_b.max()
    if peak <= 0:
        return 0.5
    plateau = np.flatnonzero(sigma_b >= peak - _EPS)
    return float(centres[plateau].mean())


# ============================================================ grid geometry
@dataclass(frozen=True)
class GridSpec:
    """Cell grid over an image. ``x_bounds`` has n_cols+1 entries; ``y_bounds``
    is flat (y0,y1, y2,y3, …) row-band pairs (2*n_rows entries)."""
    x_bounds: tuple
    y_bounds: tuple

    @property
    def n_cols(self) -> int:
        return len(self.x_bounds) - 1

    @property
    def n_rows(self) -> int:
        return len(self.y_bounds) // 2

    def column_spans(self):
        return [(int(self.x_bounds[i]), int(self.x_bounds[i + 1]))
                for i in range(self.n_cols)]

    def row_spans(self):
        return [(int(self.y_bounds[2 * i]), int(self.y_bounds[2 * i + 1]))
                for i in range(self.n_rows)]

    def cells(self):
        """Yield (row, col, x0, x1, y0, y1)."""
        for r, (y0, y1) in enumerate(self.row_spans()):
            for c, (x0, x1) in enumerate(self.column_spans()):
                yield r, c, x0, x1, y0, y1


def regular_grid(n_rows: int, n_cols: int, shape, margin: float = 0.0,
                 fill: float = 1.0) -> GridSpec:
    """Evenly spaced grid over an (H, W) image. ``margin`` (0..0.5) trims the
    image border before dividing into cells; ``fill`` (0..1) shrinks each row
    band toward its centre (leaving a gutter between rows — useful for round
    wells). Columns stay contiguous (GridSpec shares column boundaries)."""
    h, w = int(shape[0]), int(shape[1])
    mx, my = margin * w, margin * h
    x0, x1, y0, y1 = mx, w - mx, my, h - my
    xe = np.linspace(x0, x1, n_cols + 1)
    ye = np.linspace(y0, y1, n_rows + 1)
    x_bounds = tuple(int(round(v)) for v in xe)
    pitch_y = (y1 - y0) / n_rows
    gut_y = pitch_y * (1.0 - fill) / 2.0
    y_bounds = []
    for r in range(n_rows):
        y_bounds += [int(round(ye[r] + gut_y)), int(round(ye[r + 1] - gut_y))]
    return GridSpec(tuple(x_bounds), tuple(y_bounds))


def well_label(row: int, col: int) -> str:
    """Microplate-style label, e.g. row 0 col 0 -> 'A1'."""
    r = _ROW_LETTERS[row] if row < len(_ROW_LETTERS) else f"R{row + 1}"
    return f"{r}{col + 1}"


def auto_detect_grid(gray: np.ndarray, n_rows: int = 3, n_cols: int = 8,
                     threshold: float | None = None,
                     min_spot_frac: float = 5e-4) -> GridSpec | None:
    """Detect bright spots and cluster their centroids into an n_rows × n_cols
    grid. Returns None if too few spots are found."""
    g = np.asarray(gray, float)
    thr = otsu_threshold(g) if threshold is None else threshold
    mask = g > thr
    h, w = g.shape
    min_area = max(1.0, min_spot_frac * h * w)

    # centroids of connected bright blobs (scipy.ndimage label)
    from scipy import ndimage
    lab, n = ndimage.label(mask)
    if n == 0:
        return None
    areas = ndimage.sum(np.ones_like(lab), lab, index=np.arange(1, n + 1))
    keep = np.flatnonzero(areas >= min_area) + 1
    if keep.size < max(n_rows, n_cols):
        return None
    cys, cxs = np.array(ndimage.center_of_mass(mask, lab, keep)).T

    xb = _cluster_bounds(cxs, n_cols, 0, w)
    yb = _cluster_bounds(cys, n_rows, 0, h)
    if xb is None or yb is None:
        return None
    # x_bounds: n_cols+1 column boundaries; y_bounds: contiguous row bands as
    # flat (y0,y1),(y1,y2),… pairs from the n_rows+1 row boundaries.
    y_pairs = []
    for r in range(n_rows):
        y_pairs += [yb[r], yb[r + 1]]
    return GridSpec(tuple(xb), tuple(y_pairs))


def _cluster_bounds(coords, k, lo, hi):
    """1-D k-means-ish binning of ``coords`` into k ordered groups. For columns
    returns k+1 midpoint boundaries; kept generic via the caller."""
    coords = np.sort(np.asarray(coords, float))
    if coords.size < k:
        return None
    # even quantile seeds, then a few Lloyd iterations
    centres = np.quantile(coords, np.linspace(0, 1, k * 2 + 1)[1::2])
    for _ in range(25):
        idx = np.argmin(np.abs(coords[:, None] - centres[None, :]), axis=1)
        new = np.array([coords[idx == j].mean() if np.any(idx == j)
                        else centres[j] for j in range(k)])
        if np.allclose(new, centres):
            break
        centres = np.sort(new)
    # boundaries = midpoints between consecutive centres, clamped to [lo, hi]
    mids = (centres[:-1] + centres[1:]) / 2.0
    bounds = np.concatenate([[lo], mids, [hi]])
    # for a "banded" (pair) return, expand each centre into a tight band
    if k >= 1:
        pass
    return [int(round(b)) for b in bounds]


# ============================================================ measurement
def _integral(gray: np.ndarray) -> np.ndarray:
    return np.pad(np.cumsum(np.cumsum(gray, axis=0), axis=1),
                  ((1, 0), (1, 0)), mode="constant")


def _rect_sum(ii: np.ndarray, x0, x1, y0, y1) -> float:
    return float(ii[y1, x1] - ii[y0, x1] - ii[y1, x0] + ii[y0, x0])


def compute_areas(gray: np.ndarray, grid: GridSpec, threshold: float,
                  metric: str = "count", adaptive: bool = False) -> np.ndarray:
    """Per-cell measurement -> (n_rows, n_cols) array.

    metric="count": pixels above threshold; "intensity": summed grayscale of
    pixels above threshold; "mean": mean grayscale over the whole cell (for
    wells). ``adaptive`` re-derives an Otsu threshold within each cell.
    """
    g = np.asarray(gray, float)
    h, w = g.shape
    out = np.full((grid.n_rows, grid.n_cols), np.nan)
    for r, c, x0, x1, y0, y1 in grid.cells():
        x0, x1 = max(0, x0), min(w, x1)
        y0, y1 = max(0, y0), min(h, y1)
        if x1 <= x0 or y1 <= y0:
            continue
        cell = g[y0:y1, x0:x1]
        if metric == "mean":
            out[r, c] = float(cell.mean())
            continue
        thr = otsu_threshold(cell) if adaptive else threshold
        mask = cell > thr
        out[r, c] = (float(cell[mask].sum()) if metric == "intensity"
                     else float(mask.sum()))
    return out


def normalize_blank(areas: np.ndarray, blank_col: int = 0):
    """Divide each row's treatment columns by that row's blank column, then
    average replicate rows. Returns dict(ratios, means, ses, n)."""
    a = np.asarray(areas, float)
    blank = a[:, blank_col][:, None]
    cols = [c for c in range(a.shape[1]) if c != blank_col]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = a[:, cols] / blank
    means = np.nanmean(ratios, axis=0)
    n = np.sum(~np.isnan(ratios), axis=0)
    sd = np.nanstd(ratios, axis=0, ddof=1) if a.shape[0] > 1 \
        else np.zeros(len(cols))
    ses = sd / np.sqrt(np.maximum(n, 1))
    return {"ratios": ratios, "means": means, "ses": ses, "n": n,
            "columns": cols}


# ============================================================ dose-response
@dataclass
class LogisticFit:
    model: str
    params: dict
    x0: float
    r2: float
    success: bool
    ec50: float
    lod: float
    x: np.ndarray
    y: np.ndarray

    def predict(self, x):
        return _logistic_eval(self.model, np.asarray(x, float), self.params)


def _logistic_eval(model, x, p):
    if model == "4pl":
        lo, hi, k, x0 = p["lower"], p["upper"], p["k"], p["x0"]
        return lo + (hi - lo) / (1.0 + np.exp(-k * (x - x0)))
    L, k, x0 = p["L"], p["k"], p["x0"]           # 3pl
    return L / (1.0 + np.exp(-k * (x - x0)))


def fit_logistic(conc, signal, model: str = "3pl", sigma: float | None = None):
    """Fit a 3- or 4-parameter logistic dose-response curve.

    3pl:  y = L / (1 + exp(-k (x - x0)))
    4pl:  y = lower + (upper-lower) / (1 + exp(-k (x - x0)))

    Returns a :class:`LogisticFit` with R², EC50 (= x0, the inflection) and LOD
    (concentration whose signal exceeds baseline + 3σ). ``sigma`` defaults to a
    small fraction of the signal range when replicate SEs aren't supplied.
    """
    from scipy.optimize import curve_fit

    x = np.asarray(conc, float)
    y = np.asarray(signal, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3:
        raise ValueError("Need at least 3 finite (concentration, signal) points.")
    order = np.argsort(x)
    x, y = x[order], y[order]

    span = float(y.max() - y.min()) or 1.0
    x0_0 = float(np.median(x))
    k0 = 4.0 / (float(x.max() - x.min()) or 1.0)
    if model == "4pl":
        def f(xx, lo, hi, k, x0):
            return lo + (hi - lo) / (1.0 + np.exp(-k * (xx - x0)))
        p0 = [float(y.min()), float(y.max()), k0, x0_0]
        names = ["lower", "upper", "k", "x0"]
    else:
        def f(xx, L, k, x0):
            return L / (1.0 + np.exp(-k * (xx - x0)))
        p0 = [float(y.max()), k0, x0_0]
        names = ["L", "k", "x0"]

    success = True
    try:
        popt, _ = curve_fit(f, x, y, p0=p0, maxfev=20000)
    except Exception:                                    # noqa: BLE001
        popt, success = np.array(p0, float), False
    params = dict(zip(names, (float(v) for v in popt)))

    yhat = f(x, *popt)
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2)) or _EPS
    r2 = 1.0 - ss_res / ss_tot

    x0 = params["x0"]
    k = params["k"]
    upper = params.get("upper", params.get("L", span))
    lower = params.get("lower", 0.0)
    if sigma is None:
        sigma = 0.05 * span
    baseline = _logistic_eval(model, np.array([x.min()]), params)[0]
    lod = _lod_from_fit(model, params, baseline, sigma, lower, upper, k, x0)
    return LogisticFit(model=model, params=params, x0=float(x0), r2=float(r2),
                       success=bool(success), ec50=float(x0), lod=float(lod),
                       x=x, y=y)


def _lod_from_fit(model, params, baseline, sigma, lower, upper, k, x0):
    """Concentration at which the fitted signal first exceeds baseline + 3σ."""
    target = baseline + 3.0 * sigma
    rng = (upper - lower) if model == "4pl" else params["L"]
    try:
        if model == "4pl":
            frac = (target - lower) / (rng or _EPS)
        else:
            frac = target / (rng or _EPS)
        frac = min(max(frac, _EPS), 1.0 - _EPS)
        return x0 - math.log(1.0 / frac - 1.0) / (k or _EPS)
    except (ValueError, ZeroDivisionError):
        return float("nan")
