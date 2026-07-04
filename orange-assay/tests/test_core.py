"""Headless tests for orangeassay.core (no Orange / Qt needed).

Run:  python tests/test_core.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from orangeassay import core  # noqa: E402

PASS = FAIL = 0


def check(name, cond):
    global PASS, FAIL
    PASS, FAIL = PASS + cond, FAIL + (not cond)
    print(f"  {'ok  ' if cond else 'FAIL'} {name}")


def _synthetic_plate(n_rows=3, n_cols=8, cell=36, gap=24):
    """Bright spots on a dark field; spot brightness rises with column."""
    pitch = cell + gap
    h, w = n_rows * pitch, n_cols * pitch
    plate = np.full((h, w), 0.08)
    for r in range(n_rows):
        for c in range(n_cols):
            b = 0.5 + 0.4 * (c / (n_cols - 1))
            y0, x0 = r * pitch + gap // 2, c * pitch + gap // 2
            plate[y0:y0 + cell, x0:x0 + cell] = b
    return plate


print("== grayscale (ITU-R 601) ==")
rgb = np.zeros((4, 4, 3)); rgb[..., 1] = 255
check("green -> 0.587 luma", abs(core.load_grayscale(rgb)[0, 0] - 0.587) < 1e-6)
check("0..1 range clamp", core.load_grayscale(np.full((3, 3), 500.0)).max() <= 1.0)
rgb2 = np.dstack([np.full((2, 2), 255.0)] * 3)
check("white -> ~1.0", abs(core.load_grayscale(rgb2)[0, 0] - 1.0) < 1e-6)

print("== otsu threshold ==")
bimodal = np.concatenate([np.full(500, 0.15), np.full(500, 0.85)])
check("bimodal split ~0.5", 0.3 < core.otsu_threshold(bimodal) < 0.7)
check("flat image safe", 0.0 <= core.otsu_threshold(np.full(100, 0.4)) <= 1.0)

print("== regular grid ==")
grid = core.regular_grid(3, 8, (180, 480))
check("n_rows/n_cols", grid.n_rows == 3 and grid.n_cols == 8)
check("column_spans count", len(grid.column_spans()) == 8)
check("row_spans count", len(grid.row_spans()) == 3)
check("cells count", len(list(grid.cells())) == 24)
gfill = core.regular_grid(2, 2, (100, 100), fill=0.5)
(y0, y1) = gfill.row_spans()[0]
check("fill<1 gutters the row band", (y1 - y0) < 50)
gmarg = core.regular_grid(2, 2, (100, 100), margin=0.1)
check("margin trims the border", gmarg.column_spans()[0][0] >= 9)

print("== measurement ==")
plate = _synthetic_plate()
grid = core.regular_grid(3, 8, plate.shape)
counts = core.compute_areas(plate, grid, threshold=0.6, metric="count")
check("count: brighter cols have >= area", counts[0, 7] >= counts[0, 1])
check("count: dim col below thr is 0", counts[0, 0] == 0)
inten = core.compute_areas(plate, grid, threshold=0.3, metric="intensity")
check("intensity rises with brightness", inten[0, 7] > inten[0, 1])
means = core.compute_areas(plate, grid, threshold=0.3, metric="mean")
check("mean rises with brightness", means[0, 7] > means[0, 0])
adap = core.compute_areas(plate, grid, threshold=0.5, metric="count", adaptive=True)
check("adaptive returns finite", np.all(np.isfinite(adap)))

print("== integral image correctness ==")
img = np.arange(1, 13, dtype=float).reshape(3, 4) / 12.0
g1 = core.GridSpec((0, 4), (0, 3))                 # whole image, one cell
whole = core.compute_areas(img, g1, threshold=0.0, metric="intensity")[0, 0]
check("intensity == direct sum", abs(whole - img.sum()) < 1e-9)

print("== blank normalization ==")
areas = np.array([[10.0, 20.0, 30.0], [10.0, 40.0, 60.0]])
norm = core.normalize_blank(areas, blank_col=0)
check("2 treatment cols", len(norm["means"]) == 2)
check("row-wise blank ratio then mean",
      abs(norm["means"][0] - np.mean([20 / 10, 40 / 10])) < 1e-9)
check("columns skip blank", norm["columns"] == [1, 2])

print("== logistic fit (3PL / 4PL) ==")
conc = np.array([0.01, 0.09, 0.23, 0.46, 0.91, 3.63, 7.27])
lx = np.log10(conc)
true = 2.0 / (1 + np.exp(-2.5 * (lx - np.log10(0.5))))
fit = core.fit_logistic(lx, true, model="3pl")
check("3PL R^2 ~ 1", fit.r2 > 0.99)
check("3PL converged", fit.success)
check("EC50 == x0 (inflection)", abs(fit.ec50 - fit.x0) < 1e-9)
check("EC50 near log10(0.5)", abs(fit.ec50 - np.log10(0.5)) < 0.15)
check("predict matches at points",
      np.allclose(fit.predict(lx), true, atol=0.05))
fit4 = core.fit_logistic(lx, true + 0.3, model="4pl")
check("4PL R^2 ~ 1", fit4.r2 > 0.99)
check("4PL has lower/upper", "lower" in fit4.params and "upper" in fit4.params)
try:
    core.fit_logistic([1, 2], [1, 2], model="3pl")
    check("too few points raises", False)
except ValueError:
    check("too few points raises", True)
noisy = core.fit_logistic(lx, np.zeros_like(lx), model="3pl")
check("degenerate signal doesn't crash", np.isfinite(noisy.r2) or noisy.r2 == 0)

print("== auto grid detection ==")
ag = core.auto_detect_grid(plate, n_rows=3, n_cols=8, threshold=0.3)
check("auto grid shape", ag is not None and ag.n_rows == 3 and ag.n_cols == 8)
blank_img = np.full((100, 100), 0.1)
check("no spots -> None", core.auto_detect_grid(blank_img, 3, 8) is None)

print("== well labels & plate formats ==")
check("A1", core.well_label(0, 0) == "A1")
check("H12", core.well_label(7, 11) == "H12")
check("96-well format", core.PLATE_FORMATS["96 (8×12)"] == (8, 12))
check("384-well format", core.PLATE_FORMATS["384 (16×24)"] == (16, 24))

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)
