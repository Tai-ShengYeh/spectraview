"""Coffee-Ring Reader — quantify a spotted-assay plate image into a cell table.

Reproduces coffee-ring-analyzer: per-cell thresholded area/intensity over a
grid, blank-column (column 0) normalization, replicate averaging. Output feeds
Dose-Response Fit for EC50 / LOD.
"""
import os

import numpy as np

from AnyQt.QtWidgets import QFileDialog

import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets import gui, settings
from Orange.widgets.widget import Msg, Output, OWWidget

from ..core import (DEFAULT_CONCENTRATIONS, DEFAULT_THRESHOLD, auto_detect_grid,
                    compute_areas, load_grayscale, normalize_blank,
                    otsu_threshold, regular_grid, well_label)
from ._help import add_help

METRICS = [("thresholded count", "count"), ("thresholded intensity", "intensity")]
THRESHOLDS = [("manual", "manual"), ("Otsu (auto)", "otsu"),
              ("per-cell Otsu", "adaptive")]


class OWCoffeeRing(OWWidget):
    name = "Coffee-Ring Reader"
    description = ("Quantify a spotted-assay (coffee-ring) plate photo: per-cell "
                  "area/intensity, blank-normalized, into a cell table.")
    icon = "icons/coffeering.svg"
    priority = 30
    keywords = ["coffee ring", "spot", "assay", "thrombin", "plate", "咖啡環", "點樣"]

    class Outputs:
        cells = Output("Cells", Table, default=True)

    class Error(OWWidget.Error):
        load_failed = Msg("{}")
        bad_conc = Msg("Concentrations: {}")

    class Warning(OWWidget.Warning):
        grid_failed = Msg("Auto-detect failed; using a regular grid instead.")

    class Information(OWWidget.Information):
        no_image = Msg("Choose a plate image to read.")

    image_path: str = settings.Setting("")
    n_rows: int = settings.Setting(3)
    n_cols: int = settings.Setting(8)
    metric_idx: int = settings.Setting(0)
    threshold_idx: int = settings.Setting(0)
    manual_threshold: float = settings.Setting(DEFAULT_THRESHOLD)
    auto_grid: bool = settings.Setting(False)
    normalize: bool = settings.Setting(True)
    conc_text: str = settings.Setting(", ".join(f"{v:g}"
                                                 for v in DEFAULT_CONCENTRATIONS))
    want_main_area = True

    def __init__(self):
        super().__init__()
        self._gray = None

        add_help(self,
                 "選點樣/咖啡環盤照片 → 設格數（列×欄，欄0＝blank）→ 選閾值與量測 → "
                 "每格算 above-threshold 面積/強度，除以同列 blank、平均重複列 → 輸出 "
                 "Cells 表（含濃度與比值），可直接接 Dose-Response Fit 求 EC50。\n"
                 "Quantify a coffee-ring plate into a cell table.", "coffeering")

        fbox = gui.widgetBox(self.controlArea, "Image")
        gui.button(fbox, self, "Choose plate image…", callback=self._choose)
        self.file_label = gui.label(fbox, self, "(no file)")

        gbox = gui.widgetBox(self.controlArea, "Grid")
        gui.spin(gbox, self, "n_rows", 1, 64, 1, label="Rows (replicates):",
                 callback=self._recompute)
        gui.spin(gbox, self, "n_cols", 2, 64, 1, label="Columns (col 0 = blank):",
                 callback=self._recompute)
        gui.checkBox(gbox, self, "auto_grid", "Auto-detect grid from spots",
                     callback=self._recompute)

        mbox = gui.widgetBox(self.controlArea, "Measurement")
        gui.comboBox(mbox, self, "metric_idx", label="Measure:",
                     items=[n for n, _ in METRICS], callback=self._recompute,
                     orientation="horizontal")
        gui.comboBox(mbox, self, "threshold_idx", label="Threshold:",
                     items=[n for n, _ in THRESHOLDS], callback=self._recompute,
                     orientation="horizontal")
        gui.doubleSpin(mbox, self, "manual_threshold", 0.0, 1.0, 0.01,
                       label="Manual threshold:", callback=self._recompute)

        nbox = gui.widgetBox(self.controlArea, "Normalization")
        gui.checkBox(nbox, self, "normalize",
                     "Divide by blank (col 0) & average rows",
                     callback=self._recompute)
        gui.lineEdit(nbox, self, "conc_text", label="Concentrations (cols 1…):",
                     callback=self._recompute)
        self.info_label = gui.label(
            gui.widgetBox(self.controlArea, "Status"), self, "No image.")

        self.figure = Figure(figsize=(7, 5))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.mainArea.layout().addWidget(self.canvas)
        self.ax = self.figure.add_subplot(111)

        if self.image_path and os.path.exists(self.image_path):
            self._load(self.image_path)
        else:
            self.Information.no_image()

    def _choose(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Choose plate image", self.image_path or "",
            "Images (*.png *.jpg *.jpeg *.tif *.tiff *.bmp)")
        if path:
            self._load(path)

    def _load(self, path):
        self.Error.load_failed.clear()
        self.Information.no_image.clear()
        try:
            self._gray = load_grayscale(path)
        except Exception as exc:                       # noqa: BLE001
            self.Error.load_failed(str(exc))
            self._gray = None
            return
        self.image_path = path
        self.file_label.setText(os.path.basename(path))
        self._recompute()

    def _parse_conc(self):
        vals = []
        for tok in self.conc_text.replace(";", ",").split(","):
            tok = tok.strip()
            if tok:
                vals.append(float(tok))
        return vals

    def _resolve_threshold(self):
        mode = THRESHOLDS[self.threshold_idx][1]
        if mode == "adaptive":
            return self.manual_threshold, True
        if mode == "otsu":
            return otsu_threshold(self._gray), False
        return self.manual_threshold, False

    def _recompute(self):
        self.Error.bad_conc.clear()
        self.Warning.grid_failed.clear()
        self.ax.clear()
        if self._gray is None:
            self.info_label.setText("No image.")
            self.canvas.draw_idle()
            self.Outputs.cells.send(None)
            return

        grid = None
        if self.auto_grid:
            grid = auto_detect_grid(self._gray, self.n_rows, self.n_cols)
            if grid is None:
                self.Warning.grid_failed()
        if grid is None:
            grid = regular_grid(self.n_rows, self.n_cols, self._gray.shape,
                                margin=0.0)

        thr, adaptive = self._resolve_threshold()
        metric = METRICS[self.metric_idx][1]
        areas = compute_areas(self._gray, grid, threshold=thr, metric=metric,
                              adaptive=adaptive)

        self.ax.imshow(self._gray, cmap="gray", vmin=0, vmax=1)
        for r, c, x0, x1, y0, y1 in grid.cells():
            self.ax.add_patch(Rectangle(
                (x0, y0), x1 - x0, y1 - y0, fill=False,
                edgecolor="#4c72b0" if c else "#c44e52", lw=1.0))
        self.ax.set_title("col 0 (red) = blank", fontsize=9)
        self.ax.set_xticks([]); self.ax.set_yticks([])
        self.figure.tight_layout()
        self.canvas.draw_idle()

        if self.normalize:
            self._output_normalized(areas, grid)
        else:
            self._output_raw(areas, grid)

    def _output_raw(self, areas, grid):
        rows, metas = [], []
        for r in range(grid.n_rows):
            for c in range(grid.n_cols):
                rows.append([float(areas[r, c]), r + 1, c + 1])
                metas.append([well_label(r, c)])
        dom = Domain([ContinuousVariable.make(n) for n in ("value", "row", "col")],
                     metas=[StringVariable.make("cell")])
        out = Table.from_numpy(dom, np.asarray(rows, float),
                               metas=np.asarray(metas, dtype=object))
        out.name = "cells (raw)"
        self.info_label.setText(
            f"{grid.n_rows}×{grid.n_cols} cells (raw, not normalized)\n"
            f"range {np.nanmin(areas):.3g}–{np.nanmax(areas):.3g}")
        self.Outputs.cells.send(out)

    def _output_normalized(self, areas, grid):
        norm = normalize_blank(areas, blank_col=0)
        means, ses, cols = norm["means"], norm["ses"], norm["columns"]
        conc = self._parse_conc()
        if len(conc) < len(cols):
            conc = conc + [float("nan")] * (len(cols) - len(conc))
        elif len(conc) > len(cols):
            conc = conc[:len(cols)]

        rows = [[conc[i], means[i], ses[i], cols[i] + 1]
                for i in range(len(cols))]
        dom = Domain([ContinuousVariable.make(n) for n in
                      ("concentration", "ratio", "se", "col")])
        out = Table.from_numpy(dom, np.asarray(rows, float))
        out.name = "cells (blank-normalized)"
        self.info_label.setText(
            f"{grid.n_rows} replicate rows × {grid.n_cols} cols\n"
            f"{len(cols)} treatment ratios (col 0 = blank)\n"
            "→ connect to Dose-Response Fit")
        self.Outputs.cells.send(out)

    def send_report(self):
        self.report_items("Coffee-Ring Reader", [
            ("Image", os.path.basename(self.image_path) if self.image_path else "—"),
            ("Grid", f"{self.n_rows}×{self.n_cols}"),
            ("Measure", METRICS[self.metric_idx][0]),
            ("Normalized", self.normalize)])
        if self._gray is not None:
            self.report_plot(self.figure)


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWCoffeeRing).run()
