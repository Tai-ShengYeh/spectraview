"""Microplate Reader — read a microtitre-plate photo into a per-well Table."""
import os

import numpy as np

from AnyQt.QtWidgets import QFileDialog

import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from Orange.data import (ContinuousVariable, DiscreteVariable, Domain,
                         StringVariable, Table)
from Orange.widgets import gui, settings
from Orange.widgets.widget import Msg, Output, OWWidget

from ..core import (PLATE_FORMATS, compute_areas, load_grayscale,
                    otsu_threshold, regular_grid, well_label)
from ._help import add_help

METRICS = [("mean intensity", "mean"), ("thresholded count", "count"),
           ("thresholded intensity", "intensity")]
THRESHOLDS = [("Otsu (auto)", "otsu"), ("manual", "manual"),
              ("per-well Otsu", "adaptive")]


class OWMicroplate(OWWidget):
    name = "Microplate Reader"
    description = ("Read a microtitre-plate photo (6–384 wells) into a table of "
                  "per-well measurements.")
    icon = "icons/microplate.svg"
    priority = 20
    keywords = ["microplate", "well", "elisa", "96", "384", "plate", "微孔盤", "盤"]

    class Outputs:
        wells = Output("Wells", Table, default=True)

    class Error(OWWidget.Error):
        load_failed = Msg("{}")

    class Information(OWWidget.Information):
        no_image = Msg("Choose a plate image to read.")

    image_path: str = settings.Setting("")
    plate_idx: int = settings.Setting(4)          # 96 (8×12)
    metric_idx: int = settings.Setting(0)          # mean
    threshold_idx: int = settings.Setting(0)       # Otsu
    manual_threshold: float = settings.Setting(0.5)
    margin: float = settings.Setting(0.02)
    fill: float = settings.Setting(0.7)            # round-well gutter
    want_main_area = True

    def __init__(self):
        super().__init__()
        self._gray = None

        add_help(self,
                 "選一張微孔盤照片 → 選盤規格（6/12/24/48/96/384 或不需設定，直接用"
                 "格數）→ 每孔取平均亮度（或閾值面積）→ 輸出 Wells 表（well_id、row、"
                 "col、value）。Row fill 調小可避免上下相鄰孔互相干擾。\n"
                 "Read a plate photo into a per-well table.", "microplate")

        fbox = gui.widgetBox(self.controlArea, "Image")
        gui.button(fbox, self, "Choose plate image…", callback=self._choose)
        self.file_label = gui.label(fbox, self, "(no file)")

        pbox = gui.widgetBox(self.controlArea, "Plate")
        gui.comboBox(pbox, self, "plate_idx",
                     items=list(PLATE_FORMATS), label="Format:",
                     callback=self._recompute, orientation="horizontal")
        gui.comboBox(pbox, self, "metric_idx", label="Measure:",
                     items=[n for n, _ in METRICS], callback=self._recompute,
                     orientation="horizontal")
        gui.comboBox(pbox, self, "threshold_idx", label="Threshold:",
                     items=[n for n, _ in THRESHOLDS], callback=self._recompute,
                     orientation="horizontal")
        gui.doubleSpin(pbox, self, "manual_threshold", 0.0, 1.0, 0.01,
                       label="Manual threshold:", callback=self._recompute)
        gui.doubleSpin(pbox, self, "margin", 0.0, 0.4, 0.01,
                       label="Border margin:", callback=self._recompute)
        gui.doubleSpin(pbox, self, "fill", 0.2, 1.0, 0.05,
                       label="Row fill (round wells ~0.7):", callback=self._recompute)
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

    def _resolve_threshold(self, metric):
        if metric == "mean":
            return None, False
        mode = THRESHOLDS[self.threshold_idx][1]
        if mode == "adaptive":
            return self.manual_threshold, True
        if mode == "otsu":
            return otsu_threshold(self._gray), False
        return self.manual_threshold, False

    def _recompute(self):
        self.ax.clear()
        if self._gray is None:
            self.info_label.setText("No image.")
            self.canvas.draw_idle()
            self.Outputs.wells.send(None)
            return

        n_rows, n_cols = PLATE_FORMATS[list(PLATE_FORMATS)[self.plate_idx]]
        grid = regular_grid(n_rows, n_cols, self._gray.shape,
                            margin=self.margin, fill=self.fill)
        metric = METRICS[self.metric_idx][1]
        thr, adaptive = self._resolve_threshold(metric)
        areas = compute_areas(self._gray, grid,
                              threshold=thr if thr is not None else 0.5,
                              metric=metric, adaptive=adaptive)

        self.ax.imshow(self._gray, cmap="gray", vmin=0, vmax=1)
        for r, c, x0, x1, y0, y1 in grid.cells():
            self.ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False,
                                        edgecolor="#dd8452", lw=0.8))
        self.ax.set_title(f"{list(PLATE_FORMATS)[self.plate_idx]} — "
                          f"{METRICS[self.metric_idx][0]}", fontsize=10)
        self.ax.set_xticks([]); self.ax.set_yticks([])
        self.figure.tight_layout()
        self.canvas.draw_idle()

        rows, metas = [], []
        for r in range(n_rows):
            for c in range(n_cols):
                rows.append([float(areas[r, c]), r + 1, c + 1])
                metas.append([well_label(r, c)])
        dom = Domain(
            [ContinuousVariable.make("value"), ContinuousVariable.make("row"),
             ContinuousVariable.make("col")],
            metas=[StringVariable.make("well")])
        out = Table.from_numpy(dom, np.asarray(rows, float),
                               metas=np.asarray(metas, dtype=object))
        out.name = "wells"
        self.info_label.setText(
            f"{n_rows}×{n_cols} = {n_rows * n_cols} wells\n"
            f"value range {np.nanmin(areas):.3g}–{np.nanmax(areas):.3g}")
        self.Outputs.wells.send(out)

    def send_report(self):
        self.report_items("Microplate Reader", [
            ("Image", os.path.basename(self.image_path) if self.image_path else "—"),
            ("Format", list(PLATE_FORMATS)[self.plate_idx]),
            ("Measure", METRICS[self.metric_idx][0])])
        if self._gray is not None:
            self.report_plot(self.figure)


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWMicroplate).run()
