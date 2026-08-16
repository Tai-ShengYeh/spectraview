"""Peak Finder — detect, measure and label peaks in spectra."""
import matplotlib
import numpy as np

matplotlib.use("Qt5Agg")
# The Qt backend must be imported *after* matplotlib.use(); do not let an
# import sorter hoist the lines below above it.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets import gui, settings
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from .. import mplfonts  # noqa: E402, F401  (CJK-capable preview fonts)
from ..core import find_spectrum_peaks
from ..table_io import spectra_from_table
from ._help import add_help

_COLORS = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3",
           "#937860", "#da8bc3", "#8c8c8c", "#ccb974", "#64b5cd"]


class OWPeakFinder(OWWidget):
    name = "Peak Finder"
    description = ("Find peaks in spectra, measure position / height / FWHM / "
                   "prominence, and label them on the plot.")
    icon = "icons/peaks.svg"
    priority = 60
    keywords = ["peak", "peaks", "fwhm", "prominence", "find", "波峰", "尋峰"]

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        peaks = Output("Peaks", Table)

    class Error(OWWidget.Error):
        bad_table = Msg("{}")

    min_height: float = settings.Setting(5.0)       # % of full range
    min_prominence: float = settings.Setting(3.0)   # % of full range
    min_distance: float = settings.Setting(0.0)     # x-units
    smooth_window: int = settings.Setting(0)        # SavGol points (0 = off)
    show_labels: bool = settings.Setting(True)
    want_main_area = True

    def __init__(self):
        super().__init__()
        self._spectra = []

        add_help(self,
                 "接光譜（Import Spectrum URL 或任何欄名＝波數的 Table）→ 調高度/"
                 "顯著度門檻 → 圖上自動標記波峰位置，Peaks 輸出表（position, height, "
                 "FWHM, prominence, area）可接 Data Table。\n"
                 "Find & label peaks; outputs a peak table.", "peaks")

        box = gui.widgetBox(self.controlArea, "Detection")
        gui.doubleSpin(box, self, "min_height", 0.0, 100.0, 0.5,
                       label="Min height (% of range):", callback=self._recompute)
        gui.doubleSpin(box, self, "min_prominence", 0.0, 100.0, 0.5,
                       label="Min prominence (% of range):", callback=self._recompute)
        gui.doubleSpin(box, self, "min_distance", 0.0, 1e6, 1.0,
                       label="Min distance (x units):", callback=self._recompute)
        gui.spin(box, self, "smooth_window", 0, 99, 2,
                 label="Smoothing window (0 = off):", callback=self._recompute)
        gui.checkBox(gui.widgetBox(self.controlArea, "Display"), self,
                     "show_labels", "Label peak positions on plot",
                     callback=self._recompute)
        self.info_label = gui.label(
            gui.widgetBox(self.controlArea, "Status"), self, "No data.")

        self.figure = Figure(figsize=(7, 5))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.mainArea.layout().addWidget(self.canvas)
        self.ax = self.figure.add_subplot(111)

    @Inputs.data
    def set_data(self, table):
        self.Error.bad_table.clear()
        self._spectra = []
        if table is not None:
            try:
                self._spectra = spectra_from_table(table)
            except ValueError as exc:
                self.Error.bad_table(str(exc))
        self._recompute()

    def _recompute(self):
        self.ax.clear()
        if not self._spectra:
            self.info_label.setText("No data.")
            self.canvas.draw_idle()
            self.Outputs.peaks.send(None)
            return

        rows, metas, n_peaks = [], [], 0
        for i, s in enumerate(self._spectra):
            color = _COLORS[i % len(_COLORS)]
            x, y = np.asarray(s["x"], float), np.asarray(s["y"], float)
            self.ax.plot(x, y, lw=1.0, color=color, label=s["name"])
            peaks = find_spectrum_peaks(
                x, y, min_height_frac=self.min_height / 100.0,
                min_prominence_frac=self.min_prominence / 100.0,
                min_distance=self.min_distance,
                smooth_window=self.smooth_window)
            n_peaks += len(peaks)
            for p in peaks:
                self.ax.plot([p["center"]], [p["height"]], "v", ms=6,
                             color=color)
                if self.show_labels:
                    self.ax.annotate(f"{p['center']:.6g}",
                                     (p["center"], p["height"]),
                                     textcoords="offset points", xytext=(0, 7),
                                     ha="center", fontsize=8, color=color,
                                     rotation=90 if len(peaks) > 8 else 0)
                rows.append([p["center"], p["height"], p["fwhm"],
                             p["prominence"], p["area"]])
                metas.append([s["name"]])

        self.ax.set_xlabel(self._spectra[0].get("x_label", "x"))
        self.ax.set_ylabel("intensity")
        if len(self._spectra) > 1:
            self.ax.legend(fontsize=8)
        self.figure.tight_layout()
        self.canvas.draw_idle()

        if not rows:
            self.info_label.setText("No peaks found — lower the thresholds.")
            self.Outputs.peaks.send(None)
            return
        domain = Domain(
            [ContinuousVariable.make(k) for k in
             ("position", "height", "fwhm", "prominence", "area")],
            metas=[StringVariable.make("spectrum")])
        out = Table.from_numpy(domain, np.asarray(rows, float),
                               metas=np.asarray(metas, dtype=object))
        out.name = "peaks"
        self.info_label.setText(
            f"{n_peaks} peaks in {len(self._spectra)} spectra.")
        self.Outputs.peaks.send(out)

    def send_report(self):
        self.report_items("Peak Finder", [
            ("Min height %", self.min_height),
            ("Min prominence %", self.min_prominence),
            ("Min distance", self.min_distance),
            ("Smoothing window", self.smooth_window or "off")])
        if self._spectra:
            self.report_plot(self.figure)


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWPeakFinder).run()
