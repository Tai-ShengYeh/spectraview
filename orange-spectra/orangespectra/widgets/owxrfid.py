"""XRF Element ID — label element emission lines on an XRF spectrum (keV)."""
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

from ..core import find_spectrum_peaks
from ..table_io import spectra_from_table
from ..xrf import identify_energy
from ._help import add_help

_COLORS = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3",
           "#937860", "#da8bc3", "#8c8c8c", "#ccb974", "#64b5cd"]
LINE_FILTERS = [("K + L lines", None), ("K lines only", {"Ka1", "Kb1"}),
                ("L lines only", {"La1", "Lb1"})]


class OWXRFElementID(OWWidget):
    name = "XRF Element ID"
    description = ("Find peaks in an XRF spectrum (energy axis in keV) and "
                   "label them with matching element emission lines.")
    icon = "icons/xrf.svg"
    priority = 70
    keywords = ["xrf", "element", "fluorescence", "kev", "元素", "螢光"]

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        elements = Output("Elements", Table)

    class Error(OWWidget.Error):
        bad_table = Msg("{}")

    class Warning(OWWidget.Warning):
        not_kev = Msg("X range {:.6g}–{:.6g} looks unlike keV (expected "
                      "~1–40); check the energy axis / calibration.")

    tolerance: float = settings.Setting(0.10)       # keV
    line_filter: int = settings.Setting(0)
    min_height: float = settings.Setting(3.0)       # % of full range
    min_prominence: float = settings.Setting(2.0)   # % of full range
    want_main_area = True

    def __init__(self):
        super().__init__()
        self._spectra = []

        add_help(self,
                 "接 XRF 能譜（x 軸＝keV）→ 自動尋峰並比對元素特徵譜線（Kα/Kβ/"
                 "Lα/Lβ，X-ray Data Booklet 參考值）→ 圖上標記元素，Elements 輸出"
                 "比對表。tolerance 是能量容差（keV）。\n"
                 "Label XRF peaks with matching element lines.", "xrf")

        box = gui.widgetBox(self.controlArea, "Identification")
        gui.doubleSpin(box, self, "tolerance", 0.01, 1.0, 0.01,
                       label="Energy tolerance (keV):", callback=self._recompute)
        gui.comboBox(box, self, "line_filter", label="Lines:",
                     items=[name for name, _ in LINE_FILTERS],
                     orientation="horizontal", callback=self._recompute)
        pbox = gui.widgetBox(self.controlArea, "Peak detection")
        gui.doubleSpin(pbox, self, "min_height", 0.0, 100.0, 0.5,
                       label="Min height (% of range):", callback=self._recompute)
        gui.doubleSpin(pbox, self, "min_prominence", 0.0, 100.0, 0.5,
                       label="Min prominence (% of range):", callback=self._recompute)
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
        self.Warning.not_kev.clear()
        self.ax.clear()
        if not self._spectra:
            self.info_label.setText("No data.")
            self.canvas.draw_idle()
            self.Outputs.elements.send(None)
            return

        flt = LINE_FILTERS[self.line_filter][1]
        rows, metas, found = [], [], []
        for i, s in enumerate(self._spectra):
            color = _COLORS[i % len(_COLORS)]
            x, y = np.asarray(s["x"], float), np.asarray(s["y"], float)
            if x.max() > 200 or x.max() < 0.5:
                self.Warning.not_kev(float(x.min()), float(x.max()))
            self.ax.plot(x, y, lw=1.0, color=color, label=s["name"])
            peaks = find_spectrum_peaks(
                x, y, min_height_frac=self.min_height / 100.0,
                min_prominence_frac=self.min_prominence / 100.0)
            for p in peaks:
                matches = identify_energy(p["center"], tol=self.tolerance,
                                          line_filter=flt)
                best = matches[0] if matches else None
                label = (f"{best['symbol']} {best['line_label']}" if best
                         else f"{p['center']:.2f}?")
                self.ax.annotate(label, (p["center"], p["height"]),
                                 textcoords="offset points", xytext=(0, 7),
                                 ha="center", fontsize=8, rotation=90,
                                 color=color if best else "#999999")
                self.ax.plot([p["center"]], [p["height"]], "v", ms=5,
                             color=color)
                if best:
                    found.append(best["symbol"])
                for m in matches:
                    rows.append([p["center"], m["energy"], m["delta"],
                                 p["height"]])
                    metas.append([s["name"], m["symbol"], m["name"],
                                  m["line_label"]])
                if not matches:
                    rows.append([p["center"], np.nan, np.nan, p["height"]])
                    metas.append([s["name"], "?", "unidentified", ""])

        self.ax.set_xlabel("energy (keV)")
        self.ax.set_ylabel("counts")
        if len(self._spectra) > 1:
            self.ax.legend(fontsize=8)
        self.figure.tight_layout()
        self.canvas.draw_idle()

        if not rows:
            self.info_label.setText("No peaks found — lower the thresholds.")
            self.Outputs.elements.send(None)
            return
        domain = Domain(
            [ContinuousVariable.make(k) for k in
             ("peak energy (keV)", "line energy (keV)", "delta (keV)", "height")],
            metas=[StringVariable.make(k) for k in
                   ("spectrum", "symbol", "element", "line")])
        out = Table.from_numpy(domain, np.asarray(rows, float),
                               metas=np.asarray(metas, dtype=object))
        out.name = "XRF elements"
        uniq = sorted(set(found))
        self.info_label.setText(
            f"{len(rows)} line matches.\nElements: {', '.join(uniq) or '—'}")
        self.Outputs.elements.send(out)

    def send_report(self):
        self.report_items("XRF Element ID", [
            ("Tolerance (keV)", self.tolerance),
            ("Lines", LINE_FILTERS[self.line_filter][0])])
        if self._spectra:
            self.report_plot(self.figure)


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWXRFElementID).run()
