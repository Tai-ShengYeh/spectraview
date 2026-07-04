"""Merge Spectra — overlay several spectra sources into one plot & Table.

Like SpectraView's multi-file overlay: connect several spectra Tables (e.g. from
multiple Import Spectrum URL widgets, or Orange's File), see them on one plot,
and get a single combined Table (each row a spectrum on a shared grid) for
Similarity / Library / PLS-DA.
"""
import numpy as np

import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from Orange.data import Table
from Orange.widgets import gui, settings
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from ..core import _snv
from ..table_io import spectra_from_table, table_from_spectra
from ._help import add_help

NORMS = [("none", "none"), ("max = 1", "max"), ("area = 1", "area"),
         ("SNV", "snv")]
_COLORS = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3",
           "#937860", "#da8bc3", "#8c8c8c", "#ccb974", "#64b5cd"]


def _normalize(y, mode):
    y = np.asarray(y, float)
    if mode == "max":
        m = np.nanmax(np.abs(y))
        return y / m if m else y
    if mode == "area":
        a = np.trapz(np.abs(y))
        return y / a if a else y
    if mode == "snv":
        return _snv(y)
    return y


class OWMergeSpectra(OWWidget):
    name = "Merge Spectra"
    description = ("Overlay several spectra sources on one plot and output a "
                  "single combined Table (each row a spectrum).")
    icon = "icons/merge.svg"
    priority = 15
    keywords = ["merge", "overlay", "combine", "stack", "疊圖", "合併"]

    class Inputs:
        data = Input("Spectra", Table, multiple=True)

    class Outputs:
        spectra = Output("Spectra", Table)

    class Error(OWWidget.Error):
        bad_table = Msg("{}")
        no_overlap = Msg("{}")

    normalization: int = settings.Setting(0)
    stack_offset: float = settings.Setting(0.0)
    want_main_area = True

    def __init__(self):
        super().__init__()
        self._inputs = {}                     # id -> Table

        add_help(self,
                 "把多個光譜來源（多個 Import Spectrum URL、File…）接到 Spectra 輸入"
                 "→ 疊在同一張圖，並輸出成一個合併 Table（每列一條光譜、共同波段），"
                 "可再接 Similarity / Library / PLS-DA。可選正規化與顯示堆疊位移。\n"
                 "Overlay multiple spectra into one plot & combined Table.", "merge")

        box = gui.widgetBox(self.controlArea, "Display / output")
        gui.comboBox(box, self, "normalization", label="Normalize each:",
                     items=[n for n, _ in NORMS], callback=self._recompute,
                     orientation="horizontal")
        gui.doubleSpin(box, self, "stack_offset", 0.0, 1e6, 0.1,
                       label="Stack offset (plot only):", callback=self._replot)
        self.info_label = gui.label(
            gui.widgetBox(self.controlArea, "Status"), self, "No inputs.")

        self.figure = Figure(figsize=(7, 5))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.mainArea.layout().addWidget(self.canvas)
        self.ax = self.figure.add_subplot(111)
        self._spectra = []

    @Inputs.data
    def set_data(self, table, id):
        if table is None:
            self._inputs.pop(id, None)
        else:
            self._inputs[id] = table

    def handleNewSignals(self):
        self._recompute()

    def _collect(self):
        spectra = []
        for table in self._inputs.values():
            spectra.extend(spectra_from_table(table))
        return spectra

    def _recompute(self):
        self.Error.bad_table.clear()
        self.Error.no_overlap.clear()
        try:
            self._spectra = self._collect()
        except ValueError as exc:
            self.Error.bad_table(str(exc))
            self._spectra = []
            self._replot()
            self.Outputs.spectra.send(None)
            return

        if not self._spectra:
            self.info_label.setText("No inputs.")
            self._replot()
            self.Outputs.spectra.send(None)
            return

        norm = NORMS[self.normalization][1]
        out_spectra = [dict(s, y=_normalize(s["y"], norm)) for s in self._spectra]
        try:
            out = table_from_spectra(out_spectra)
        except ValueError as exc:
            self.Error.no_overlap(str(exc))
            self._replot()
            self.Outputs.spectra.send(None)
            return
        out.name = "merged spectra"
        self.info_label.setText(
            f"{len(self._spectra)} spectra from {len(self._inputs)} inputs\n"
            f"combined on {len(out.domain.attributes)} shared points")
        self._replot()
        self.Outputs.spectra.send(out)

    def _replot(self):
        self.ax.clear()
        norm = NORMS[self.normalization][1]
        for i, s in enumerate(self._spectra):
            y = _normalize(s["y"], norm) + i * self.stack_offset
            self.ax.plot(s["x"], y, lw=1.0, color=_COLORS[i % len(_COLORS)],
                         label=s.get("name", f"spectrum {i}"))
        if self._spectra:
            self.ax.set_xlabel(self._spectra[0].get("x_label", "x"))
            self.ax.set_ylabel("intensity"
                               + (" (norm.)" if norm != "none" else ""))
            if len(self._spectra) <= 12:
                self.ax.legend(fontsize=8)
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def send_report(self):
        self.report_items("Merge Spectra", [
            ("Inputs", len(self._inputs)),
            ("Spectra", len(self._spectra)),
            ("Normalize", NORMS[self.normalization][0])])
        if self._spectra:
            self.report_plot(self.figure)


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWMergeSpectra).run()
