"""Aquagram — aquaphotomics radar plot over the 12 water bands (WAMACs).

Sources / method:
  * Aquaphotomics (Tsenkova) — the 12 water matrix coordinates (WAMACs).
  * https://nirpyresearch.com/aquagrams-python-matplotlib/ — SNV +
    across-sample normalization, then a 12-axis radar plot.
"""
import matplotlib
import numpy as np
from AnyQt.QtWidgets import QLineEdit

matplotlib.use("Qt5Agg")
# The Qt backend must be imported *after* matplotlib.use(); do not let an
# import sorter hoist the lines below above it.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets import gui, settings
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from .. import mplfonts  # noqa: E402, F401  (CJK-capable preview fonts)
from ..core import AQUAGRAM_NORMS, WAMACS, aquagram_coordinates
from ..table_io import spectra_from_table
from ._help import add_help

_COLORS = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3",
           "#937860", "#da8bc3", "#8c8c8c", "#ccb974", "#64b5cd"]


class OWAquagram(OWWidget):
    name = "Aquagram"
    description = ("Aquaphotomics aquagram: normalized absorbance at the 12 "
                   "water bands (WAMACs) drawn as a radar plot.")
    icon = "icons/aquagram.svg"
    priority = 50
    keywords = ["aquagram", "aquaphotomics", "water", "wamacs", "nir", "水光譜"]

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        coordinates = Output("Aquagram Coordinates", Table)

    class Error(OWWidget.Error):
        bad_table = Msg("{}")
        compute_failed = Msg("{}")

    class Warning(OWWidget.Warning):
        out_of_range = Msg("Some WAMACs lie outside the spectra's range; "
                           "edge values were used. Check the band list / units (nm).")

    normalization: int = settings.Setting(2)   # default "aquagram"
    bands_text: str = settings.Setting(", ".join(f"{v:g}" for v in WAMACS))
    want_main_area = True

    def __init__(self):
        super().__init__()
        self._spectra = []

        add_help(self,
                 "接 NIR 光譜（需涵蓋 ~1300–1600 nm 的水吸收區）→ 在水的 12 個特徵帶"
                 "（WAMACs）取正規化吸光度，畫 12 軸雷達圖。正規化選 aquagram 時 0＝組平均。"
                 "\nAquaphotomics radar over the 12 water bands.", "aquagram")

        box = gui.widgetBox(self.controlArea, "Normalization")
        gui.comboBox(box, self, "normalization", items=AQUAGRAM_NORMS,
                     callback=self._recompute)
        gui.label(box, self, "raw：原始吸光度｜snv：各譜 SNV｜"
                             "aquagram：SNV＋跨樣品標準化（0＝組平均）")

        bbox = gui.widgetBox(self.controlArea, "WAMACs bands (nm)")
        self.bands_edit = QLineEdit(self.bands_text)
        self.bands_edit.editingFinished.connect(self._bands_changed)
        bbox.layout().addWidget(self.bands_edit)
        gui.button(bbox, self, "Reset to standard 12", callback=self._reset_bands)

        self.info_label = gui.label(
            gui.widgetBox(self.controlArea, "Status"), self, "No data.")

        self.figure = Figure(figsize=(5, 5))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.mainArea.layout().addWidget(self.canvas)
        self.ax = self.figure.add_subplot(111, projection="polar")

    # ------------------------------------------------------------- input
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

    # ------------------------------------------------------------- bands
    def _bands_changed(self):
        self.bands_text = self.bands_edit.text()
        self._recompute()

    def _reset_bands(self):
        self.bands_text = ", ".join(f"{v:g}" for v in WAMACS)
        self.bands_edit.setText(self.bands_text)
        self._recompute()

    def _parse_bands(self):
        vals = []
        for tok in self.bands_text.replace(";", ",").split(","):
            tok = tok.strip()
            if tok:
                vals.append(float(tok))
        return vals

    # ---------------------------------------------------------- compute
    def _recompute(self):
        self.Error.compute_failed.clear()
        self.Warning.out_of_range.clear()
        self.ax.clear()
        if not self._spectra:
            self.info_label.setText("No data.")
            self.canvas.draw_idle()
            self.Outputs.coordinates.send(None)
            return
        try:
            bands = self._parse_bands()
            res = aquagram_coordinates(
                self._spectra, wamacs=bands,
                normalization=AQUAGRAM_NORMS[self.normalization])
        except Exception as exc:  # noqa: BLE001
            self.Error.compute_failed(str(exc))
            self.canvas.draw_idle()
            self.Outputs.coordinates.send(None)
            return
        if not res["covered"]:
            self.Warning.out_of_range()

        self._draw(res)
        self.info_label.setText(
            f"{len(res['names'])} spectra × {len(res['wamacs'])} WAMACs\n"
            f"normalization: {res['normalization']}")
        self._send(res)

    def _draw(self, res):
        bands = res["wamacs"]
        n = len(bands)
        angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
        closed = np.concatenate([angles, angles[:1]])
        for i, (name, row) in enumerate(zip(res["names"], res["values"])):
            vals = np.concatenate([row, row[:1]])
            color = _COLORS[i % len(_COLORS)]
            self.ax.plot(closed, vals, "-o", ms=3, lw=1.5, color=color,
                         label=name)
            self.ax.fill(closed, vals, color=color, alpha=0.06)
        self.ax.set_xticks(angles)
        self.ax.set_xticklabels([f"{b:g}" for b in bands], fontsize=8)
        self.ax.set_theta_offset(np.pi / 2)
        self.ax.set_theta_direction(-1)
        self.ax.set_title(f"Aquagram ({res['normalization']})", fontsize=11)
        if len(res["names"]) <= 12:
            self.ax.legend(loc="upper right", bbox_to_anchor=(1.28, 1.10),
                           fontsize=7)
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def _send(self, res):
        domain = Domain(
            [ContinuousVariable.make(f"{b:g}nm") for b in res["wamacs"]],
            metas=[StringVariable.make("name")])
        table = Table.from_numpy(
            domain, res["values"],
            metas=np.asarray([[n] for n in res["names"]], dtype=object))
        table.name = f"aquagram ({res['normalization']})"
        table.attributes["wamacs"] = list(res["wamacs"])
        table.attributes["normalization"] = res["normalization"]
        self.Outputs.coordinates.send(table)

    def send_report(self):
        self.report_items("Aquagram", [
            ("Normalization", AQUAGRAM_NORMS[self.normalization]),
            ("WAMACs", self.bands_text),
            ("Spectra", len(self._spectra)),
        ])
        self.report_plot(self.figure)


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWAquagram).run()
