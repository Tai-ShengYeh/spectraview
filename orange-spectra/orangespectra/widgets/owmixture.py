"""Mixture Analysis — decompose a mixed spectrum into reference components (NNLS)."""
import numpy as np

import pyqtgraph as pg
from AnyQt.QtWidgets import QTableWidget, QTableWidgetItem

from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets import gui, settings
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from ..core import make_spectrum, mixture_nnls
from ..table_io import spectra_from_table, table_from_spectra


class OWMixtureAnalysis(OWWidget):
    name = "Mixture Analysis"
    description = ("Estimate component proportions of a mixed spectrum from "
                   "pure reference spectra (non-negative least squares).")
    icon = "icons/mixture.svg"
    priority = 40
    keywords = ["mixture", "nnls", "unmixing", "components", "混合", "成分"]

    class Inputs:
        mixture = Input("Mixture", Table)
        references = Input("References", Table)

    class Outputs:
        composition = Output("Composition", Table)
        fit = Output("Fit", Table)

    class Error(OWWidget.Error):
        analysis_failed = Msg("{}")

    class Warning(OWWidget.Warning):
        multiple_rows = Msg("Mixture input has several rows; only the first is used.")

    fit_offset: bool = settings.Setting(True)
    want_main_area = True

    def __init__(self):
        super().__init__()
        self._mixture = None
        self._refs = []

        box = gui.widgetBox(self.controlArea, "Options")
        gui.checkBox(box, self, "fit_offset", "Fit a constant baseline offset",
                     callback=self.commit)
        self.info_label = gui.label(
            gui.widgetBox(self.controlArea, "Result"), self,
            "Connect a Mixture and References.")

        self.comp_table = QTableWidget(0, 3)
        self.comp_table.setHorizontalHeaderLabels(
            ["component", "coefficient", "fraction"])
        self.comp_table.horizontalHeader().setStretchLastSection(True)
        self.controlArea.layout().addWidget(self.comp_table)

        self.plot = pg.PlotWidget(background="w")
        self.plot.showGrid(x=True, y=True, alpha=0.3)
        self.plot.addLegend()
        self.mainArea.layout().addWidget(self.plot)

    @Inputs.mixture
    def set_mixture(self, table):
        self.Warning.multiple_rows.clear()
        self._mixture = None
        if table is not None:
            try:
                specs = spectra_from_table(table)
            except ValueError as exc:
                self.Error.analysis_failed(str(exc))
                specs = []
            if len(specs) > 1:
                self.Warning.multiple_rows()
            if specs:
                self._mixture = specs[0]
        self.commit()

    @Inputs.references
    def set_references(self, table):
        self._refs = []
        if table is not None:
            try:
                self._refs = spectra_from_table(table)
            except ValueError as exc:
                self.Error.analysis_failed(str(exc))
        self.commit()

    def commit(self):
        self.Error.analysis_failed.clear()
        self.plot.clear()
        self.comp_table.setRowCount(0)
        if self._mixture is None or not self._refs:
            self.info_label.setText("Connect a Mixture and References.")
            self.Outputs.composition.send(None)
            self.Outputs.fit.send(None)
            return
        try:
            res = mixture_nnls(self._mixture, self._refs,
                               fit_offset=self.fit_offset)
        except Exception as exc:  # noqa: BLE001
            self.Error.analysis_failed(str(exc))
            self.Outputs.composition.send(None)
            self.Outputs.fit.send(None)
            return

        # ---- composition table (widget + output)
        order = np.argsort(res["fractions"])[::-1]
        self.comp_table.setRowCount(len(order))
        for row, i in enumerate(order):
            for col, val in enumerate([res["names"][i], f"{res['coeffs'][i]:.4g}",
                                       f"{res['fractions'][i] * 100:.1f} %"]):
                self.comp_table.setItem(row, col, QTableWidgetItem(str(val)))
        self.info_label.setText(f"R² = {res['r_squared']:.4f}"
                                + (f"   offset = {res['offset']:.3g}"
                                   if self.fit_offset else ""))

        domain = Domain([ContinuousVariable.make("coefficient"),
                         ContinuousVariable.make("fraction")],
                        metas=[StringVariable.make("component")])
        comp = Table.from_numpy(
            domain,
            np.column_stack([res["coeffs"][order], res["fractions"][order]]),
            metas=np.asarray([[res["names"][i]] for i in order], dtype=object))
        comp.name = "mixture composition"
        comp.attributes["r_squared"] = res["r_squared"]
        self.Outputs.composition.send(comp)

        # ---- fit output: mixture (on fit grid), fit, residual as 3 rows
        gx = res["x"]
        mix_on_grid = np.interp(gx, self._mixture["x"], self._mixture["y"])
        x_label = self._mixture.get("x_label", "x")
        fit_table = table_from_spectra([
            make_spectrum(gx, mix_on_grid, name="mixture", x_label=x_label),
            make_spectrum(gx, res["fit"], name="NNLS fit", x_label=x_label),
            make_spectrum(gx, res["residual"], name="residual", x_label=x_label),
        ])
        self.Outputs.fit.send(fit_table)

        # ---- plot
        self.plot.plot(gx, mix_on_grid, pen=pg.mkPen("#4c72b0", width=1.6),
                       name="mixture")
        self.plot.plot(gx, res["fit"], pen=pg.mkPen("#c44e52", width=1.6,
                       style=pg.QtCore.Qt.DashLine), name="fit")
        self.plot.plot(gx, res["residual"], pen=pg.mkPen("#8c8c8c", width=1.0),
                       name="residual")
        self.plot.setLabel("bottom", x_label)
        self.plot.getViewBox().invertX("cm" in x_label)

    def send_report(self):
        self.report_items("Mixture Analysis", [
            ("References", len(self._refs)),
            ("Fit offset", self.fit_offset),
        ])


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWMixtureAnalysis).run()
