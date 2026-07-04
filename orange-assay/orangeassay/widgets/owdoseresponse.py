"""Dose-Response Fit — 3PL/4PL logistic curve, EC50, LOD, R² from any Table."""
import numpy as np

import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets import gui, settings
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from ..core import fit_logistic
from ._help import add_help

MODELS = [("3-parameter (3PL)", "3pl"), ("4-parameter (4PL)", "4pl")]


def _pow10(v):
    """10**v, but never overflow/NaN-crash a poorly-constrained LOD/EC50."""
    if not np.isfinite(v):
        return float("nan")
    if v > 300:
        return float("inf")
    if v < -300:
        return 0.0
    return 10.0 ** v


class OWDoseResponse(OWWidget):
    name = "Dose-Response Fit"
    description = ("Fit a 3- or 4-parameter logistic dose-response curve to a "
                  "concentration/signal table; report EC50, LOD and R².")
    icon = "icons/doseresponse.svg"
    priority = 10
    keywords = ["dose", "response", "logistic", "ec50", "4pl", "3pl", "lod",
                "劑量", "檢量線", "curve"]

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        fit = Output("Fit summary", Table, default=True)
        curve = Output("Fitted curve", Table)

    class Error(OWWidget.Error):
        fit_failed = Msg("{}")

    class Warning(OWWidget.Warning):
        not_converged = Msg("Curve fit did not converge — showing initial guess.")

    conc_var = settings.ContextSetting(None)
    signal_var = settings.ContextSetting(None)
    model_idx: int = settings.Setting(0)
    log_x: bool = settings.Setting(True)
    settingsHandler = settings.DomainContextHandler()
    want_main_area = True

    def __init__(self):
        super().__init__()
        self._data = None

        add_help(self,
                 "接任何含『濃度』與『訊號』兩個數值欄的 Table（例如 Microplate / "
                 "Coffee-Ring Reader 的輸出，或比色法數據）→ 選欄位、選 3PL/4PL → "
                 "擬合出 EC50、LOD、R² 與曲線。log X 適合濃度跨數量級。\n"
                 "Fit a logistic dose-response curve; outputs EC50 / LOD / R².",
                 "doseresponse")

        box = gui.widgetBox(self.controlArea, "Columns")
        self.var_model = DomainModel(valid_types=ContinuousVariable)
        gui.comboBox(box, self, "conc_var", label="Concentration:",
                     model=self.var_model, callback=self._recompute,
                     orientation="horizontal")
        gui.comboBox(box, self, "signal_var", label="Signal:",
                     model=self.var_model, callback=self._recompute,
                     orientation="horizontal")
        gui.checkBox(box, self, "log_x", "Fit on log10(concentration)",
                     callback=self._recompute)

        mbox = gui.widgetBox(self.controlArea, "Model")
        gui.comboBox(mbox, self, "model_idx",
                     items=[name for name, _ in MODELS],
                     callback=self._recompute)
        self.info_label = gui.label(
            gui.widgetBox(self.controlArea, "Result"), self, "No data.")

        self.figure = Figure(figsize=(6, 5))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.mainArea.layout().addWidget(self.canvas)
        self.ax = self.figure.add_subplot(111)

    @Inputs.data
    def set_data(self, table):
        self.closeContext()
        self._data = table
        self.conc_var = self.signal_var = None
        if table is not None:
            self.var_model.set_domain(table.domain)
            cont = [v for v in table.domain.variables + table.domain.metas
                    if v.is_continuous]
            if len(cont) >= 2:
                self.conc_var, self.signal_var = cont[0], cont[1]
            self.openContext(table.domain)
        else:
            self.var_model.set_domain(None)
        self._recompute()

    def _column(self, var):
        return self._data.get_column(var).astype(float)

    def _recompute(self):
        self.Error.fit_failed.clear()
        self.Warning.not_converged.clear()
        self.ax.clear()
        if self._data is None or self.conc_var is None or self.signal_var is None:
            self.info_label.setText("No data.")
            self.canvas.draw_idle()
            self.Outputs.fit.send(None)
            self.Outputs.curve.send(None)
            return

        conc = self._column(self.conc_var)
        signal = self._column(self.signal_var)
        x = np.log10(np.where(conc > 0, conc, np.nan)) if self.log_x else conc
        model = MODELS[self.model_idx][1]
        try:
            fit = fit_logistic(x, signal, model=model)
        except ValueError as exc:
            self.Error.fit_failed(str(exc))
            self.canvas.draw_idle()
            self.Outputs.fit.send(None)
            self.Outputs.curve.send(None)
            return
        if not fit.success:
            self.Warning.not_converged()

        self._draw(fit, conc)
        ec50 = _pow10(fit.ec50) if self.log_x else fit.ec50
        lod = _pow10(fit.lod) if self.log_x else fit.lod
        pstr = ", ".join(f"{k}={v:.3g}" for k, v in fit.params.items())
        xlbl = "log10 " + self.conc_var.name if self.log_x else self.conc_var.name
        self.info_label.setText(
            f"{model.upper()}  R²={fit.r2:.4f}\nEC50={ec50:.4g}"
            f"  LOD={lod:.4g}\n{pstr}")
        self._send_outputs(fit, ec50, lod, xlbl)

    def _draw(self, fit, conc):
        x, y = fit.x, fit.y
        self.ax.scatter(x, y, s=32, color="#c44e52", zorder=3, label="data")
        gx = np.linspace(x.min(), x.max(), 200)
        self.ax.plot(gx, fit.predict(gx), color="#4c72b0", lw=2, label="fit")
        self.ax.axvline(fit.ec50, color="#55a868", ls="--", lw=1,
                        label=f"EC50 (x0={fit.ec50:.3g})")
        self.ax.set_xlabel(("log10 " if self.log_x else "")
                           + self.conc_var.name)
        self.ax.set_ylabel(self.signal_var.name)
        self.ax.set_title(f"{MODELS[self.model_idx][1].upper()}  R²={fit.r2:.4f}",
                          fontsize=10)
        self.ax.legend(fontsize=8)
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def _send_outputs(self, fit, ec50, lod, xlbl):
        keys = list(fit.params) + ["EC50", "LOD", "R2"]
        vals = list(fit.params.values()) + [ec50, lod, fit.r2]
        dom = Domain([ContinuousVariable.make(k) for k in keys],
                     metas=[StringVariable.make("model"),
                            StringVariable.make("x")])
        out = Table.from_numpy(
            dom, np.asarray([vals], float),
            metas=np.array([[fit.model, xlbl]], dtype=object))
        out.name = "dose-response fit"
        self.Outputs.fit.send(out)

        gx = np.linspace(fit.x.min(), fit.x.max(), 200)
        cdom = Domain([ContinuousVariable.make(xlbl),
                       ContinuousVariable.make("fit")])
        curve = Table.from_numpy(cdom, np.column_stack([gx, fit.predict(gx)]))
        curve.name = "fitted curve"
        self.Outputs.curve.send(curve)

    def send_report(self):
        self.report_items("Dose-Response Fit", [
            ("Model", MODELS[self.model_idx][0]),
            ("log X", self.log_x)])
        if self._data is not None:
            self.report_plot(self.figure)


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWDoseResponse).run()
