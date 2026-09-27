"""PLS Regression — quantitative PLS calibration of spectra, with VIP."""
import matplotlib
import numpy as np

matplotlib.use("Qt5Agg")
# The Qt backend must be imported *after* matplotlib.use(); do not let an
# import sorter hoist the lines below above it.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from Orange.data import (
    ContinuousVariable,
    Domain,
    StringVariable,
    Table,
)
from Orange.widgets import gui, settings
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from .. import mplfonts  # noqa: E402, F401  (CJK-capable preview fonts)
from ..core import pls_regression_fit
from ._help import add_help

_COLORS = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3",
           "#937860", "#da8bc3", "#8c8c8c", "#ccb974", "#64b5cd"]


def _numeric_wavelengths(names):
    """Parse every attribute name as a float, or return None if any fails."""
    out = []
    for n in names:
        try:
            out.append(float(n))
        except ValueError:
            return None
    return out


class OWPLSRegression(OWWidget):
    name = "PLS Regression"
    description = ("Quantitative partial least squares calibration: scores, "
                   "loadings, coefficients, VIP variable importance, RMSEC "
                   "and R² for a continuous target (e.g. concentration, %).")
    icon = "icons/plsvip.svg"
    priority = 81
    keywords = ["pls", "pls regression", "regression", "calibration",
                "chemometrics", "vip", "quantitative", "迴歸", "定量", "校正"]

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        scores = Output("Scores", Table, default=True)
        loadings = Output("Loadings", Table)
        vip = Output("VIP", Table)
        coefficients = Output("Coefficients", Table)
        predictions = Output("Predictions", Table)

    class Error(OWWidget.Error):
        no_target = Msg("Data needs at least one continuous target (class) "
                        "variable — use Select Columns to set one.")
        discrete_target = Msg("Target '{}' is categorical (discrete) — PLS "
                              "Regression needs a continuous target. Use "
                              "PLS-DA for classification instead.")
        compute_failed = Msg("{}")

    class Warning(OWWidget.Warning):
        components_clipped = Msg("Components reduced to {} (max = "
                                 "min(labelled samples - 1, features)).")
        rows_dropped = Msg("{} row(s) with a missing target were excluded "
                           "from fitting (predictions are still computed for "
                           "them).")

    n_components: int = settings.Setting(2)
    scale: bool = settings.Setting(False)
    want_main_area = True

    def __init__(self):
        super().__init__()
        self._data = None

        add_help(self,
                 "Connect a spectra Table with one or more continuous target "
                 "columns (e.g. % aspartame, acetic acid concentration) - use "
                 "Select Columns to set the target(s). PLS Regression fits "
                 "PLS2 (NIPALS); VIP > 1 is the usual 'important wavelength' "
                 "rule of thumb — check that high-VIP wavelengths match the "
                 "chemical bands you expect. Outputs: Scores / Loadings / "
                 "VIP / Coefficients / Predictions.",
                 "plsvip")

        box = gui.widgetBox(self.controlArea, "Model")
        gui.spin(box, self, "n_components", 1, 20, 1,
                 label="Components:", callback=self._recompute)
        gui.checkBox(box, self, "scale", "Autoscale (divide by std. dev.)",
                     callback=self._recompute)
        self.info_label = gui.label(
            gui.widgetBox(self.controlArea, "Status"), self, "No data.")

        self.figure = Figure(figsize=(6, 5))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.mainArea.layout().addWidget(self.canvas)
        self.ax = self.figure.add_subplot(111)

    @Inputs.data
    def set_data(self, table):
        self._data = table
        self._recompute()

    def _send_none(self):
        for out in (self.Outputs.scores, self.Outputs.loadings,
                    self.Outputs.vip, self.Outputs.coefficients,
                    self.Outputs.predictions):
            out.send(None)

    def _recompute(self):
        self.Error.clear()
        self.Warning.clear()
        self.ax.clear()
        if self._data is None or len(self._data) == 0:
            self.info_label.setText("No data.")
            self.canvas.draw_idle()
            self._send_none()
            return

        class_vars = self._data.domain.class_vars
        if not class_vars:
            self.Error.no_target()
            self.canvas.draw_idle()
            self._send_none()
            return
        bad = [v for v in class_vars if not v.is_continuous]
        if bad:
            self.Error.discrete_target(bad[0].name)
            self.canvas.draw_idle()
            self._send_none()
            return

        attrs = [a for a in self._data.domain.attributes if a.is_continuous]
        if not attrs:
            self.Error.compute_failed("No continuous attribute (spectral) columns.")
            self.canvas.draw_idle()
            self._send_none()
            return

        X_all = np.nan_to_num(np.asarray(self._data.transform(Domain(attrs)).X, float))
        Y_all = np.asarray(self._data.Y, float)
        if Y_all.ndim == 1:
            Y_all = Y_all[:, None]
        finite_mask = np.isfinite(Y_all).all(axis=1)
        n_fit = int(finite_mask.sum())
        n_dropped = len(finite_mask) - n_fit
        if n_dropped:
            self.Warning.rows_dropped(n_dropped)
        if n_fit < 2:
            self.Error.compute_failed(
                "Need at least 2 samples with a known target value.")
            self.canvas.draw_idle()
            self._send_none()
            return

        max_A = min(n_fit - 1, X_all.shape[1])
        A = self.n_components
        if A > max_A:
            A = max(1, max_A)
            self.Warning.components_clipped(A)

        try:
            res = pls_regression_fit(X_all[finite_mask], Y_all[finite_mask],
                                     A, self.scale)
        except Exception as exc:  # noqa: BLE001
            self.Error.compute_failed(str(exc))
            self.canvas.draw_idle()
            self._send_none()
            return

        y_pred_all = X_all @ res["coefficients"] + res["intercept"]

        self._draw(res, Y_all[finite_mask])
        A = res["n_components"]
        lines = [f"{n_fit} labelled samples "
                 f"({len(class_vars)} target(s)), {A} components."]
        for k, cv in enumerate(class_vars):
            n_vip = int((res["vip_per_target"][:, k] > 1).sum())
            lines.append(f"  {cv.name}: RMSEC={res['rmsec'][k]:.4g}  "
                         f"R²={res['r2'][k]:.4f}  VIP>1: {n_vip}")
        n_vip_combined = int((res["vip"] > 1).sum())
        lines.append(f"Combined VIP>1: {n_vip_combined} / {len(attrs)} variables.")
        self.info_label.setText("\n".join(lines))
        self._send_outputs(res, attrs, class_vars, X_all, y_pred_all, finite_mask)

    def _draw(self, res, Y_fit):
        T = res["scores"]
        self.ax.scatter(T[:, 0], T[:, 1] if T.shape[1] > 1 else np.zeros(T.shape[0]),
                        s=28, c=Y_fit[:, 0], cmap="viridis")
        self.ax.axhline(0, color="#cccccc", lw=0.7)
        self.ax.axvline(0, color="#cccccc", lw=0.7)
        self.ax.set_xlabel("t1")
        self.ax.set_ylabel("t2" if T.shape[1] > 1 else "")
        self.ax.set_title(f"PLS Regression scores — R² {res['r2'][0]:.3f}",
                          fontsize=10)
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def _send_outputs(self, res, attrs, class_vars, X_all, y_pred_all, finite_mask):
        A = res["n_components"]
        comp_names = [f"t{a + 1}" for a in range(A)]
        names = [a.name for a in attrs]
        wavelengths = _numeric_wavelengths(names)

        # Scores: fit-subset rows only (scores are only defined for samples
        # actually used in the fit).
        sdom = Domain([ContinuousVariable.make(n) for n in comp_names], class_vars)
        Y_fit = np.asarray(self._data.Y, float)
        if Y_fit.ndim == 1:
            Y_fit = Y_fit[:, None]
        scores = Table.from_numpy(sdom, res["scores"], Y_fit[finite_mask])
        scores.name = "PLS Regression scores"
        self.Outputs.scores.send(scores)

        # Loadings.
        ldom = Domain([ContinuousVariable.make(f"p{a + 1}") for a in range(A)],
                      metas=[StringVariable.make("variable")])
        loadings = Table.from_numpy(
            ldom, res["x_loadings"],
            metas=np.array([[n] for n in names], dtype=object))
        loadings.name = "PLS Regression loadings"
        self.Outputs.loadings.send(loadings)

        # VIP: kept in wavelength (input) order — a numeric 'wavelength' meta
        # (when attribute names parse as numbers) lets downstream widgets
        # such as Scatter Plot plot VIP against wavelength directly, and an
        # integer 'rank' meta (1 = highest combined VIP) gives the
        # sorted-by-importance view (as in PLS-DA's VIP output) on demand,
        # e.g. by sorting the Data Table on that column, without reordering
        # rows away from the spectral axis.
        vip_vars = [ContinuousVariable.make("VIP")] + \
            [ContinuousVariable.make(f"VIP_{cv.name}") for cv in class_vars]
        meta_vars = [StringVariable.make("variable")]
        if wavelengths is not None:
            meta_vars.append(ContinuousVariable.make("wavelength"))
        meta_vars.append(ContinuousVariable.make("rank"))
        vdom = Domain(vip_vars, metas=meta_vars)
        rank = np.empty(len(names), dtype=float)
        order = np.argsort(res["vip"])[::-1]
        rank[order] = np.arange(1, len(names) + 1)
        vip_x = np.column_stack([res["vip"], res["vip_per_target"]])
        meta_cols = [np.array([[n] for n in names], dtype=object)]
        if wavelengths is not None:
            meta_cols.append(np.array([[w] for w in wavelengths], dtype=float))
        meta_cols.append(rank[:, None])
        vip_metas = np.hstack(meta_cols)
        vip = Table.from_numpy(vdom, vip_x, metas=vip_metas)
        vip.name = "PLS Regression VIP"
        self.Outputs.vip.send(vip)

        # Coefficients: same feature-row layout as VIP; the intercept (one
        # scalar per target, not per-feature) is attached to each
        # coefficient variable's .attributes and echoed in the status box,
        # rather than given its own (meaningless) feature row.
        coef_vars = []
        for k, cv in enumerate(class_vars):
            v = ContinuousVariable.make(f"Coefficient_{cv.name}")
            v.attributes["intercept"] = float(res["intercept"][k])
            coef_vars.append(v)
        cdom = Domain(coef_vars, metas=meta_vars)
        coef = Table.from_numpy(cdom, res["coefficients"], metas=vip_metas)
        coef.name = "PLS Regression coefficients"
        self.Outputs.coefficients.send(coef)

        # Predictions: original data (all rows, including any with a
        # missing/unknown target) plus one predicted-value meta per target.
        pred_vars = [ContinuousVariable.make(f"Predicted {cv.name}")
                    for cv in class_vars]
        pdom = Domain(self._data.domain.attributes, self._data.domain.class_vars,
                      metas=self._data.domain.metas + tuple(pred_vars))
        pred = self._data.transform(pdom)
        with pred.unlocked(pred.metas):
            pred.metas[:, -len(pred_vars):] = y_pred_all
        pred.name = "PLS Regression predictions"
        self.Outputs.predictions.send(pred)

    def send_report(self):
        self.report_items("PLS Regression", [("Components", self.n_components),
                                              ("Autoscale", self.scale)])
        if self._data is not None:
            self.report_plot(self.figure)


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWPLSRegression).run(Table("housing"))
