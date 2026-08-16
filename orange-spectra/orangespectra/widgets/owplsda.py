"""PLS-DA — partial least squares discriminant analysis of spectra."""
import matplotlib
import numpy as np

matplotlib.use("Qt5Agg")
# The Qt backend must be imported *after* matplotlib.use(); do not let an
# import sorter hoist the lines below above it.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from Orange.data import (
    ContinuousVariable,
    DiscreteVariable,
    Domain,
    StringVariable,
    Table,
)
from Orange.widgets import gui, settings
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from .. import mplfonts  # noqa: E402, F401  (CJK-capable preview fonts)
from ..core import plsda_fit
from ._help import add_help

_COLORS = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3",
           "#937860", "#da8bc3", "#8c8c8c", "#ccb974", "#64b5cd"]


class OWPLSDA(OWWidget):
    name = "PLS-DA"
    description = ("Partial least squares discriminant analysis: class-aware "
                   "projection of spectra with scores, loadings and VIP.")
    icon = "icons/plsda.svg"
    priority = 80
    keywords = ["pls", "plsda", "discriminant", "classification",
                "chemometrics", "vip", "判別"]

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        scores = Output("Scores", Table, default=True)
        loadings = Output("Loadings", Table)
        vip = Output("VIP", Table)
        predictions = Output("Predictions", Table)

    class Error(OWWidget.Error):
        no_class = Msg("Data needs a categorical target (class) variable — "
                       "use Select Columns to set one.")
        compute_failed = Msg("{}")

    n_components: int = settings.Setting(2)
    want_main_area = True

    def __init__(self):
        super().__init__()
        self._data = None

        add_help(self,
                 "接含類別（target/class）的光譜 Table（用 Select Columns 指定 "
                 "target）→ PLS-DA 把類別 one-hot 後做 PLS2（NIPALS）→ 分數圖依"
                 "類別上色。輸出 Scores / Loadings / VIP（>1 常視為重要波數）/ "
                 "Predictions。\nClass-aware projection; VIP > 1 ≈ important.",
                 "plsda")

        box = gui.widgetBox(self.controlArea, "Model")
        gui.spin(box, self, "n_components", 1, 20, 1,
                 label="Components:", callback=self._recompute)
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
                    self.Outputs.vip, self.Outputs.predictions):
            out.send(None)

    def _recompute(self):
        self.Error.clear()
        self.ax.clear()
        if self._data is None or len(self._data) == 0:
            self.info_label.setText("No data.")
            self.canvas.draw_idle()
            self._send_none()
            return
        cls = self._data.domain.class_var
        if cls is None or not cls.is_discrete:
            self.Error.no_class()
            self.canvas.draw_idle()
            self._send_none()
            return

        attrs = [a for a in self._data.domain.attributes if a.is_continuous]
        X = self._data.transform(Domain(attrs)).X
        labels = [cls.values[int(v)] for v in self._data.Y]
        try:
            res = plsda_fit(np.nan_to_num(np.asarray(X, float)), labels,
                            self.n_components)
        except Exception as exc:  # noqa: BLE001
            self.Error.compute_failed(str(exc))
            self.canvas.draw_idle()
            self._send_none()
            return

        self._draw(res, labels)
        A = res["n_components"]
        conf = "\n".join(
            "  " + " ".join(f"{v:4d}" for v in row) for row in res["confusion"])
        self.info_label.setText(
            f"{len(labels)} samples, {len(res['classes'])} classes, "
            f"{A} components.\nTraining accuracy: {res['accuracy']:.1%}\n"
            f"Confusion (rows = true):\n{conf}")
        self._send_outputs(res, labels, attrs)

    def _draw(self, res, labels):
        T = res["scores"]
        classes = res["classes"]
        xv = res["explained_x_variance"]
        for ci, c in enumerate(classes):
            mask = np.array([lab == c for lab in labels])
            if T.shape[1] > 1:
                self.ax.scatter(T[mask, 0], T[mask, 1], s=28,
                                color=_COLORS[ci % len(_COLORS)], label=c)
            else:
                self.ax.scatter(T[mask, 0], np.zeros(mask.sum()), s=28,
                                color=_COLORS[ci % len(_COLORS)], label=c)
        self.ax.axhline(0, color="#cccccc", lw=0.7)
        self.ax.axvline(0, color="#cccccc", lw=0.7)
        self.ax.set_xlabel(f"t1 ({xv[0]:.1%} X var)")
        self.ax.set_ylabel(f"t2 ({xv[1]:.1%} X var)" if len(xv) > 1 else "")
        self.ax.legend(fontsize=8)
        self.ax.set_title(f"PLS-DA scores — accuracy {res['accuracy']:.1%}",
                          fontsize=10)
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def _send_outputs(self, res, labels, attrs):
        A = res["n_components"]
        classes = res["classes"]
        comp_names = [f"t{a + 1}" for a in range(A)]
        cls_var = DiscreteVariable.make("class", values=tuple(classes))

        sdom = Domain([ContinuousVariable.make(n) for n in comp_names], cls_var)
        scores = Table.from_numpy(
            sdom, res["scores"],
            np.array([[classes.index(lab)] for lab in labels], float))
        scores.name = "PLS-DA scores"
        self.Outputs.scores.send(scores)

        ldom = Domain([ContinuousVariable.make(f"p{a + 1}") for a in range(A)],
                      metas=[StringVariable.make("variable")])
        loadings = Table.from_numpy(
            ldom, res["loadings"],
            metas=np.array([[a.name] for a in attrs], dtype=object))
        loadings.name = "PLS-DA loadings"
        self.Outputs.loadings.send(loadings)

        vdom = Domain([ContinuousVariable.make("VIP")],
                      metas=[StringVariable.make("variable")])
        order = np.argsort(res["vip"])[::-1]
        vip = Table.from_numpy(
            vdom, res["vip"][order, None],
            metas=np.array([[attrs[i].name] for i in order], dtype=object))
        vip.name = "PLS-DA VIP"
        self.Outputs.vip.send(vip)

        pred_var = DiscreteVariable.make("PLS-DA prediction",
                                         values=tuple(classes))
        pdom = Domain(self._data.domain.attributes, self._data.domain.class_var,
                      metas=self._data.domain.metas + (pred_var,))
        pred = self._data.transform(pdom)
        with pred.unlocked(pred.metas):
            pred.metas[:, -1] = [classes.index(p) for p in res["predicted"]]
        pred.name = "PLS-DA predictions"
        self.Outputs.predictions.send(pred)

    def send_report(self):
        self.report_items("PLS-DA", [("Components", self.n_components)])
        if self._data is not None:
            self.report_plot(self.figure)


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWPLSDA).run(Table("iris"))
