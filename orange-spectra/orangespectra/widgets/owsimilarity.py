"""Spectra Similarity — pairwise similarity scores between two sets of spectra."""
import numpy as np

from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets import gui, settings
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from ..core import similarity_scores
from ..table_io import spectra_from_table

RANK_KEYS = ["correlation", "cosine", "sam", "euclid"]


class OWSpectraSimilarity(OWWidget):
    name = "Spectra Similarity"
    description = ("Correlation, cosine, spectral angle (SAM) and Euclidean "
                   "similarity between spectra.")
    icon = "icons/similarity.svg"
    priority = 20
    keywords = ["similarity", "correlation", "cosine", "sam", "比對", "相似度"]

    class Inputs:
        data = Input("Data", Table)
        references = Input("References", Table)

    class Outputs:
        scores = Output("Scores", Table)

    class Error(OWWidget.Error):
        bad_table = Msg("{}")

    rank_by: int = settings.Setting(0)
    want_main_area = False
    resizing_enabled = False

    def __init__(self):
        super().__init__()
        self._data = None
        self._refs = None
        box = gui.widgetBox(self.controlArea, "Options")
        gui.comboBox(box, self, "rank_by", label="Sort scores by:",
                     items=[k for k in RANK_KEYS], orientation="horizontal",
                     callback=self.commit)
        gui.label(box, self, "If no References input is connected,\n"
                             "all pairs within Data are compared.")
        self.info_label = gui.label(gui.widgetBox(self.controlArea, "Result"),
                                    self, "No data.")

    @Inputs.data
    def set_data(self, table):
        self._data = table
        self.commit()

    @Inputs.references
    def set_references(self, table):
        self._refs = table
        self.commit()

    def commit(self):
        self.Error.bad_table.clear()
        if self._data is None:
            self.info_label.setText("No data.")
            self.Outputs.scores.send(None)
            return
        try:
            queries = spectra_from_table(self._data)
            refs = (spectra_from_table(self._refs)
                    if self._refs is not None else queries)
        except ValueError as exc:
            self.Error.bad_table(str(exc))
            self.Outputs.scores.send(None)
            return

        within = self._refs is None
        rows, metas = [], []
        for i, q in enumerate(queries):
            for j, r in enumerate(refs):
                if within and j <= i:      # skip self- and duplicate pairs
                    continue
                s = similarity_scores(q["x"], q["y"], r["x"], r["y"])
                rows.append([s[k] for k in RANK_KEYS])
                metas.append([q["name"], r["name"]])
        if not rows:
            self.info_label.setText("Need at least 2 spectra.")
            self.Outputs.scores.send(None)
            return

        key = RANK_KEYS[self.rank_by]
        order = np.argsort([r[self.rank_by] for r in rows])
        if key in ("correlation", "cosine"):
            order = order[::-1]
        rows = [rows[i] for i in order]
        metas = [metas[i] for i in order]

        domain = Domain(
            [ContinuousVariable.make(k) for k in RANK_KEYS],
            metas=[StringVariable.make("query"), StringVariable.make("reference")])
        out = Table.from_numpy(domain, np.asarray(rows, float),
                               metas=np.asarray(metas, dtype=object))
        out.name = "similarity scores"
        best = metas[0]
        self.info_label.setText(
            f"{len(rows)} pairs scored.\nBest ({key}): "
            f"{best[0]} ↔ {best[1]}  ({rows[0][self.rank_by]:.4f})")
        self.Outputs.scores.send(out)

    def send_report(self):
        self.report_items("Spectra Similarity",
                          [("Ranked by", RANK_KEYS[self.rank_by])])


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWSpectraSimilarity).run()
