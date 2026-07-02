"""Spectral Library — build/save/load a reference library and search against it.

The library file format (.speclib, JSON) is identical to SpectraView's, so
libraries built here open in SpectraView and vice versa.
"""
import numpy as np

from AnyQt.QtWidgets import QFileDialog, QListWidget, QSizePolicy

from Orange.data import ContinuousVariable, Domain, StringVariable, Table
from Orange.widgets import gui, settings
from Orange.widgets.widget import Input, Msg, Output, OWWidget

from ..core import load_library, save_library, search_library
from ..table_io import spectra_from_table, table_from_spectra

RANK_KEYS = ["correlation", "cosine", "sam", "euclid"]


class OWSpectralLibrary(OWWidget):
    name = "Spectral Library"
    description = ("Build a reference spectral library (.speclib, shared with "
                   "SpectraView) and search unknown spectra against it.")
    icon = "icons/library.svg"
    priority = 30
    keywords = ["library", "speclib", "search", "reference", "光譜庫", "比對"]

    class Inputs:
        spectra = Input("Spectra", Table)
        query = Input("Query", Table)

    class Outputs:
        hits = Output("Hits", Table)
        best = Output("Best Match", Table)
        library = Output("Library", Table)

    class Error(OWWidget.Error):
        bad_table = Msg("{}")
        file_error = Msg("{}")

    rank_by: int = settings.Setting(0)
    last_dir: str = settings.Setting("")
    want_main_area = False

    def __init__(self):
        super().__init__()
        self._library: list[dict] = []
        self._pending: list[dict] = []
        self._queries: list[dict] = []

        inbox = gui.widgetBox(self.controlArea, "Build")
        self.in_label = gui.label(inbox, self, "No input spectra.")
        gui.button(inbox, self, "Add input spectra to library",
                   callback=self.add_pending)

        lbox = gui.widgetBox(self.controlArea, "Library")
        self.listing = QListWidget()
        self.listing.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        lbox.layout().addWidget(self.listing)
        hb = gui.hBox(lbox)
        gui.button(hb, self, "Load…", callback=self.load_file)
        gui.button(hb, self, "Save…", callback=self.save_file)
        hb2 = gui.hBox(lbox)
        gui.button(hb2, self, "Remove selected", callback=self.remove_selected)
        gui.button(hb2, self, "Clear", callback=self.clear)

        sbox = gui.widgetBox(self.controlArea, "Search")
        gui.comboBox(sbox, self, "rank_by", label="Rank hits by:",
                     items=RANK_KEYS, orientation="horizontal",
                     callback=self.commit)
        self.result_label = gui.label(sbox, self, "Connect a Query input.")

    # ------------------------------------------------------------- inputs
    @Inputs.spectra
    def set_spectra(self, table):
        self.Error.bad_table.clear()
        self._pending = []
        if table is not None:
            try:
                self._pending = spectra_from_table(table)
            except ValueError as exc:
                self.Error.bad_table(str(exc))
        self.in_label.setText(f"{len(self._pending)} input spectrum(s) ready."
                              if self._pending else "No input spectra.")

    @Inputs.query
    def set_query(self, table):
        self.Error.bad_table.clear()
        self._queries = []
        if table is not None:
            try:
                self._queries = spectra_from_table(table)
            except ValueError as exc:
                self.Error.bad_table(str(exc))
        self.commit()

    # ------------------------------------------------------------ library
    def add_pending(self):
        for s in self._pending:
            self._library.append(s)
            self.listing.addItem(f"{s['name']}  ({s['x'].size} pts)")
        self.commit()

    def remove_selected(self):
        row = self.listing.currentRow()
        if 0 <= row < len(self._library):
            del self._library[row]
            self.listing.takeItem(row)
            self.commit()

    def clear(self):
        self._library.clear()
        self.listing.clear()
        self.commit()

    def load_file(self):
        self.Error.file_error.clear()
        path, _ = QFileDialog.getOpenFileName(
            self, "Load spectral library", self.last_dir,
            "Spectral library (*.speclib *.json);;All files (*.*)")
        if not path:
            return
        try:
            entries = load_library(path)
        except Exception as exc:  # noqa: BLE001
            self.Error.file_error(f"Could not load library: {exc}")
            return
        import os
        self.last_dir = os.path.dirname(path)
        for s in entries:
            self._library.append(s)
            self.listing.addItem(f"{s['name']}  ({s['x'].size} pts)")
        self.commit()

    def save_file(self):
        self.Error.file_error.clear()
        if not self._library:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save spectral library", self.last_dir,
            "Spectral library (*.speclib);;All files (*.*)")
        if not path:
            return
        if not path.lower().endswith(".speclib"):
            path += ".speclib"
        try:
            save_library(self._library, path)
        except Exception as exc:  # noqa: BLE001
            self.Error.file_error(f"Could not save library: {exc}")
            return
        import os
        self.last_dir = os.path.dirname(path)

    # ------------------------------------------------------------- search
    def commit(self):
        lib_table = (table_from_spectra(self._library) if self._library else None)
        self.Outputs.library.send(lib_table)

        if not self._queries or not self._library:
            self.result_label.setText(
                f"Library: {len(self._library)} entries. "
                + ("Connect a Query input." if not self._queries else ""))
            self.Outputs.hits.send(None)
            self.Outputs.best.send(None)
            return

        key = RANK_KEYS[self.rank_by]
        rows, metas, best_specs, lines = [], [], [], []
        for q in self._queries:
            hits = search_library(q, self._library, rank_by=key)
            for rank, h in enumerate(hits, start=1):
                rows.append([rank] + [h["scores"][k] for k in RANK_KEYS])
                metas.append([q["name"], h["name"]])
            if hits:
                best_specs.append(hits[0]["entry"])
                lines.append(f"{q['name']} → {hits[0]['name']} "
                             f"({key}={hits[0]['scores'][key]:.4f})")
        domain = Domain(
            [ContinuousVariable.make("rank")]
            + [ContinuousVariable.make(k) for k in RANK_KEYS],
            metas=[StringVariable.make("query"), StringVariable.make("match")])
        out = Table.from_numpy(domain, np.asarray(rows, float),
                               metas=np.asarray(metas, dtype=object))
        out.name = "library hits"
        self.result_label.setText("\n".join(lines))
        self.Outputs.hits.send(out)
        self.Outputs.best.send(table_from_spectra(best_specs) if best_specs else None)

    def send_report(self):
        self.report_items("Spectral Library", [
            ("Entries", len(self._library)),
            ("Ranked by", RANK_KEYS[self.rank_by]),
        ])


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWSpectralLibrary).run()
