"""Import Spectrum URL — fetch a spectrum from IRUG / SOPRANO / any URL."""
from AnyQt.QtCore import Qt
from AnyQt.QtWidgets import QLineEdit, QListWidget, QSizePolicy

import pyqtgraph as pg

from Orange.data import Table
from Orange.widgets import gui, settings
from Orange.widgets.widget import Msg, Output, OWWidget

from ..core import load_spectrum_url
from ..table_io import table_from_spectra

_COLORS = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3",
           "#937860", "#da8bc3", "#8c8c8c", "#ccb974", "#64b5cd"]


class OWImportSpectrumURL(OWWidget):
    name = "Import Spectrum URL"
    description = ("Fetch a spectrum from the web: IRUG id or page, "
                   "SOPRANO page, JCAMP-DX or CSV URL.")
    icon = "icons/importurl.svg"
    priority = 10
    keywords = ["irug", "soprano", "url", "download", "import", "光譜"]

    class Outputs:
        data = Output("Spectra", Table)

    class Error(OWWidget.Error):
        fetch_failed = Msg("{}")

    class Warning(OWWidget.Warning):
        no_overlap = Msg("Spectra share no common x-range; only the latest is output.")

    url: str = settings.Setting("")
    want_main_area = True

    def __init__(self):
        super().__init__()
        self._spectra: list[dict] = []

        box = gui.widgetBox(self.controlArea, "Source")
        gui.label(box, self, "IRUG id (e.g. 4119), IRUG/SOPRANO page URL,\n"
                             "or a direct JCAMP-DX / CSV URL:")
        self.url_edit = QLineEdit(self.url)
        self.url_edit.setPlaceholderText("4119  or  https://…")
        self.url_edit.returnPressed.connect(self.fetch)
        box.layout().addWidget(self.url_edit)
        gui.button(box, self, "Fetch && add", callback=self.fetch)

        lbox = gui.widgetBox(self.controlArea, "Fetched spectra")
        self.listing = QListWidget()
        self.listing.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        lbox.layout().addWidget(self.listing)
        gui.button(lbox, self, "Remove selected", callback=self.remove_selected)
        gui.button(lbox, self, "Clear all", callback=self.clear_all)

        self.plot = pg.PlotWidget(background="w")
        self.plot.showGrid(x=True, y=True, alpha=0.3)
        self.mainArea.layout().addWidget(self.plot)

    # ------------------------------------------------------------- actions
    def fetch(self):
        self.Error.fetch_failed.clear()
        self.url = self.url_edit.text().strip()
        if not self.url:
            return
        try:
            spec = load_spectrum_url(self.url)
        except Exception as exc:  # noqa: BLE001
            self.Error.fetch_failed(str(exc))
            return
        self._spectra.append(spec)
        self.listing.addItem(f"{spec['name']}  ({spec['x'].size} pts)")
        self._refresh()

    def remove_selected(self):
        row = self.listing.currentRow()
        if 0 <= row < len(self._spectra):
            del self._spectra[row]
            self.listing.takeItem(row)
            self._refresh()

    def clear_all(self):
        self._spectra.clear()
        self.listing.clear()
        self._refresh()

    # ------------------------------------------------------------ plumbing
    def _refresh(self):
        self.Warning.no_overlap.clear()
        self.plot.clear()
        for i, s in enumerate(self._spectra):
            pen = pg.mkPen(_COLORS[i % len(_COLORS)], width=1.5)
            self.plot.plot(s["x"], s["y"], pen=pen, name=s["name"])
        if self._spectra:
            x_label = self._spectra[-1].get("x_label", "x")
            self.plot.setLabel("bottom", x_label)
            self.plot.setLabel("left", "intensity")
            # IR/Raman convention: wavenumber decreases left-to-right.
            self.plot.getViewBox().invertX("cm-1" in x_label or "cm" in x_label)
        self.commit()

    def commit(self):
        if not self._spectra:
            self.Outputs.data.send(None)
            return
        try:
            table = table_from_spectra(self._spectra)
        except ValueError:
            self.Warning.no_overlap()
            table = table_from_spectra(self._spectra[-1:])
        self.Outputs.data.send(table)

    def send_report(self):
        self.report_items("Import Spectrum URL", [
            ("Spectra", len(self._spectra)),
            ("Names", ", ".join(s["name"] for s in self._spectra)),
        ])


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWImportSpectrumURL).run()
