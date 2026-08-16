"""Load Spectra Files — read a folder / .zip / files of spectra into one Table.

Reads JCAMP-DX (AFFN), two-column CSV/TSV, matrix CSV (both SpectraView
"combined export" layouts) and classic NetCDF (.cdf/.nc, e.g. chemometrics
datasets) — outputs everything as one merged Table (one spectrum per row).
"""
import os

import matplotlib
from AnyQt.QtWidgets import QFileDialog, QListWidget

matplotlib.use("Qt5Agg")
# The Qt backend must be imported *after* matplotlib.use(); do not let an
# import sorter hoist the lines below above it.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from Orange.data import Table
from Orange.widgets import gui, settings
from Orange.widgets.widget import Msg, Output, OWWidget

from .. import mplfonts  # noqa: E402, F401  (CJK-capable preview fonts)
from ..files import SPECTRA_EXTS, load_spectra_folder, load_spectra_path
from ..table_io import table_from_spectra
from ._help import add_help

_COLORS = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3",
           "#937860", "#da8bc3", "#8c8c8c", "#ccb974", "#64b5cd"]
_FILTER = ("Spectra files (*" + " *".join(sorted(SPECTRA_EXTS)) + " *.zip)"
           ";;All files (*)")


class OWLoadSpectraFiles(OWWidget):
    name = "Load Spectra Files"
    description = ("Load every spectrum in chosen files, a folder or a .zip "
                  "(JCAMP-DX, CSV, matrix CSV, NetCDF .cdf) into one Table.")
    icon = "icons/loadfiles.svg"
    priority = 12
    keywords = ["load", "file", "folder", "zip", "netcdf", "cdf", "import",
                "批次", "資料夾", "載入"]

    class Outputs:
        spectra = Output("Spectra", Table)

    class Error(OWWidget.Error):
        load_failed = Msg("{}")
        no_overlap = Msg("{}")

    class Warning(OWWidget.Warning):
        some_failed = Msg("{} file(s) could not be parsed (skipped).")

    class Information(OWWidget.Information):
        nothing = Msg("Add spectra files, a folder, or a .zip archive.")

    sources = settings.Setting([])          # list of file/folder paths
    recursive: bool = settings.Setting(True)
    want_main_area = True

    def __init__(self):
        super().__init__()
        self._spectra = []

        add_help(self,
                 "一鍵吃整包光譜：加入檔案、整個資料夾或 .zip（免解壓）→ 全部讀進來"
                 "疊圖，輸出合併 Table（每列一條光譜）。支援 JCAMP-DX、兩欄 CSV、"
                 "矩陣 CSV（欄名＝波長）與 NetCDF .cdf（如 applewine 等 chemometrics "
                 "資料集）。\nLoad a whole folder/zip of spectra into one Table.",
                 "loadfiles")

        box = gui.widgetBox(self.controlArea, "Sources")
        gui.button(box, self, "Add files…", callback=self._add_files)
        gui.button(box, self, "Add folder…", callback=self._add_folder)
        gui.button(box, self, "Add .zip…", callback=self._add_zip)
        gui.checkBox(box, self, "recursive", "Include subfolders",
                     callback=self._reload)
        self.src_list = QListWidget()
        self.src_list.setMaximumHeight(120)
        box.layout().addWidget(self.src_list)
        gui.button(box, self, "Remove selected", callback=self._remove)
        gui.button(box, self, "Clear all", callback=self._clear)

        self.info_label = gui.label(
            gui.widgetBox(self.controlArea, "Status"), self, "No sources.")

        self.figure = Figure(figsize=(7, 5))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.mainArea.layout().addWidget(self.canvas)
        self.ax = self.figure.add_subplot(111)

        self._sync_list()
        self._reload()

    # ------------------------------------------------------------- sources
    def _add_files(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Choose spectra files", "", _FILTER)
        self._extend(paths)

    def _add_folder(self):
        path = QFileDialog.getExistingDirectory(self, "Choose a folder")
        self._extend([path] if path else [])

    def _add_zip(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Choose a .zip archive", "", "Zip archives (*.zip)")
        self._extend([path] if path else [])

    def _extend(self, paths):
        added = [p for p in paths if p and p not in self.sources]
        if added:
            self.sources = self.sources + added
            self._sync_list()
            self._reload()

    def _remove(self):
        keep = [p for i, p in enumerate(self.sources)
                if not self.src_list.item(i).isSelected()]
        self.sources = keep
        self._sync_list()
        self._reload()

    def _clear(self):
        self.sources = []
        self._sync_list()
        self._reload()

    def _sync_list(self):
        self.src_list.clear()
        for p in self.sources:
            self.src_list.addItem(os.path.basename(p.rstrip("/\\")) or p)

    # --------------------------------------------------------------- load
    def _reload(self):
        self.Error.clear()
        self.Warning.some_failed.clear()
        self.Information.nothing.clear()
        self._spectra, failed = [], 0
        for p in self.sources:
            try:
                if os.path.isdir(p):
                    self._spectra.extend(
                        load_spectra_folder(p, recursive=self.recursive))
                else:
                    self._spectra.extend(load_spectra_path(p))
            except (ValueError, OSError):
                failed += 1
        if failed:
            self.Warning.some_failed(failed)
        self._replot()
        self._send()

    def _replot(self):
        self.ax.clear()
        for i, s in enumerate(self._spectra):
            self.ax.plot(s["x"], s["y"], lw=0.9,
                         color=_COLORS[i % len(_COLORS)],
                         label=s["name"] if len(self._spectra) <= 12 else None)
        if self._spectra:
            self.ax.set_xlabel(self._spectra[0].get("x_label", "x"))
            self.ax.set_ylabel("intensity")
            if len(self._spectra) <= 12:
                self.ax.legend(fontsize=8)
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def _send(self):
        if not self._spectra:
            if not self.sources:
                self.Information.nothing()
                self.info_label.setText("No sources.")
            else:
                self.info_label.setText("No spectra found in the sources.")
            self.Outputs.spectra.send(None)
            return
        try:
            out = table_from_spectra(self._spectra)
        except ValueError as exc:
            self.Error.no_overlap(str(exc))
            self.Outputs.spectra.send(None)
            return
        out.name = "loaded spectra"
        self.info_label.setText(
            f"{len(self._spectra)} spectra from {len(self.sources)} source(s)\n"
            f"{len(out.domain.attributes)} shared points")
        self.Outputs.spectra.send(out)

    def send_report(self):
        self.report_items("Load Spectra Files", [
            ("Sources", len(self.sources)),
            ("Spectra", len(self._spectra))])
        if self._spectra:
            self.report_plot(self.figure)


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWLoadSpectraFiles).run()
