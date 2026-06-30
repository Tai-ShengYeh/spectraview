"""The application main window."""
from __future__ import annotations

import os
import traceback

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets

from .. import (__app_name__, __version__, analysis, axes, calibration, cos2d,
                eem, processing, xrf)
from ..demo import demo_cos_series, demo_eem, demo_eem_stack, load_demo_set
from ..library import SpectralLibrary
from ..formats import (OPEN_FILTER, MissingDependency, load_any, load_online,
                       load_soprano_url, save_combined_csv, save_csv, save_jcamp,
                       save_json)
from ..spectrum import Spectrum, SpectrumSet, X_UNIT_LABELS, Y_UNIT_LABELS
from .calibration_view import CalibrationDialog, CalibrationWindow
from .dialogs import FormDialog, TableDialog
from .mapwindow import EEMWindow, MapWindow, ParafacWindow
from .plotview import PlotView

_X_TARGETS = [("Wavelength (nm)", "nm"), ("Wavenumber (cm⁻¹)", "cm-1"),
              ("Wavelength (µm)", "um"), ("Raman shift (cm⁻¹)", "raman_cm-1"),
              ("Energy (eV)", "eV"), ("Energy (keV)", "keV"),
              ("Frequency (THz)", "THz")]
_Y_TARGETS = [("Transmittance", "transmittance"), ("Transmittance %", "%T"),
              ("Absorbance", "absorbance"), ("Reflectance", "reflectance"),
              ("Kubelka-Munk", "KM"), ("log(1/R)", "log1R")]


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"{__app_name__} {__version__}")
        self.resize(1200, 760)
        self.setAcceptDrops(True)

        self.document = SpectrumSet()
        self._undo: list[list[Spectrum]] = []
        self._redo: list[list[Spectrum]] = []
        self._dark = False
        self._last_dir = ""   # remembered folder for open/save dialogs
        self._dialogs: list = []   # keep non-modal result dialogs alive
        self.library = SpectralLibrary()   # reference spectral library
        self.eems: list = []       # opened EEMs (for PARAFAC stacks)

        self.plotview = PlotView(self.document)
        self.plotview.cursorMoved.connect(self._on_cursor)

        self._build_spectrum_panel()
        self._build_central()
        self._build_actions()
        self._build_menus()
        self._build_toolbar()

        self.status = self.statusBar()
        self._cursor_label = QtWidgets.QLabel("  ")
        self.status.addPermanentWidget(self._cursor_label)
        self.status.showMessage("Open a file, drag-and-drop, or use File ▸ Load demo spectra.")

    # ===================================================== layout
    def _build_spectrum_panel(self) -> None:
        self.table = QtWidgets.QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["", "Colour", "Spectrum"])
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setColumnWidth(0, 28)
        self.table.setColumnWidth(1, 56)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.DoubleClicked
            | QtWidgets.QAbstractItemView.EditTrigger.EditKeyPressed
        )
        self.table.itemChanged.connect(self._on_item_changed)
        self.table.cellDoubleClicked.connect(self._on_cell_double_clicked)

    def _build_central(self) -> None:
        splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        left = QtWidgets.QWidget()
        lv = QtWidgets.QVBoxLayout(left)
        lv.setContentsMargins(4, 4, 4, 4)
        lv.addWidget(QtWidgets.QLabel("Spectra"))
        lv.addWidget(self.table)
        splitter.addWidget(left)
        splitter.addWidget(self.plotview)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([300, 900])
        self.setCentralWidget(splitter)
        self.left_panel = left   # the spectra-list sidebar (hideable via View menu)

    def _set_list_visible(self, visible: bool) -> None:
        """Show/hide the spectra-list sidebar so it stops taking up space."""
        self.left_panel.setVisible(visible)

    def _build_actions(self) -> None:
        SP = QtWidgets.QStyle.StandardPixmap

        def icon(sp):
            return self.style().standardIcon(sp)
        self.act_open = QtGui.QAction(icon(SP.SP_DialogOpenButton), "&Open…", self,
                                      shortcut="Ctrl+O", triggered=self.open_files)
        self.act_import_url = QtGui.QAction("Import from &URL / IRUG…", self,
                                            triggered=self.import_from_url)
        self.act_open_soprano = QtGui.QAction("Open SOPRANO URL…", self,
                                              triggered=self.open_soprano_url)
        self.act_demo = QtGui.QAction("Load &demo spectra", self, triggered=self.load_demo)
        self.act_save = QtGui.QAction(icon(SP.SP_DialogSaveButton), "&Save spectrum…",
                                      self, shortcut="Ctrl+S", triggered=self.save_spectrum)
        self.act_export = QtGui.QAction("&Export plot image…", self,
                                        triggered=self.export_image)
        self.act_export_data = QtGui.QAction("Export &combined data…", self,
                                             shortcut="Ctrl+Shift+S",
                                             triggered=self.export_combined)
        self.act_quit = QtGui.QAction("&Quit", self, shortcut="Ctrl+Q",
                                      triggered=self.close)
        self.act_undo = QtGui.QAction(icon(SP.SP_ArrowBack), "&Undo", self,
                                      shortcut="Ctrl+Z", triggered=self.undo)
        self.act_redo = QtGui.QAction(icon(SP.SP_ArrowForward), "&Redo", self,
                                      shortcut="Ctrl+Y", triggered=self.redo)
        self.act_remove = QtGui.QAction("Remove selected", self, shortcut="Delete",
                                        triggered=self.remove_selected)
        self.act_clear = QtGui.QAction("Remove all", self, triggered=self.remove_all)
        self.act_autoscale = QtGui.QAction(icon(SP.SP_FileDialogContentsView),
                                           "Auto-scale", self, shortcut="Ctrl+0",
                                           triggered=self.plotview.autoscale)
        self.act_toggle_list = QtGui.QAction("Show spectra list", self,
                                             checkable=True, checked=True,
                                             shortcut="F9",
                                             triggered=self._set_list_visible)

    def _build_menus(self) -> None:
        mb = self.menuBar()
        m_file = mb.addMenu("&File")
        m_file.addActions([self.act_open, self.act_import_url, self.act_open_soprano,
                           self.act_demo])
        m_file.addSeparator()
        m_file.addActions([self.act_save, self.act_export_data, self.act_export])
        m_file.addSeparator()
        m_file.addAction(self.act_quit)

        m_edit = mb.addMenu("&Edit")
        m_edit.addActions([self.act_undo, self.act_redo])
        m_edit.addSeparator()
        m_edit.addActions([self.act_remove, self.act_clear])

        m_view = mb.addMenu("&View")
        mx = m_view.addMenu("X axis ▸ convert to")
        for label, unit in _X_TARGETS:
            mx.addAction(QtGui.QAction(label, self,
                         triggered=lambda _=False, u=unit: self.convert_x(u)))
        my = m_view.addMenu("Y axis ▸ convert to")
        for label, unit in _Y_TARGETS:
            my.addAction(QtGui.QAction(label, self,
                         triggered=lambda _=False, u=unit: self.convert_y(u)))
        m_view.addSeparator()
        self.act_flip = QtGui.QAction("Flip X axis", self, checkable=True,
                                      triggered=self.plotview.set_flip_x)
        self.act_grid = QtGui.QAction("Grid", self, checkable=True, checked=True,
                                      triggered=self.plotview.set_grid)
        self.act_logy = QtGui.QAction("Log Y", self, checkable=True,
                                      triggered=self.plotview.set_log_y)
        self.act_dark = QtGui.QAction("Dark background", self, checkable=True,
                                      triggered=self.set_dark)
        m_view.addActions([self.act_flip, self.act_grid, self.act_logy, self.act_dark])
        m_view.addAction(QtGui.QAction("Stack / offset…", self, triggered=self.set_stack))
        m_view.addAction(self.act_autoscale)
        m_view.addSeparator()
        m_view.addAction(self.act_toggle_list)

        m_proc = mb.addMenu("&Process")
        ms = m_proc.addMenu("Smoothing")
        ms.addAction(QtGui.QAction("Savitzky-Golay…", self, triggered=self.smooth_sg))
        ms.addAction(QtGui.QAction("Moving average…", self, triggered=self.smooth_ma))
        mb_ = m_proc.addMenu("Baseline")
        mb_.addAction(QtGui.QAction("Rubberband", self, triggered=self.baseline_rb))
        mb_.addAction(QtGui.QAction("Polynomial…", self, triggered=self.baseline_poly))
        mb_.addAction(QtGui.QAction("Asymmetric Least Squares…", self,
                                    triggered=self.baseline_als))
        mb_.addAction(QtGui.QAction("airPLS…", self, triggered=self.baseline_airpls))
        m_proc.addAction(QtGui.QAction("Derivative…", self, triggered=self.derivative))
        mn = m_proc.addMenu("Normalize")
        for label, meth in [("Peak = 1 (max)", "max"), ("Min–Max (0..1)", "minmax"),
                            ("Area = 1", "area"), ("Vector (L2)", "vector")]:
            mn.addAction(QtGui.QAction(label, self,
                         triggered=lambda _=False, m=meth: self.normalize(m)))
        mn.addAction(QtGui.QAction("At value x₀…", self, triggered=self.normalize_value))
        msc = m_proc.addMenu("Scatter correction")
        msc.addAction(QtGui.QAction("SNV", self, triggered=self.snv))
        msc.addAction(QtGui.QAction("MSC (vs mean)", self, triggered=self.msc))
        msc.addAction(QtGui.QAction("Detrend…", self, triggered=self.detrend))
        ma = m_proc.addMenu("Arithmetic")
        for label, op in [("Subtract (A − B)", "sub"), ("Add (A + B)", "add"),
                          ("Multiply (A × B)", "mul"), ("Divide (A ÷ B)", "div")]:
            ma.addAction(QtGui.QAction(label, self,
                         triggered=lambda _=False, o=op: self.arithmetic(o)))
        ma.addAction(QtGui.QAction("Average selected (mean)", self,
                                   triggered=lambda: self.average("mean")))
        ma.addAction(QtGui.QAction("Average selected (median)", self,
                                   triggered=lambda: self.average("median")))
        m_proc.addSeparator()
        m_proc.addAction(QtGui.QAction("Cut range…", self, triggered=self.cut_range))
        m_proc.addAction(QtGui.QAction("Interpolate / resample…", self,
                                       triggered=self.interpolate))

        m_an = mb.addMenu("&Analyze")
        m_an.addAction(QtGui.QAction("Find peaks…", self, triggered=self.find_peaks))
        m_an.addAction(QtGui.QAction("Peak fit / deconvolution…", self,
                                     triggered=self.fit_peaks))
        m_an.addAction(QtGui.QAction("Integrate range…", self,
                                     triggered=self.integrate_range))
        m_an.addSeparator()
        m_an.addAction(QtGui.QAction("Mixture analysis (NNLS)…", self,
                                     triggered=self.mixture_analysis))
        m_an.addAction(QtGui.QAction("Identify XRF elements…", self,
                                     triggered=self.identify_xrf))
        m_an.addAction(QtGui.QAction("Calibration curve (檢量線)…", self,
                                     triggered=self.calibration_curve))
        m_an.addSeparator()
        m_an.addAction(QtGui.QAction("2D correlation (2D-COS / 2T2D)…", self,
                                     triggered=self.cos2d_analysis))
        m_eem = m_an.addMenu("Fluorescence EEM")
        m_eem.addAction(QtGui.QAction("Open EEM matrix file…", self,
                                      triggered=self.open_eem_file))
        m_eem.addAction(QtGui.QAction("Build EEM from loaded spectra…", self,
                                      triggered=self.build_eem_from_spectra))
        m_eem.addAction(QtGui.QAction("Open demo EEM", self, triggered=self.open_demo_eem))
        m_eem.addSeparator()
        m_eem.addAction(QtGui.QAction("PARAFAC on loaded EEM stack…", self,
                                      triggered=self.parafac_analysis))
        m_eem.addAction(QtGui.QAction("Open demo EEM stack (for PARAFAC)", self,
                                      triggered=self.open_demo_eem_stack))
        m_an.addSeparator()
        m_an.addAction(QtGui.QAction("Load demo perturbation series (for 2D-COS)", self,
                                     triggered=self.load_demo_series))
        m_an.addAction(QtGui.QAction("Clear analysis overlays", self,
                                     triggered=self.clear_analysis))

        m_lib = mb.addMenu("&Library")
        m_lib.addAction(QtGui.QAction("Add selected to library", self,
                                      triggered=self.library_add))
        m_lib.addAction(QtGui.QAction("Search selected against library…", self,
                                      triggered=self.library_search))
        m_lib.addSeparator()
        m_lib.addAction(QtGui.QAction("View library…", self, triggered=self.library_view))
        m_lib.addAction(QtGui.QAction("Load library…", self, triggered=self.library_load))
        m_lib.addAction(QtGui.QAction("Save library…", self, triggered=self.library_save))
        m_lib.addAction(QtGui.QAction("Clear library", self, triggered=self.library_clear))

        m_help = mb.addMenu("&Help")
        m_help.addAction(QtGui.QAction("About", self, triggered=self.about))

    def _build_toolbar(self) -> None:
        tb = self.addToolBar("Main")
        tb.setIconSize(QtCore.QSize(18, 18))
        tb.addActions([self.act_open, self.act_save])
        tb.addSeparator()
        tb.addActions([self.act_undo, self.act_redo])
        tb.addSeparator()
        tb.addAction(self.act_autoscale)

    # ===================================================== spectrum table
    def _rebuild_table(self) -> None:
        self.table.blockSignals(True)
        self.table.setRowCount(len(self.document))
        for i, spec in enumerate(self.document):
            chk = QtWidgets.QTableWidgetItem()
            chk.setFlags(QtCore.Qt.ItemFlag.ItemIsUserCheckable
                         | QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable)
            chk.setCheckState(QtCore.Qt.CheckState.Checked if spec.visible
                              else QtCore.Qt.CheckState.Unchecked)
            self.table.setItem(i, 0, chk)

            col = QtWidgets.QTableWidgetItem()
            col.setBackground(QtGui.QColor(spec.color or "#888888"))
            col.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable)
            self.table.setItem(i, 1, col)

            name = QtWidgets.QTableWidgetItem(spec.name)
            name.setToolTip(f"{spec.npoints} pts · {spec.x_label} · {spec.y_label}")
            self.table.setItem(i, 2, name)
        self.table.blockSignals(False)

    def _on_item_changed(self, item: QtWidgets.QTableWidgetItem) -> None:
        i = item.row()
        if not (0 <= i < len(self.document)):
            return
        spec = self.document[i]
        if item.column() == 0:
            spec.visible = item.checkState() == QtCore.Qt.CheckState.Checked
            self.plotview.refresh()
        elif item.column() == 2:
            spec.name = item.text()
            self.plotview.refresh()

    def _on_cell_double_clicked(self, row: int, col: int) -> None:
        if col != 1 or not (0 <= row < len(self.document)):
            return
        spec = self.document[row]
        c = QtWidgets.QColorDialog.getColor(QtGui.QColor(spec.color or "#888888"), self,
                                            "Pick spectrum colour")
        if c.isValid():
            spec.color = c.name()
            self._rebuild_table()
            self.plotview.refresh()

    def _selected_indices(self) -> list[int]:
        rows = sorted({idx.row() for idx in self.table.selectionModel().selectedRows()})
        return [r for r in rows if 0 <= r < len(self.document)]

    def _targets(self) -> list[Spectrum]:
        """Selected spectra, or all spectra if the selection is empty."""
        idx = self._selected_indices()
        return [self.document[i] for i in idx] if idx else list(self.document)

    # ===================================================== file IO
    def open_files(self) -> None:
        paths, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self, "Open spectra", self._last_dir, OPEN_FILTER)
        if paths:
            self._load_paths(paths)

    def import_from_url(self) -> None:
        """Download a spectrum from IRUG (by id/URL) or any direct URL."""
        text, ok = QtWidgets.QInputDialog.getText(
            self, "Import from URL / IRUG",
            "IRUG spectrum id, IRUG page URL, or a direct spectrum URL:\n"
            "(e.g.  3537   or   http://www.irug.org/jcamp-details?id=3537 )")
        if not ok or not text.strip():
            return
        self.status.showMessage(f"Downloading {text.strip()}…", 0)
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
        try:
            specs = load_online(text.strip())
            self.document.add_many(specs)
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QApplication.restoreOverrideCursor()
            self.status.clearMessage()
            QtWidgets.QMessageBox.warning(self, "Import failed", str(exc))
            return
        QtWidgets.QApplication.restoreOverrideCursor()
        self._rebuild_table()
        self.plotview.suggest_flip()
        self.plotview.refresh()
        self.plotview.autoscale()
        self.status.showMessage(f"Imported {len(specs)} spectrum(s) from the web.", 5000)

    def _load_paths(self, paths: list[str]) -> None:
        if paths:
            self._last_dir = os.path.dirname(os.path.abspath(paths[0]))
        loaded, errors = 0, []
        for p in paths:
            try:
                specs = load_any(p)
                self.document.add_many(specs)
                loaded += len(specs)
            except MissingDependency as exc:
                errors.append(str(exc))
            except Exception as exc:  # noqa: BLE001
                errors.append(f"{os.path.basename(p)}: {exc}")
        self._rebuild_table()
        self.plotview.suggest_flip()
        self.plotview.refresh()
        self.plotview.autoscale()
        if loaded:
            self.status.showMessage(f"Loaded {loaded} spectrum(s).", 5000)
        if errors:
            QtWidgets.QMessageBox.warning(self, "Some files could not be read",
                                          "\n\n".join(errors))

    def open_soprano_url(self) -> None:
        default = "https://soprano.kikirpa.be/index.php?lib=sop&id=PR1_E_785_kikirpa"
        url, ok = QtWidgets.QInputDialog.getText(
            self, "Open SOPRANO URL", "SOPRANO spectrum URL:", text=default)
        if not ok or not url.strip():
            return
        try:
            specs = load_soprano_url(url.strip())
            self.document.add_many(specs)
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.critical(self, "SOPRANO read failed", str(exc))
            return
        self._rebuild_table()
        self.plotview.suggest_flip()
        self.plotview.refresh()
        self.plotview.autoscale()
        self.status.showMessage(f"Loaded {len(specs)} SOPRANO spectrum(s).", 5000)

    def load_demo(self) -> None:
        self.document.add_many(load_demo_set())
        self._rebuild_table()
        self.plotview.refresh()
        self.plotview.autoscale()
        self.status.showMessage("Loaded demo spectra (FTIR, Raman, UV/Vis, NIR).", 5000)

    def _default_save_path(self, filename: str) -> str:
        """Default for save dialogs: last-used folder + a sensible filename."""
        return os.path.join(self._last_dir, filename) if self._last_dir else filename

    def _report_saved(self, path: str, detail: str = "") -> None:
        """Show the FULL saved path and offer to open the containing folder."""
        full = os.path.abspath(path)
        self._last_dir = os.path.dirname(full)
        self.status.showMessage(f"Saved to {full}", 0)
        box = QtWidgets.QMessageBox(self)
        box.setIcon(QtWidgets.QMessageBox.Icon.Information)
        box.setWindowTitle("Saved")
        box.setText((detail + "\n\n" if detail else "") + f"File saved to:\n{full}")
        open_btn = box.addButton("Open folder",
                                 QtWidgets.QMessageBox.ButtonRole.ActionRole)
        box.addButton(QtWidgets.QMessageBox.StandardButton.Ok)
        box.setDefaultButton(QtWidgets.QMessageBox.StandardButton.Ok)
        box.exec()
        if box.clickedButton() is open_btn:
            QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(self._last_dir))

    def save_spectrum(self) -> None:
        idx = self._selected_indices()
        if not idx:
            QtWidgets.QMessageBox.information(self, "Save", "Select a spectrum first.")
            return
        spec = self.document[idx[0]]
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save spectrum", self._default_save_path(f"{spec.name}.csv"),
            "CSV (*.csv);;JCAMP-DX (*.dx);;JSON (*.json)")
        if not path:
            return
        try:
            low = path.lower()
            if low.endswith(".dx"):
                save_jcamp(spec, path)
            elif low.endswith(".json"):
                save_json(spec, path)
            else:
                save_csv(spec, path)
            self._report_saved(path)
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.critical(self, "Save failed", str(exc))

    def export_image(self) -> None:
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export plot", self._default_save_path("spectrum.png"),
            "PNG image (*.png);;SVG vector (*.svg);;PDF document (*.pdf);;"
            "EPS vector (*.eps)")
        if not path:
            return
        try:
            self.plotview.export_image(path)
            self._report_saved(path)
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.critical(self, "Export failed", str(exc))

    def export_combined(self) -> None:
        """Export several spectra merged onto one common x-grid as a single CSV."""
        if not len(self.document):
            QtWidgets.QMessageBox.information(self, "Export combined data",
                                              "There are no spectra to export.")
            return
        v = FormDialog.exec_form("Export combined data", [
            {"key": "scope", "label": "Spectra to export", "type": "choice",
             "options": [("All loaded", "all"), ("Visible only", "visible"),
                         ("Selected only", "selected")], "default": "all"},
            {"key": "layout", "label": "Layout", "type": "choice",
             "options": [("x-column + one column per spectrum", "columns"),
                         ("one row per spectrum (X-matrix for PLS/PCA)", "rows")],
             "default": "columns"},
        ], self, "All spectra are put on a shared wavelength grid (interpolated if "
                 "their grids differ).")
        if not v:
            return
        if v["scope"] == "visible":
            specs = self.document.visible_spectra()
        elif v["scope"] == "selected":
            specs = [self.document[i] for i in self._selected_indices()]
        else:
            specs = list(self.document)
        specs = [s for s in specs if s.npoints]
        if not specs:
            QtWidgets.QMessageBox.information(self, "Export combined data",
                                              "That selection contains no spectra.")
            return

        units = {s.y_unit for s in specs}
        if len(units) > 1:
            ans = QtWidgets.QMessageBox.question(
                self, "Mixed Y units",
                "The chosen spectra have different Y units:\n  "
                + ", ".join(sorted(units))
                + "\n\nMerge them into one file anyway?")
            if ans != QtWidgets.QMessageBox.StandardButton.Yes:
                return

        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export combined data",
            self._default_save_path("combined_spectra.csv"), "CSV (*.csv)")
        if not path:
            return
        try:
            info = save_combined_csv(specs, path, layout=v["layout"])
            detail = (f"Exported {info['n_spectra']} spectra × "
                      f"{info['n_points']} points.")
            if info["resampled"]:
                detail += "\n(grids differed → interpolated onto a common grid)"
            self._report_saved(path, detail)
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.critical(self, "Export failed", str(exc))

    # ===================================================== undo / edit
    def _push_undo(self) -> None:
        self._undo.append(self.document.snapshot())
        self._redo.clear()
        if len(self._undo) > 50:
            self._undo.pop(0)

    def undo(self) -> None:
        if not self._undo:
            return
        self._redo.append(self.document.snapshot())
        self.document.restore(self._undo.pop())
        self._rebuild_table()
        self.plotview.refresh()

    def redo(self) -> None:
        if not self._redo:
            return
        self._undo.append(self.document.snapshot())
        self.document.restore(self._redo.pop())
        self._rebuild_table()
        self.plotview.refresh()

    def remove_selected(self) -> None:
        idx = self._selected_indices()
        if not idx:
            return
        self._push_undo()
        for i in reversed(idx):
            self.document.remove_at(i)
        self._rebuild_table()
        self.plotview.refresh()

    def remove_all(self) -> None:
        if not len(self.document):
            return
        self._push_undo()
        self.document.clear()
        self._rebuild_table()
        self.plotview.refresh()

    # ===================================================== apply helpers
    def _apply_inplace(self, func, label: str) -> None:
        """Run ``func(spec) -> Spectrum`` on each target, replacing it in place."""
        targets = self._targets()
        if not targets:
            self.status.showMessage("No spectra to process.", 4000)
            return
        self._push_undo()
        ok, errs = 0, []
        for spec in list(targets):
            try:
                self.document.replace(spec, func(spec))
                ok += 1
            except Exception as exc:  # noqa: BLE001
                errs.append(f"{spec.name}: {exc}")
        self._rebuild_table()
        self.plotview.refresh()
        msg = f"{label}: {ok} spectrum(s)."
        if errs:
            msg += f" {len(errs)} skipped."
            self.status.showMessage(msg, 6000)
            QtWidgets.QMessageBox.warning(self, label, "\n".join(errs))
        else:
            self.status.showMessage(msg, 4000)

    def _add_result(self, spec: Spectrum, label: str) -> None:
        self._push_undo()
        self.document.add(spec)
        self._rebuild_table()
        self.plotview.refresh()
        self.status.showMessage(label, 4000)

    # ===================================================== axis conversions
    def convert_x(self, unit: str) -> None:
        laser = None
        if unit == "raman_cm-1" or any(s.x_unit == "raman_cm-1" for s in self._targets()):
            laser = self._ask_laser()
            if laser is None:
                return
        self._apply_inplace(lambda s: axes.convert_x(s, unit, laser),
                            f"Convert X → {X_UNIT_LABELS.get(unit, unit)}")

    def convert_y(self, unit: str) -> None:
        self._apply_inplace(lambda s: axes.convert_y(s, unit),
                            f"Convert Y → {Y_UNIT_LABELS.get(unit, unit)}")

    def _ask_laser(self) -> float | None:
        val, ok = QtWidgets.QInputDialog.getDouble(
            self, "Raman excitation", "Laser wavelength (nm):", 785.0, 100.0, 12000.0, 1)
        return val if ok else None

    # ===================================================== processing ops
    def smooth_sg(self) -> None:
        v = FormDialog.exec_form("Savitzky-Golay smoothing", [
            {"key": "window", "label": "Window length (odd)", "type": "int",
             "default": 11, "min": 3, "max": 999, "step": 2},
            {"key": "polyorder", "label": "Polynomial order", "type": "int",
             "default": 3, "min": 1, "max": 9},
        ], self)
        if v:
            self._apply_inplace(lambda s: processing.savitzky_golay(
                s, v["window"], v["polyorder"]), "Savitzky-Golay")

    def smooth_ma(self) -> None:
        v = FormDialog.exec_form("Moving average", [
            {"key": "window", "label": "Window length", "type": "int",
             "default": 5, "min": 1, "max": 999}], self)
        if v:
            self._apply_inplace(lambda s: processing.moving_average(s, v["window"]),
                                "Moving average")

    def baseline_rb(self) -> None:
        self._apply_inplace(processing.baseline_rubberband, "Baseline rubberband")

    def baseline_poly(self) -> None:
        v = FormDialog.exec_form("Polynomial baseline", [
            {"key": "order", "label": "Polynomial order", "type": "int",
             "default": 3, "min": 1, "max": 12}], self,
            "Iteratively fits a polynomial below the spectrum (ModPoly).")
        if v:
            self._apply_inplace(lambda s: processing.baseline_polynomial(s, v["order"]),
                                "Baseline polynomial")

    def baseline_als(self) -> None:
        v = FormDialog.exec_form("Asymmetric Least Squares baseline", [
            {"key": "lam", "label": "Smoothness λ", "type": "float",
             "default": 1e5, "min": 1.0, "max": 1e12, "decimals": 0, "step": 1e4},
            {"key": "p", "label": "Asymmetry p", "type": "float",
             "default": 0.01, "min": 0.0001, "max": 0.5, "decimals": 4, "step": 0.005},
        ], self, "Eilers & Boelens ALS. Larger λ = smoother; smaller p = follows valleys.")
        if v:
            self._apply_inplace(lambda s: processing.baseline_als(s, v["lam"], v["p"]),
                                "Baseline ALS")

    def baseline_airpls(self) -> None:
        v = FormDialog.exec_form("airPLS baseline", [
            {"key": "lam", "label": "Smoothness λ", "type": "float",
             "default": 1e5, "min": 1.0, "max": 1e9, "decimals": 0, "step": 1e4},
            {"key": "porder", "label": "Difference order (2 = robust)", "type": "int",
             "default": 2, "min": 1, "max": 3},
            {"key": "n_iter", "label": "Max iterations", "type": "int",
             "default": 15, "min": 1, "max": 100},
        ], self, "Adaptive iteratively reweighted penalised least squares "
                 "(Zhang et al. 2010). Great for fluorescence-style backgrounds.")
        if v:
            self._apply_inplace(
                lambda s: processing.baseline_airpls(s, v["lam"], int(v["porder"]),
                                                     int(v["n_iter"])),
                "Baseline airPLS")

    def derivative(self) -> None:
        v = FormDialog.exec_form("Derivative (Savitzky-Golay)", [
            {"key": "order", "label": "Derivative order", "type": "int",
             "default": 1, "min": 1, "max": 4},
            {"key": "window", "label": "Window length (odd)", "type": "int",
             "default": 11, "min": 5, "max": 999, "step": 2},
            {"key": "polyorder", "label": "Polynomial order", "type": "int",
             "default": 3, "min": 2, "max": 9},
        ], self)
        if v:
            self._apply_inplace(lambda s: processing.derivative(
                s, v["order"], v["window"], v["polyorder"]), "Derivative")

    def normalize(self, method: str) -> None:
        self._apply_inplace(lambda s: processing.normalize(s, method),
                            f"Normalize ({method})")

    def normalize_value(self) -> None:
        v = FormDialog.exec_form("Normalize at value", [
            {"key": "x0", "label": "x position", "type": "float", "default": 0.0,
             "decimals": 3}], self, "Scales each spectrum so y(x₀) = 1.")
        if v:
            self._apply_inplace(lambda s: processing.normalize(s, "value", v["x0"]),
                                "Normalize at value")

    def snv(self) -> None:
        self._apply_inplace(processing.snv, "SNV")

    def msc(self) -> None:
        targets = self._targets()
        if len(targets) < 2:
            QtWidgets.QMessageBox.information(
                self, "MSC", "MSC needs at least two spectra (uses their mean as reference).")
            return
        ref = processing.average(targets, "mean")
        self._apply_inplace(lambda s: processing.msc(s, ref), "MSC")

    def detrend(self) -> None:
        v = FormDialog.exec_form("Detrend", [
            {"key": "order", "label": "Trend polynomial order", "type": "int",
             "default": 1, "min": 1, "max": 6}], self)
        if v:
            self._apply_inplace(lambda s: processing.detrend(s, v["order"]), "Detrend")

    def arithmetic(self, op: str) -> None:
        idx = self._selected_indices()
        if len(idx) != 2:
            QtWidgets.QMessageBox.information(
                self, "Arithmetic", "Select exactly two spectra (A then B).")
            return
        a, b = self.document[idx[0]], self.document[idx[1]]
        try:
            res = processing.combine(a, b, op)
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.critical(self, "Arithmetic failed", str(exc))
            return
        self._add_result(res, f"Created {res.name}")

    def average(self, method: str) -> None:
        targets = self._targets()
        if len(targets) < 2:
            QtWidgets.QMessageBox.information(self, "Average", "Select at least two spectra.")
            return
        try:
            res = processing.average(targets, method)
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.critical(self, "Average failed", str(exc))
            return
        self._add_result(res, f"Created {res.name}")

    def cut_range(self) -> None:
        v = FormDialog.exec_form("Cut range", [
            {"key": "lo", "label": "From x", "type": "float", "default": 0.0, "decimals": 3},
            {"key": "hi", "label": "To x", "type": "float", "default": 1000.0, "decimals": 3},
            {"key": "keep", "label": "Keep inside (else remove)", "type": "bool",
             "default": True},
        ], self)
        if v:
            self._apply_inplace(lambda s: processing.cut(s, v["lo"], v["hi"], v["keep"]),
                                "Cut range")

    def interpolate(self) -> None:
        v = FormDialog.exec_form("Interpolate / resample", [
            {"key": "n", "label": "Number of points", "type": "int",
             "default": 1024, "min": 8, "max": 1_000_000}], self,
            "Resample onto a uniform grid with this many points.")
        if v:
            self._apply_inplace(lambda s: processing.interpolate(s, n=v["n"]),
                                "Interpolate")

    # ===================================================== view helpers
    def set_stack(self) -> None:
        v = FormDialog.exec_form("Stack / offset", [
            {"key": "offset", "label": "Offset per spectrum (× y-span)", "type": "float",
             "default": self.plotview.stack_offset, "min": 0.0, "max": 5.0,
             "decimals": 2, "step": 0.1}], self,
            "0 = overlay. Try 0.2–1.0 to separate spectra vertically.")
        if v:
            self.plotview.set_stack_offset(v["offset"])

    def set_dark(self, on: bool) -> None:
        self._dark = on
        self.plotview.plot_widget.setBackground("#202225" if on else "w")

    def _on_cursor(self, x: float, y: float) -> None:
        vis = self.document.visible_spectra()
        xlab = vis[0].x_label if vis else "x"
        ylab = vis[0].y_label if vis else "y"
        self._cursor_label.setText(f"{xlab}: {x:.4g}    {ylab}: {y:.4g}")

    def about(self) -> None:
        QtWidgets.QMessageBox.about(
            self, f"About {__app_name__}",
            f"<b>{__app_name__} {__version__}</b><br>"
            "A SpectraGryph-inspired spectroscopy viewer.<br><br>"
            "Built with PySide6, pyqtgraph, NumPy and SciPy.<br>"
            "Load FTIR / Raman / UV-Vis / NIR spectra, convert axes, "
            "process and analyse.")

    # ===================================================== analysis
    def _analysis_target(self) -> Spectrum | None:
        """Spectrum to analyse: the first selected one, else the first visible."""
        idx = self._selected_indices()
        if idx:
            return self.document[idx[0]]
        vis = self.document.visible_spectra()
        return vis[0] if vis else None

    def _show_dialog(self, dlg) -> None:
        self._dialogs = [d for d in self._dialogs if d.isVisible()]
        self._dialogs.append(dlg)
        dlg.show()

    def find_peaks(self) -> None:
        spec = self._analysis_target()
        if spec is None:
            QtWidgets.QMessageBox.information(self, "Find peaks", "Load a spectrum first.")
            return
        v = FormDialog.exec_form("Find peaks", [
            {"key": "height", "label": "Min height (% of range)", "type": "float",
             "default": 5.0, "min": 0.0, "max": 100.0, "decimals": 1, "step": 1.0},
            {"key": "prom", "label": "Min prominence (% of range)", "type": "float",
             "default": 3.0, "min": 0.0, "max": 100.0, "decimals": 1, "step": 1.0},
            {"key": "dist", "label": "Min distance (x units, 0 = off)", "type": "float",
             "default": 0.0, "min": 0.0, "decimals": 3},
            {"key": "smooth", "label": "Pre-smooth window (0 = off)", "type": "int",
             "default": 0, "min": 0, "max": 199, "step": 2},
            {"key": "valleys", "label": "Find valleys (dips) instead", "type": "bool",
             "default": False},
        ], self, f"Analysing: {spec.name}")
        if not v:
            return
        peaks = analysis.find_peaks(spec, v["height"] / 100.0, v["prom"] / 100.0,
                                    v["dist"], v["valleys"], int(v["smooth"]))
        self.plotview.mark_peaks(peaks)
        if not peaks:
            self.status.showMessage("No peaks found — try lowering the thresholds.", 6000)
            return
        rows = [[i + 1, f"{p.center:.6g}", f"{p.height:.6g}", f"{p.fwhm:.6g}",
                 f"{p.area:.6g}"] for i, p in enumerate(peaks)]
        self._show_dialog(TableDialog(
            f"Peaks — {spec.name}",
            ["#", f"Position ({spec.x_unit})", f"Height ({spec.y_unit})", "FWHM",
             "Area (approx)"], rows, self,
            summary=f"{len(peaks)} peak(s) found in “{spec.name}”.",
            default_dir=self._last_dir))
        self.status.showMessage(f"Found {len(peaks)} peaks in {spec.name}.", 5000)

    def fit_peaks(self) -> None:
        spec = self._analysis_target()
        if spec is None:
            QtWidgets.QMessageBox.information(self, "Peak fit", "Load a spectrum first.")
            return
        v = FormDialog.exec_form("Peak fit / deconvolution", [
            {"key": "model", "label": "Peak shape", "type": "choice",
             "options": [("Gaussian", "gaussian"), ("Lorentzian", "lorentzian"),
                         ("Pseudo-Voigt", "pseudovoigt")], "default": "gaussian"},
            {"key": "max_peaks", "label": "Max peaks", "type": "int",
             "default": 6, "min": 1, "max": 30},
            {"key": "baseline", "label": "Also fit a linear baseline", "type": "bool",
             "default": False},
        ], self, f"Analysing: {spec.name}. Initial peaks come from auto detection.")
        if not v:
            return
        try:
            fit = analysis.fit_peaks(spec, v["model"], int(v["max_peaks"]), v["baseline"])
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.warning(self, "Peak fit failed", str(exc))
            return
        self.plotview.show_fit(fit.x, fit.total, fit.comp_curves, fit.baseline)
        rows = [[i + 1, f"{c.center:.6g}", f"{c.amplitude:.6g}", f"{c.fwhm:.6g}",
                 f"{c.area:.6g}", "" if c.eta is None else f"{c.eta:.3f}"]
                for i, c in enumerate(fit.components)]

        def add_components():
            self._push_undo()
            for i, curve in enumerate(fit.comp_curves):
                self.document.add(analysis.component_to_spectrum(
                    spec, fit.x, curve, f"{spec.name} · {v['model']} {i + 1}"))
            self._rebuild_table()
            self.plotview.refresh()
            self.status.showMessage(
                f"Added {len(fit.comp_curves)} fit components to the list.", 5000)

        self._show_dialog(TableDialog(
            f"Peak fit ({v['model']}) — {spec.name}",
            ["#", f"Center ({spec.x_unit})", "Height", "FWHM", "Area", "η"], rows, self,
            summary=f"{len(fit.components)} components · R² = {fit.r_squared:.5f}",
            extra_buttons=[("Add components to list", add_components)],
            default_dir=self._last_dir))
        self.status.showMessage(
            f"Fitted {len(fit.components)} {v['model']} peaks, R²={fit.r_squared:.4f}.",
            6000)

    def integrate_range(self) -> None:
        spec = self._analysis_target()
        if spec is None:
            QtWidgets.QMessageBox.information(self, "Integrate", "Load a spectrum first.")
            return
        x0, x1 = spec.xrange
        v = FormDialog.exec_form("Integrate range", [
            {"key": "lo", "label": "From x", "type": "float", "default": x0,
             "decimals": 4},
            {"key": "hi", "label": "To x", "type": "float", "default": x1,
             "decimals": 4},
            {"key": "linear", "label": "Subtract linear baseline (endpoints)",
             "type": "bool", "default": True},
        ], self, f"Analysing: {spec.name}")
        if not v:
            return
        try:
            res = analysis.integrate(spec, v["lo"], v["hi"],
                                     "linear" if v["linear"] else "none")
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.warning(self, "Integration failed", str(exc))
            return
        self.plotview.show_region(res["x_lo"], res["x_hi"])
        QtWidgets.QMessageBox.information(
            self, "Integration result",
            f"Spectrum: {spec.name}\n"
            f"Range: {res['x_lo']:.6g} – {res['x_hi']:.6g} {spec.x_unit}\n"
            f"Baseline: {res['baseline']}\n\n"
            f"Area = {res['area']:.6g}\n"
            f"Centroid = {res['centroid']:.6g} {spec.x_unit}")
        self.status.showMessage(
            f"Area [{res['x_lo']:.4g}, {res['x_hi']:.4g}] = {res['area']:.6g}", 0)

    def clear_analysis(self) -> None:
        self.plotview.clear_analysis()
        self.status.showMessage("Cleared analysis overlays.", 3000)

    # ---- 2D correlation spectroscopy ------------------------------------
    def cos2d_analysis(self) -> None:
        idx = self._selected_indices()
        targets = [self.document[i] for i in idx] if idx else list(self.document)
        targets = [s for s in targets if s.npoints]
        if len(targets) < 2:
            QtWidgets.QMessageBox.information(
                self, "2D correlation",
                "Select at least two spectra (a perturbation series).\n"
                "Tip: Analyze ▸ Load demo perturbation series.")
            return
        methods = [("Generalized 2D-COS (Noda)", "generalized")]
        if len(targets) == 2:
            methods.append(("Two-trace 2D (2T2D)", "2t2d"))
        if len(targets) >= 4 and len(targets) % 2 == 0:
            methods.append(("Hetero-correlation (1st half × 2nd half)", "hetero"))
        v = FormDialog.exec_form("2D correlation analysis", [
            {"key": "method", "label": "Method", "type": "choice",
             "options": methods, "default": "generalized"},
            {"key": "ref", "label": "Reference (generalized / hetero)", "type": "choice",
             "options": [("Mean → dynamic spectra", "mean"), ("None", "none")],
             "default": "mean"},
        ], self, f"{len(targets)} spectra selected. (Hetero splits them in half: the "
                 "first set is one technique, the second the other — pair by sample order.)")
        if not v:
            return
        try:
            if v["method"] == "2t2d":
                x, sync, asyn = cos2d.two_trace_from_spectra(targets[0], targets[1])
                xl = yl = targets[0].x_label
                tag = "2T2D"
            elif v["method"] == "hetero":
                h = len(targets) // 2
                g1, g2 = targets[:h], targets[h:]
                x, y, sync, asyn = cos2d.hetero_from_spectra(g1, g2, ref=v["ref"])
                xl, yl = g1[0].x_label, g2[0].x_label
                tag = "Hetero 2D-COS"
            else:
                x, sync, asyn = cos2d.correlation_from_spectra(targets, ref=v["ref"])
                xl = yl = targets[0].x_label
                tag = "2D-COS"
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.warning(self, "2D correlation failed", str(exc))
            return
        y = x if v["method"] != "hetero" else y
        panels = [
            {"x": x, "y": y, "Z": sync, "label": "Synchronous  Φ",
             "xlabel": xl, "ylabel": yl, "diverging": True},
            {"x": x, "y": y, "Z": asyn, "label": "Asynchronous  Ψ",
             "xlabel": xl, "ylabel": yl, "diverging": True},
        ]
        self._show_dialog(MapWindow(panels, title=f"{tag} — {len(targets)} spectra",
                                    parent=self))
        self.status.showMessage(f"{tag}: synchronous + asynchronous map.", 6000)

    # ---- fluorescence EEM -----------------------------------------------
    def _show_eem(self, e) -> None:
        self.eems.append(e)
        self._show_dialog(EEMWindow(e, parent=self))

    def open_demo_eem_stack(self) -> None:
        stack = demo_eem_stack(7)
        self.eems = list(stack)        # a self-contained demo stack (shared grid)
        self._show_dialog(EEMWindow(stack[0], parent=self))
        self.status.showMessage(
            f"Loaded a {len(stack)}-EEM demo stack (shared grid). "
            "Run EEM ▸ PARAFAC on loaded EEM stack.", 9000)

    def parafac_analysis(self) -> None:
        if len(self.eems) < 2:
            QtWidgets.QMessageBox.information(
                self, "PARAFAC",
                "PARAFAC needs at least two EEMs sharing the same grid. Open several "
                "EEM files, or use EEM ▸ Open demo EEM stack.")
            return
        maxr = max(1, min(8, len(self.eems)))
        v = FormDialog.exec_form("PARAFAC decomposition", [
            {"key": "rank", "label": "Number of components", "type": "int",
             "default": min(3, maxr), "min": 1, "max": maxr},
            {"key": "nonneg", "label": "Non-negative loadings (fluorescence)",
             "type": "bool", "default": True},
        ], self, f"{len(self.eems)} EEMs loaded.")
        if not v:
            return
        try:
            res = eem.parafac_from_eems(self.eems, int(v["rank"]), nonneg=v["nonneg"])
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.warning(self, "PARAFAC failed", str(exc))
            return
        self._show_dialog(ParafacWindow(res, parent=self))
        self.status.showMessage(f"PARAFAC: {res.rank} components, fit {res.fit:.4f} "
                                f"from {len(self.eems)} EEMs.", 8000)

    def open_eem_file(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open EEM matrix", self._last_dir,
            "EEM matrix (*.csv *.txt *.dat *.tsv);;All files (*.*)")
        if not path:
            return
        v = FormDialog.exec_form("EEM layout", [
            {"key": "exc", "label": "Excitation axis is in…", "type": "choice",
             "options": [("Columns (first row = excitation)", "cols"),
                         ("Rows (first column = excitation)", "rows")], "default": "cols"},
        ], self, "How the matrix file is laid out.")
        if not v:
            return
        try:
            e = eem.read_eem_matrix(path, ex_in_columns=(v["exc"] == "cols"))
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.critical(self, "EEM read failed", str(exc))
            return
        self._last_dir = os.path.dirname(os.path.abspath(path))
        self._show_eem(e)

    def build_eem_from_spectra(self) -> None:
        targets = [self.document[i] for i in self._selected_indices()] or \
            list(self.document)
        targets = [s for s in targets if s.npoints]
        if len(targets) < 2:
            QtWidgets.QMessageBox.information(
                self, "Build EEM",
                "Select at least two emission spectra (each at one excitation "
                "wavelength; excitation read from the name like 'ex280' or meta).")
            return
        try:
            e = eem.eem_from_spectra(targets)
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.warning(self, "Build EEM failed", str(exc))
            return
        self._show_eem(e)

    def open_demo_eem(self) -> None:
        self._show_eem(demo_eem())

    def load_demo_series(self) -> None:
        self.document.add_many(demo_cos_series(12))
        self._rebuild_table()
        self.plotview.refresh()
        self.plotview.autoscale()
        self.status.showMessage("Loaded a 12-spectrum demo perturbation series — "
                                "select them all, then Analyze ▸ 2D correlation.", 9000)

    def mixture_analysis(self) -> None:
        idx = self._selected_indices()
        if len(idx) >= 2:
            mixture, refs, src = self.document[idx[0]], \
                [self.document[i] for i in idx[1:]], "selected"
        elif len(idx) == 1 and len(self.library):
            mixture, refs, src = self.document[idx[0]], list(self.library.entries), \
                "library"
        else:
            QtWidgets.QMessageBox.information(
                self, "Mixture analysis",
                "Select the mixture plus its reference components (the FIRST "
                "selected is the mixture), or select just the mixture and load a "
                "reference library.")
            return
        try:
            res = analysis.mixture_nnls(mixture, refs, fit_offset=True)
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.warning(self, "Mixture analysis failed", str(exc))
            return
        self.plotview.show_fit(res.x, res.fit, [])
        order = np.argsort(res.fractions)[::-1]
        rows = [[i + 1, res.names[j], f"{res.coeffs[j]:.4g}",
                 f"{res.fractions[j] * 100:.2f}%"] for i, j in enumerate(order)]
        summary = (f"Mixture “{mixture.name}” · {len(refs)} references ({src}) · "
                   f"R²={res.r_squared:.4f} · offset={res.offset:.4g}")
        self._show_dialog(TableDialog(
            f"Mixture composition — {mixture.name}",
            ["#", "Component", "Coefficient", "Fraction"], rows, self,
            summary=summary, default_dir=self._last_dir))
        self.status.showMessage(summary, 9000)

    def identify_xrf(self) -> None:
        spec = self._analysis_target()
        if spec is None:
            QtWidgets.QMessageBox.information(self, "Identify XRF", "Load a spectrum first.")
            return
        if spec.x_unit not in ("keV", "eV"):
            QtWidgets.QMessageBox.warning(
                self, "Identify XRF",
                f"XRF element ID needs an energy axis in keV or eV; this spectrum is "
                f"in '{spec.x_unit}'. Convert/calibrate the x-axis to energy first.")
            return
        v = FormDialog.exec_form("Identify XRF elements", [
            {"key": "height", "label": "Min peak height (% of range)", "type": "float",
             "default": 3.0, "min": 0.0, "max": 100.0, "decimals": 1},
            {"key": "prom", "label": "Min prominence (% of range)", "type": "float",
             "default": 2.0, "min": 0.0, "max": 100.0, "decimals": 1},
            {"key": "tol", "label": "Match tolerance (keV)", "type": "float",
             "default": 0.10, "min": 0.01, "max": 1.0, "decimals": 3},
        ], self, f"Analysing: {spec.name}")
        if not v:
            return
        peaks = analysis.find_peaks(spec, v["height"] / 100.0, v["prom"] / 100.0)
        if not peaks:
            self.status.showMessage("No XRF peaks found — lower the thresholds.", 6000)
            return
        scale = 1.0 if spec.x_unit == "keV" else 1e-3   # eV -> keV
        ident = xrf.identify_peaks([p.center * scale for p in peaks], tol=v["tol"])
        labels = [f"{r['best']['symbol']} {r['best']['line_label']}" if r["best"]
                  else "?" for r in ident]
        self.plotview.mark_peaks(peaks, color="#1f77b4", labels=labels)
        rows = []
        for i, (p, r) in enumerate(zip(peaks, ident)):
            b = r["best"]
            rows.append([i + 1, f"{p.center * scale:.4g}",
                         f"{b['symbol']} ({b['name']})" if b else "—",
                         b["line_label"] if b else "—",
                         f"{b['energy']:.4g}" if b else "—",
                         f"{b['delta'] * 1000:+.0f}" if b else "—"])
        n_match = sum(1 for r in ident if r["best"])
        self._show_dialog(TableDialog(
            f"XRF elements — {spec.name}",
            ["#", "Peak (keV)", "Element", "Line", "Line (keV)", "ΔE (eV)"], rows, self,
            summary=f"{len(peaks)} peaks, {n_match} matched to elements.",
            default_dir=self._last_dir))
        self.status.showMessage(f"XRF: matched {n_match}/{len(peaks)} peaks.", 6000)

    # ---- calibration curve (檢量線) -------------------------------------
    def calibration_curve(self) -> None:
        specs = [s for s in self._targets() if s.npoints]
        if len(specs) < 2:
            QtWidgets.QMessageBox.information(
                self, "Calibration curve",
                "Load the standard spectra (and any unknown samples) and select "
                "them first. Each standard needs a known concentration — put it in "
                "the spectrum name (e.g. 'std 5 ppm') or type it into the table.")
            return
        res = CalibrationDialog.exec_dialog(specs, self, self._last_dir)
        if not res:
            return
        stds = res["standards"]
        need = 2 if res["degree"] == 1 else 3
        if len(stds) < need:
            QtWidgets.QMessageBox.warning(
                self, "Calibration curve",
                f"Need at least {need} standards with a concentration and a signal "
                f"for a {'quadratic' if res['degree'] == 2 else 'linear'} fit; "
                f"got {len(stds)}.")
            return
        try:
            model = calibration.fit_calibration(
                [c for c, _, _ in stds], [sig for _, sig, _ in stds],
                degree=res["degree"], through_origin=res["through_origin"],
                conc_unit=res["conc_unit"], signal_label=res["signal_label"])
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.warning(self, "Calibration failed", str(exc))
            return
        preds = []
        for name, sig in res["unknowns"]:
            preds.append({"name": name, "signal": sig,
                          **model.predict_concentration(sig, res["replicates"])})
        self._show_dialog(CalibrationWindow(model, preds, parent=self,
                                            last_dir=self._last_dir))
        unit = f" {model.conc_unit}" if model.conc_unit else ""
        msg = f"Calibration: {model.equation()}  R²={model.r_squared:.4f}"
        if model.lod is not None:
            msg += f"  LOD={model.lod:.3g}{unit}"
        self.status.showMessage(msg, 0)

    # ---- spectral library -----------------------------------------------
    def library_add(self) -> None:
        targets = [self.document[i] for i in self._selected_indices()] or \
            list(self.document)
        targets = [s for s in targets if s.npoints]
        if not targets:
            QtWidgets.QMessageBox.information(self, "Library", "No spectra to add.")
            return
        for s in targets:
            self.library.add(s)
        self.status.showMessage(
            f"Added {len(targets)} → library now has {len(self.library)} entries.", 6000)

    def library_search(self) -> None:
        spec = self._analysis_target()
        if spec is None:
            QtWidgets.QMessageBox.information(self, "Library search", "Load a spectrum first.")
            return
        if not len(self.library):
            QtWidgets.QMessageBox.information(
                self, "Library search",
                "The library is empty. Add spectra (Library ▸ Add selected) or load a "
                "library file first.")
            return
        hits = self.library.search(spec, top_n=15)
        rows = [[i + 1, h["name"], f"{h['scores']['correlation']:.4f}",
                 f"{h['scores']['cosine']:.4f}", f"{h['scores']['sam']:.4f}",
                 f"{h['scores']['euclid']:.4f}"] for i, h in enumerate(hits)]

        def overlay_top():
            if hits:
                top = hits[0]["entry"].copy()
                top.name, top.color = f"[ref] {hits[0]['name']}", None
                self.document.add(top)
                self._rebuild_table()
                self.plotview.refresh()
                self.status.showMessage(f"Overlaid: {hits[0]['name']}", 5000)

        top_txt = (f"Best match: {hits[0]['name']} "
                   f"(correlation {hits[0]['scores']['correlation']:.4f})") if hits else ""
        self._show_dialog(TableDialog(
            f"Library search — {spec.name}",
            ["#", "Reference", "Correlation", "Cosine", "SAM (rad)", "Euclid"], rows, self,
            summary=top_txt, extra_buttons=[("Overlay top hit", overlay_top)],
            default_dir=self._last_dir))

    def library_view(self) -> None:
        if not len(self.library):
            QtWidgets.QMessageBox.information(self, "Library", "The library is empty.")
            return
        rows = [[i + 1, s.name, s.npoints, s.x_unit, s.y_unit]
                for i, s in enumerate(self.library.entries)]
        self._show_dialog(TableDialog(
            f"Library: {self.library.name} ({len(self.library)} entries)",
            ["#", "Name", "Points", "X unit", "Y unit"], rows, self,
            default_dir=self._last_dir))

    def library_load(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Load library", self._last_dir,
            "Spectral library (*.speclib *.json);;All files (*.*)")
        if not path:
            return
        try:
            self.library = SpectralLibrary.load(path)
            self._last_dir = os.path.dirname(os.path.abspath(path))
            self.status.showMessage(
                f"Loaded library '{self.library.name}' ({len(self.library)} entries).", 7000)
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.critical(self, "Load failed", str(exc))

    def library_save(self) -> None:
        if not len(self.library):
            QtWidgets.QMessageBox.information(self, "Library", "The library is empty.")
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save library", self._default_save_path("references.speclib"),
            "Spectral library (*.speclib);;JSON (*.json)")
        if not path:
            return
        try:
            self.library.save(path)
            self._report_saved(path)
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.critical(self, "Save failed", str(exc))

    def library_clear(self) -> None:
        self.library.clear()
        self.status.showMessage("Library cleared.", 3000)

    # ===================================================== drag & drop
    def dragEnterEvent(self, event: QtGui.QDragEnterEvent) -> None:
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event: QtGui.QDropEvent) -> None:
        paths = [u.toLocalFile() for u in event.mimeData().urls() if u.isLocalFile()]
        if paths:
            self._load_paths(paths)
