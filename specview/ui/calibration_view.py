"""Calibration-curve UI: a setup dialog and the result window with the plot.

``CalibrationDialog`` collects, for a set of selected spectra, how to read each
one's signal and which are standards (known concentration) vs. unknowns. The
main window then fits a :class:`~specview.calibration.CalibrationModel` and shows
``CalibrationWindow`` — a scatter of the standards, the fitted line, and the
predicted unknowns with confidence-interval error bars.
"""
from __future__ import annotations

import csv
import os

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

from .. import calibration
from .dialogs import TableDialog

_ROLES = ["Standard", "Unknown", "Ignore"]


class CalibrationDialog(QtWidgets.QDialog):
    """Choose the signal source and tag each spectrum as standard / unknown."""

    def __init__(self, spectra, parent=None, last_dir: str = ""):
        super().__init__(parent)
        self.setWindowTitle("Calibration curve (檢量線)")
        self.resize(620, 460)
        self.spectra = list(spectra)
        self._signals = [float("nan")] * len(self.spectra)

        # sensible default x positions from the loaded data
        xs = np.concatenate([s.x for s in self.spectra if s.npoints]) \
            if any(s.npoints for s in self.spectra) else np.array([0.0, 1.0])
        xmin, xmax = float(xs.min()), float(xs.max())
        first = next((s for s in self.spectra if s.npoints), None)
        x0_default = float(first.x[int(np.argmax(first.y))]) if first is not None else xmin

        layout = QtWidgets.QVBoxLayout(self)
        intro = QtWidgets.QLabel(
            "Read one signal from each spectrum, then mark the standards (known "
            "concentration) and the unknowns to predict. Concentrations are guessed "
            "from the spectrum name — edit any cell to correct them.")
        intro.setWordWrap(True)
        intro.setStyleSheet("color: palette(mid);")
        layout.addWidget(intro)

        # ---- signal-source controls ----
        form = QtWidgets.QFormLayout()
        self.mode = QtWidgets.QComboBox()
        for label, val in calibration.SIGNAL_MODES:
            self.mode.addItem(label, val)
        self.mode.currentIndexChanged.connect(self._on_mode_changed)
        form.addRow("Signal", self.mode)

        self.x0 = self._spin(x0_default, xmin, xmax)
        form.addRow("x₀ (value at)", self.x0)
        self.lo = self._spin(xmin, xmin, xmax)
        self.hi = self._spin(xmax, xmin, xmax)
        rng = QtWidgets.QHBoxLayout()
        rng.addWidget(self.lo)
        rng.addWidget(QtWidgets.QLabel("to"))
        rng.addWidget(self.hi)
        rng_w = QtWidgets.QWidget()
        rng_w.setLayout(rng)
        form.addRow("Range (height / area)", rng_w)
        self.baseline = QtWidgets.QCheckBox("Subtract endpoint baseline (height / area)")
        form.addRow("", self.baseline)

        self.degree = QtWidgets.QComboBox()
        self.degree.addItem("Linear", 1)
        self.degree.addItem("Quadratic", 2)
        form.addRow("Model", self.degree)
        self.origin = QtWidgets.QCheckBox("Force through origin (linear)")
        form.addRow("", self.origin)
        self.unit = QtWidgets.QLineEdit()
        self.unit.setPlaceholderText("e.g. ppm, mg/L")
        form.addRow("Concentration unit", self.unit)
        self.replicates = QtWidgets.QSpinBox()
        self.replicates.setRange(1, 999)
        form.addRow("Unknown replicates (for CI)", self.replicates)
        layout.addLayout(form)

        read_btn = QtWidgets.QPushButton("Read signals")
        read_btn.clicked.connect(self._read_signals)
        layout.addWidget(read_btn, alignment=QtCore.Qt.AlignmentFlag.AlignLeft)

        # ---- per-spectrum table ----
        self.table = QtWidgets.QTableWidget(len(self.spectra), 4)
        self.table.setHorizontalHeaderLabels(["Role", "Spectrum", "Concentration", "Signal"])
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)
        self._roles: list[QtWidgets.QComboBox] = []
        for row, spec in enumerate(self.spectra):
            conc = calibration.parse_concentration(spec.name)
            combo = QtWidgets.QComboBox()
            combo.addItems(_ROLES)
            combo.setCurrentText("Standard" if conc is not None else "Unknown")
            self.table.setCellWidget(row, 0, combo)
            self._roles.append(combo)

            name = QtWidgets.QTableWidgetItem(spec.name)
            name.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable)
            self.table.setItem(row, 1, name)

            citem = QtWidgets.QTableWidgetItem("" if conc is None else f"{conc:g}")
            self.table.setItem(row, 2, citem)

            sitem = QtWidgets.QTableWidgetItem("—")
            sitem.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable)
            self.table.setItem(row, 3, sitem)
        self.table.resizeColumnsToContents()
        layout.addWidget(self.table, 1)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self._on_mode_changed()
        self._read_signals()

    @staticmethod
    def _spin(value: float, lo: float, hi: float) -> QtWidgets.QDoubleSpinBox:
        w = QtWidgets.QDoubleSpinBox()
        w.setDecimals(4)
        span = (hi - lo) or 1.0
        w.setRange(lo - span, hi + span)
        w.setValue(value)
        return w

    def _on_mode_changed(self, *_) -> None:
        is_value = self.mode.currentData() == "value"
        self.x0.setEnabled(is_value)
        self.lo.setEnabled(not is_value)
        self.hi.setEnabled(not is_value)
        self.baseline.setEnabled(not is_value)

    def _read_signals(self) -> None:
        mode = self.mode.currentData()
        x0, lo, hi = self.x0.value(), self.lo.value(), self.hi.value()
        sb = self.baseline.isChecked()
        for row, spec in enumerate(self.spectra):
            try:
                sig = calibration.read_signal(spec, mode, x0=x0, lo=lo, hi=hi,
                                              subtract_baseline=sb)
            except Exception:  # noqa: BLE001
                sig = float("nan")
            self._signals[row] = sig
            self.table.item(row, 3).setText("—" if not np.isfinite(sig) else f"{sig:.6g}")

    def _signal_label(self) -> str:
        mode = self.mode.currentData()
        if mode == "value":
            return f"signal @ {self.x0.value():g}"
        kind = "area" if mode == "area" else "height"
        return f"{kind} {self.lo.value():g}–{self.hi.value():g}"

    def result(self) -> dict:
        """Collect the user's choices into the dict consumed by the main window."""
        self._read_signals()
        standards, unknowns = [], []
        for row, spec in enumerate(self.spectra):
            role = self._roles[row].currentText()
            sig = self._signals[row]
            if role == "Standard":
                try:
                    conc = float(self.table.item(row, 2).text())
                except ValueError:
                    conc = float("nan")
                if np.isfinite(conc) and np.isfinite(sig):
                    standards.append((conc, sig, spec.name))
            elif role == "Unknown" and np.isfinite(sig):
                unknowns.append((spec.name, sig))
        return {
            "standards": standards,
            "unknowns": unknowns,
            "degree": int(self.degree.currentData()),
            "through_origin": self.origin.isChecked(),
            "replicates": int(self.replicates.value()),
            "conc_unit": self.unit.text().strip(),
            "signal_label": self._signal_label(),
        }

    @staticmethod
    def exec_dialog(spectra, parent=None, last_dir: str = "") -> dict | None:
        dlg = CalibrationDialog(spectra, parent, last_dir)
        if dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            return dlg.result()
        return None


class CalibrationWindow(QtWidgets.QMainWindow):
    """Scatter of standards + fitted line + predicted unknowns, with stats."""

    def __init__(self, model, predictions=None, parent=None, last_dir: str = ""):
        super().__init__(parent)
        self.model = model
        self.predictions = list(predictions or [])
        self._last_dir = last_dir
        self._extra: list = []
        self.setWindowTitle(f"Calibration curve — {model.equation()}")
        self.resize(720, 600)

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        v = QtWidgets.QVBoxLayout(central)

        self.pw = pg.PlotWidget()
        v.addWidget(self.pw, 1)
        self._draw_plot()

        self.stats = QtWidgets.QLabel(self._stats_html())
        self.stats.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        self.stats.setWordWrap(True)
        v.addWidget(self.stats)

        bar = QtWidgets.QHBoxLayout()
        b_std = QtWidgets.QPushButton("Standards table…")
        b_std.clicked.connect(self._show_standards)
        b_pred = QtWidgets.QPushButton("Predictions table…")
        b_pred.clicked.connect(self._show_predictions)
        b_pred.setEnabled(bool(self.predictions))
        bar.addWidget(b_std)
        bar.addWidget(b_pred)
        bar.addStretch(1)
        v.addLayout(bar)

        tb = self.addToolBar("Calibration")
        tb.addAction(self.style().standardIcon(
            QtWidgets.QStyle.StandardPixmap.SP_DialogSaveButton), "Export image",
            self._export_image)
        tb.addAction("Export data (CSV)…", self._export_csv)
        self.status = self.statusBar()

    # ---- plot -----------------------------------------------------------
    def _unit_suffix(self) -> str:
        return f" {self.model.conc_unit}" if self.model.conc_unit else ""

    def _draw_plot(self) -> None:
        m = self.model
        plot = self.pw.getPlotItem()
        plot.addLegend(offset=(10, 10))
        unit = f" ({m.conc_unit})" if m.conc_unit else ""
        plot.setLabel("bottom", f"Concentration{unit}")
        plot.setLabel("left", m.signal_label)
        plot.showGrid(x=True, y=True, alpha=0.3)

        plot.plot(m.conc, m.signal, pen=None, symbol="o", symbolSize=9,
                  symbolBrush=pg.mkBrush("#1f77b4"), symbolPen=pg.mkPen("#0d3d66"),
                  name="standards")

        lo, hi = m.conc_range
        pad = 0.05 * ((hi - lo) or 1.0)
        finite_unk = [p for p in self.predictions if np.isfinite(p.get("conc", np.nan))]
        right = max([hi] + [p["conc"] for p in finite_unk]) + pad
        left = min([lo, 0.0] + [p["conc"] for p in finite_unk]) - pad
        xs = np.linspace(left, right, 200)
        plot.plot(xs, m.predict_signal(xs), pen=pg.mkPen("#d62728", width=2), name="fit")

        if finite_unk:
            cu = np.array([p["conc"] for p in finite_unk], float)
            su = np.array([p["signal"] for p in finite_unk], float)
            plot.plot(cu, su, pen=None, symbol="t", symbolSize=11,
                      symbolBrush=pg.mkBrush("#2ca02c"), symbolPen=pg.mkPen("#145214"),
                      name="unknowns")
            ci = np.array([p["ci"] if p.get("ci") else 0.0 for p in finite_unk], float)
            if np.any(ci > 0):
                err = pg.ErrorBarItem(x=cu, y=su, left=ci, right=ci,
                                      beam=0.02 * ((hi - lo) or 1.0),
                                      pen=pg.mkPen("#2ca02c"))
                plot.addItem(err)

    def _stats_html(self) -> str:
        m = self.model
        kind = ("quadratic" if m.degree == 2
                else "linear through origin" if m.through_origin else "linear")
        u = self._unit_suffix()
        rows = [f"<b>{m.equation()}</b>",
                f"Model: {kind} · n = {m.conc.size} standards · "
                f"signal = {m.signal_label}"]
        if m.degree == 1:
            rows.append(f"R² = {m.r_squared:.5f}&nbsp;&nbsp;(r = {m.r:.5f})")
            rows.append(f"slope = {m.slope:.6g} ± {m.se_slope:.3g}")
            if not m.through_origin:
                rows.append(f"intercept = {m.intercept:.6g} ± {m.se_intercept:.3g}")
            rows.append(f"s(y/x) = {m.s_yx:.4g}")
            if m.lod is not None:
                rows.append(f"LOD = {m.lod:.4g}{u}&nbsp;&nbsp;LOQ = {m.loq:.4g}{u}")
        else:
            rows.append(f"R² = {m.r_squared:.5f}")
        return "<br>".join(rows)

    # ---- tables ---------------------------------------------------------
    def _show_standards(self) -> None:
        m = self.model
        rows = [[i + 1, f"{c:.6g}", f"{s:.6g}", f"{f:.6g}", f"{r:+.4g}"]
                for i, (c, s, f, r) in enumerate(
                    zip(m.conc, m.signal, m.fitted, m.residuals))]
        d = TableDialog("Calibration standards",
                        ["#", f"Concentration{self._unit_suffix()}", "Signal",
                         "Fitted", "Residual"], rows, self,
                        summary=f"{m.equation()} · R² = {m.r_squared:.5f}",
                        default_dir=self._last_dir)
        d.show()
        self._extra.append(d)

    def _show_predictions(self) -> None:
        if not self.predictions:
            return
        u = self._unit_suffix()
        rows = []
        for i, p in enumerate(self.predictions):
            conc = p.get("conc", float("nan"))
            ci = p.get("ci")
            ci_txt = f"± {ci:.4g}" if ci else "—"
            flag = "" if p.get("in_range") else "  ⚠ out of range"
            rows.append([i + 1, p.get("name", ""), f"{p.get('signal', float('nan')):.6g}",
                         "—" if not np.isfinite(conc) else f"{conc:.6g}{u}", ci_txt, flag.strip()])
        d = TableDialog("Predicted concentrations",
                        ["#", "Sample", "Signal", "Concentration", "95% CI", "Note"],
                        rows, self,
                        summary="Unknowns predicted by inverse calibration.",
                        default_dir=self._last_dir)
        d.show()
        self._extra.append(d)

    # ---- export ---------------------------------------------------------
    def _export_image(self) -> None:
        default = os.path.join(self._last_dir, "calibration.png") if self._last_dir \
            else "calibration.png"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export image", default, "PNG (*.png)")
        if path:
            self.centralWidget().grab().save(path)
            self.status.showMessage(f"Saved {os.path.abspath(path)}", 6000)

    def _export_csv(self) -> None:
        m = self.model
        default = os.path.join(self._last_dir, "calibration.csv") if self._last_dir \
            else "calibration.csv"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export calibration data", default, "CSV (*.csv)")
        if not path:
            return
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["# equation", m.equation()])
            w.writerow(["# R_squared", f"{m.r_squared:.6g}"])
            w.writerow(["# s(y/x)", f"{m.s_yx:.6g}"])
            if m.degree == 1:
                w.writerow(["# slope", f"{m.slope:.6g}", "se", f"{m.se_slope:.6g}"])
                if not m.through_origin:
                    w.writerow(["# intercept", f"{m.intercept:.6g}", "se",
                                f"{m.se_intercept:.6g}"])
                if m.lod is not None:
                    w.writerow(["# LOD", f"{m.lod:.6g}", "LOQ", f"{m.loq:.6g}"])
            w.writerow([])
            w.writerow(["standard_concentration", "signal", "fitted", "residual"])
            for c, s, f, r in zip(m.conc, m.signal, m.fitted, m.residuals):
                w.writerow([f"{c:.8g}", f"{s:.8g}", f"{f:.8g}", f"{r:.8g}"])
            if self.predictions:
                w.writerow([])
                w.writerow(["unknown_sample", "signal", "predicted_concentration",
                            "ci_95", "in_range"])
                for p in self.predictions:
                    w.writerow([p.get("name", ""), f"{p.get('signal', float('nan')):.8g}",
                                f"{p.get('conc', float('nan')):.8g}",
                                "" if not p.get("ci") else f"{p['ci']:.8g}",
                                bool(p.get("in_range"))])
        self.status.showMessage(f"Saved {os.path.abspath(path)}", 6000)
