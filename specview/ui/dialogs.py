"""Reusable parameter dialogs for processing operations."""
from __future__ import annotations

from PySide6 import QtCore, QtWidgets


class FormDialog(QtWidgets.QDialog):
    """A small modal form built from a list of field specifications.

    Each field is a dict:
        {"key", "label", "type": int|float|choice|text|bool,
         "default", "min", "max", "decimals", "options"}
    ``exec_form`` returns a dict of values, or ``None`` if cancelled.
    """

    def __init__(self, title: str, fields: list[dict], parent=None,
                 description: str | None = None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self._fields = fields
        self._widgets: dict[str, QtWidgets.QWidget] = {}

        layout = QtWidgets.QVBoxLayout(self)
        if description:
            lbl = QtWidgets.QLabel(description)
            lbl.setWordWrap(True)
            lbl.setStyleSheet("color: palette(mid); margin-bottom: 6px;")
            layout.addWidget(lbl)

        form = QtWidgets.QFormLayout()
        form.setLabelAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
        for f in fields:
            w = self._make_widget(f)
            self._widgets[f["key"]] = w
            form.addRow(f["label"], w)
        layout.addLayout(form)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _make_widget(self, f: dict) -> QtWidgets.QWidget:
        t = f.get("type", "float")
        if t == "int":
            w = QtWidgets.QSpinBox()
            w.setRange(int(f.get("min", 0)), int(f.get("max", 10_000_000)))
            w.setValue(int(f.get("default", 0)))
            w.setSingleStep(int(f.get("step", 1)))
        elif t == "float":
            w = QtWidgets.QDoubleSpinBox()
            w.setRange(float(f.get("min", -1e12)), float(f.get("max", 1e12)))
            w.setDecimals(int(f.get("decimals", 4)))
            w.setValue(float(f.get("default", 0.0)))
            w.setSingleStep(float(f.get("step", 0.1)))
        elif t == "choice":
            w = QtWidgets.QComboBox()
            for opt in f["options"]:
                # options may be (label, value) tuples or plain strings
                if isinstance(opt, (tuple, list)):
                    w.addItem(opt[0], opt[1])
                else:
                    w.addItem(str(opt), opt)
            if "default" in f:
                idx = w.findData(f["default"])
                if idx >= 0:
                    w.setCurrentIndex(idx)
        elif t == "bool":
            w = QtWidgets.QCheckBox()
            w.setChecked(bool(f.get("default", False)))
        else:  # text
            w = QtWidgets.QLineEdit(str(f.get("default", "")))
        return w

    def _value(self, key: str):
        w = self._widgets[key]
        if isinstance(w, QtWidgets.QComboBox):
            return w.currentData()
        if isinstance(w, (QtWidgets.QSpinBox, QtWidgets.QDoubleSpinBox)):
            return w.value()
        if isinstance(w, QtWidgets.QCheckBox):
            return w.isChecked()
        return w.text()

    def values(self) -> dict:
        return {f["key"]: self._value(f["key"]) for f in self._fields}

    @staticmethod
    def exec_form(title: str, fields: list[dict], parent=None,
                  description: str | None = None) -> dict | None:
        dlg = FormDialog(title, fields, parent, description)
        if dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            return dlg.values()
        return None


class TableDialog(QtWidgets.QDialog):
    """A non-modal results table with an 'Export CSV' button.

    ``extra_buttons`` is a list of (label, callback) added next to Export/Close.
    """

    def __init__(self, title: str, headers: list[str], rows: list[list],
                 parent=None, summary: str | None = None,
                 extra_buttons=None, default_dir: str = ""):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(560, 420)
        self._headers = headers
        self._rows = rows
        self._default_dir = default_dir

        layout = QtWidgets.QVBoxLayout(self)
        if summary:
            lbl = QtWidgets.QLabel(summary)
            lbl.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
            layout.addWidget(lbl)

        self.table = QtWidgets.QTableWidget(len(rows), len(headers))
        self.table.setHorizontalHeaderLabels(headers)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        for r, row in enumerate(rows):
            for c, val in enumerate(row):
                self.table.setItem(r, c, QtWidgets.QTableWidgetItem(str(val)))
        self.table.resizeColumnsToContents()
        layout.addWidget(self.table)

        bar = QtWidgets.QHBoxLayout()
        export = QtWidgets.QPushButton("Export CSV…")
        export.clicked.connect(self._export_csv)
        bar.addWidget(export)
        for label, callback in (extra_buttons or []):
            btn = QtWidgets.QPushButton(label)
            btn.clicked.connect(callback)
            bar.addWidget(btn)
        bar.addStretch(1)
        close = QtWidgets.QPushButton("Close")
        close.clicked.connect(self.accept)
        bar.addWidget(close)
        layout.addLayout(bar)

    def _export_csv(self) -> None:
        import csv
        import os
        default = os.path.join(self._default_dir, "table.csv") if self._default_dir \
            else "table.csv"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export table", default, "CSV (*.csv)")
        if not path:
            return
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(self._headers)
            w.writerows(self._rows)
        QtWidgets.QMessageBox.information(
            self, "Exported", f"Saved to:\n{os.path.abspath(path)}")
