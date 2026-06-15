"""Headless (offscreen) GUI smoke test: builds the window and drives it.

Run:  python tests/gui_smoke_test.py
Uses the Qt 'offscreen' platform so it needs no display.
"""
import os
import sys
import tempfile

import numpy as np

os.environ["QT_QPA_PLATFORM"] = "offscreen"
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from PySide6 import QtCore, QtGui, QtWidgets  # noqa: E402

from specview.app import configure_pyqtgraph  # noqa: E402
from specview.ui import dialogs as dlg  # noqa: E402

PASS, FAIL = 0, 0


def check(name, cond):
    global PASS, FAIL
    PASS, FAIL = PASS + bool(cond), FAIL + (not cond)
    print(f"  {'ok  ' if cond else 'FAIL'} {name}")


# --- neutralise modal dialogs so nothing blocks -----------------------------
dlg.FormDialog.exec_form = staticmethod(
    lambda title, fields, parent=None, description=None:
    {f["key"]: f.get("default") for f in fields})
QtWidgets.QInputDialog.getDouble = staticmethod(lambda *a, **k: (785.0, True))
QtWidgets.QColorDialog.getColor = staticmethod(lambda *a, **k: QtGui.QColor("#ff0000"))
for _m in ("warning", "information", "critical", "about"):
    setattr(QtWidgets.QMessageBox, _m, staticmethod(lambda *a, **k: None))
QtWidgets.QMessageBox.question = staticmethod(
    lambda *a, **k: QtWidgets.QMessageBox.StandardButton.Yes)
QtWidgets.QMessageBox.exec = lambda self: 0          # _report_saved info box
QtWidgets.QMessageBox.clickedButton = lambda self: None

configure_pyqtgraph()
app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

from specview.ui import MainWindow  # noqa: E402

win = MainWindow()
win.show()

print("== window & demo ==")
win.load_demo()
check("5 spectra loaded into table", win.table.rowCount() == 5)
check("plot has 5 curves", len(win.plotview._curves) == 5)


def select_rows(*rows):
    sm = win.table.selectionModel()
    sm.clearSelection()
    for r in rows:
        sm.select(win.table.model().index(r, 0),
                  QtCore.QItemSelectionModel.SelectionFlag.Select
                  | QtCore.QItemSelectionModel.SelectionFlag.Rows)


print("== processing actions (no-arg) ==")
select_rows(0)
win.baseline_rb()
check("rubberband applied, still 5 spectra", len(win.document) == 5)
win.snv()
check("SNV applied", len(win.document) == 5)
win.normalize("max")
check("normalize max applied", abs(max(abs(win.document[0].y.min()),
                                       abs(win.document[0].y.max())) - 1.0) < 1e-6)

print("== dialog-driven actions ==")
select_rows(0)
n_before = win.document[0].npoints
win.smooth_sg()
check("Savitzky-Golay ran", win.document[0].npoints == n_before)
win.derivative()
check("derivative ran", len(win.document) == 5)
win.baseline_als()
check("ALS ran", len(win.document) == 5)
win.baseline_airpls()
check("airPLS ran", len(win.document) == 5)

print("== axis conversions ==")
win.remove_all()
win.load_demo()           # fresh [FTIR, Raman, UV/Vis, NIR]
select_rows(0)            # FTIR (cm-1, absorbance)
win.convert_x("nm")
check("FTIR converted to nm", win.document[0].x_unit == "nm")
win.convert_y("transmittance")
check("FTIR converted to transmittance", win.document[0].y_unit == "transmittance")

print("== arithmetic & average ==")
win.remove_all()
win.load_demo()           # fresh set; rows 0 (FTIR) and 1 (Raman) overlap 400–3200
n0 = len(win.document)
select_rows(0, 1)
win.arithmetic("sub")
check("arithmetic added one spectrum", len(win.document) == n0 + 1)
select_rows(0, 1)
win.average("mean")
check("average added one spectrum", len(win.document) == n0 + 2)

print("== undo / redo ==")
n_now = len(win.document)
win.undo()
check("undo removed the average", len(win.document) == n_now - 1)
win.redo()
check("redo restored it", len(win.document) == n_now)

print("== display toggles & export ==")
win.plotview.set_stack_offset(0.5)
win.plotview.set_flip_x(True)
win.set_dark(True)
check("stack offset set", win.plotview.stack_offset == 0.5)
tmp = tempfile.mkdtemp()
png = os.path.join(tmp, "plot.png")
win.plotview.export_image(png)
check("PNG export created file", os.path.getsize(png) > 1000)
svg = os.path.join(tmp, "plot.svg")
win.plotview.export_image(svg)
check("SVG export created file", os.path.getsize(svg) > 200)

print("== remove ==")
select_rows(0)
m = len(win.document)
win.remove_selected()
check("remove_selected dropped one", len(win.document) == m - 1)
win.remove_all()
check("remove_all clears document", len(win.document) == 0)

print("== real FormDialog construction (every field type) ==")
# Build the ACTUAL dialog (not the monkeypatched exec_form) so __init__ runs for
# real — this catches binding bugs like QtWidgets.Qt vs QtCore.Qt.
_real = dlg.FormDialog("t", [
    {"key": "i", "label": "i", "type": "int", "default": 11, "min": 3, "max": 99},
    {"key": "f", "label": "f", "type": "float", "default": 1.5, "min": 0, "max": 9},
    {"key": "c", "label": "c", "type": "choice",
     "options": [("A", "a"), ("B", "b")], "default": "b"},
    {"key": "b", "label": "b", "type": "bool", "default": True},
    {"key": "t", "label": "t", "type": "text", "default": "hi"},
])
check("FormDialog builds every field type",
      _real.values() == {"i": 11, "f": 1.5, "c": "b", "b": True, "t": "hi"})
check("Export-combined menu action exists",
      win.act_export_data.text().replace("&", "") == "Export combined data…")

print("== combined data export (handler) ==")
from specview.spectrum import Spectrum  # noqa: E402
from specview.formats import load_any  # noqa: E402

out_csv = os.path.join(tempfile.mkdtemp(), "combo.csv")
QtWidgets.QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: (out_csv, "CSV (*.csv)"))
win.remove_all()
for k in range(3):
    win.document.add(Spectrum(np.linspace(1000, 2000, 256),
                              np.linspace(0.1, 1.0, 256) + k * 0.05,
                              name=f"s{k}", x_unit="nm", y_unit="absorbance"))
win._rebuild_table()
win.table.clearSelection()
win.export_combined()                                # FormDialog/getSaveFileName patched
check("combined export via handler wrote a file",
      os.path.exists(out_csv) and os.path.getsize(out_csv) > 100)
check("combined export re-imports 3 spectra", len(load_any(out_csv)) == 3)
check("last_dir remembered after save", win._last_dir == os.path.dirname(out_csv))

print("== analysis handlers (find / fit / integrate / clear) ==")
win.remove_all()
xa = np.linspace(400, 800, 1500)
ya = np.exp(-(xa - 500) ** 2 / (2 * 8 ** 2)) + 0.6 * np.exp(-(xa - 650) ** 2 / (2 * 12 ** 2))
win.document.add(Spectrum(xa, ya, name="syn", x_unit="nm", y_unit="absorbance"))
win._rebuild_table()
select_rows(0)
win.find_peaks()
check("find_peaks drew markers", len(win.plotview._peak_items) > 0)
win.fit_peaks()
check("fit drew overlays", len(win.plotview._fit_items) > 0)
# the fit dialog's "Add components to list" button callback
fit_dlg = win._dialogs[-1]
add_btn = [b for b in fit_dlg.findChildren(QtWidgets.QPushButton)
           if b.text() == "Add components to list"][0]
n_before = len(win.document)
add_btn.click()
check("'Add components to list' adds spectra", len(win.document) > n_before)
win.integrate_range()
check("integrate drew a region", len(win.plotview._fit_items) > 0)
win.clear_analysis()
check("clear_analysis removes all overlays",
      not win.plotview._peak_items and not win.plotview._fit_items)

print("== peak label staggering (no overlap for close peaks) ==")
import pyqtgraph as pg  # noqa: E402
from specview.analysis import Peak  # noqa: E402

win.plotview.plot.getViewBox().setRange(xRange=(0, 3800), yRange=(0, 40000), padding=0)
win.plotview.mark_peaks([Peak(1052, 13000, 10, 1, 1, 0), Peak(1093, 12500, 10, 1, 1, 1)])
labels = [it for it in win.plotview._peak_items if isinstance(it, pg.TextItem)]
ys = sorted(it.pos().y() for it in labels)
check("two close peaks -> 2 labels", len(labels) == 2)
check("close-peak labels staggered vertically", ys[1] - ys[0] > 500)
win.plotview.mark_peaks([Peak(500, 8000, 10, 1, 1, 0), Peak(3000, 8000, 10, 1, 1, 1)])
ys2 = sorted(it.pos().y() for it in win.plotview._peak_items
             if isinstance(it, pg.TextItem))
check("far-apart labels share a tier (not staggered)", abs(ys2[1] - ys2[0]) < 1e-6)

print("== library / mixture / XRF handlers ==")
from specview.demo import demo_xrf  # noqa: E402

win.remove_all()
xl = np.linspace(400, 1800, 600)
def _bb(c, a, w):  # noqa: E306
    return a * np.exp(-(xl - c) ** 2 / (2 * w ** 2))
rA = Spectrum(xl, _bb(700, 1, 30) + _bb(1200, .5, 40), name="refA", x_unit="cm-1")
rB = Spectrum(xl, _bb(900, .8, 25) + _bb(1500, .6, 35), name="refB", x_unit="cm-1")
mix = Spectrum(xl, 0.7 * rA.y + 0.3 * rB.y, name="mix", x_unit="cm-1")
for s in (mix, rA, rB):
    win.document.add(s)
win._rebuild_table()
select_rows(1, 2)
win.library_add()
check("library_add stored 2 references", len(win.library) == 2)
select_rows(0)
win.library_search()
check("library_search opened a results dialog",
      win._dialogs and isinstance(win._dialogs[-1], dlg.TableDialog))
select_rows(0, 1, 2)
win.mixture_analysis()
check("mixture_analysis drew the fit overlay", len(win.plotview._fit_items) > 0)

win.remove_all()
win.document.add(demo_xrf())
win._rebuild_table()
select_rows(0)
win.identify_xrf()
check("identify_xrf marked element peaks", len(win.plotview._peak_items) > 0)

lib_path = os.path.join(tempfile.mkdtemp(), "lib.speclib")
QtWidgets.QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: (lib_path, ""))
QtWidgets.QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: (lib_path, ""))
win.library_save()
win.library_clear()
check("library_clear empties the library", len(win.library) == 0)
win.library_load()
check("library_load restores entries from file", len(win.library) == 2)

print("== 2D-COS / EEM windows ==")
from specview.ui.mapwindow import EEMWindow, MapWindow, Surface3DDialog  # noqa: E402
from specview.demo import demo_eem  # noqa: E402

win.remove_all()
win.load_demo_series()
sm2 = win.table.selectionModel()
for r in range(len(win.document)):
    sm2.select(win.table.model().index(r, 0),
               QtCore.QItemSelectionModel.SelectionFlag.Select
               | QtCore.QItemSelectionModel.SelectionFlag.Rows)
win.cos2d_analysis()
maps = [d for d in win._dialogs if isinstance(d, MapWindow)]
check("2D-COS opened a 2-panel map window", bool(maps) and len(maps[-1]._items) == 2)
win.open_demo_eem()
eemw = [d for d in win._dialogs if isinstance(d, EEMWindow)]
check("EEM window opened", bool(eemw))
eemw[-1]._toggle_scatter(True)
check("EEM scatter toggle masks NaN", bool(np.isnan(eemw[-1]._items[0]["Zxy"]).any()))
check("EEM 3D surface dialog builds", Surface3DDialog(demo_eem()) is not None)

print(f"\n{PASS} passed, {FAIL} failed")
app.quit()
sys.exit(1 if FAIL else 0)
