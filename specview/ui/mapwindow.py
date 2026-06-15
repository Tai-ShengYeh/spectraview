"""Windows for 2-D maps: heatmaps (EEM, 2D-COS) and a 3-D EEM surface."""
from __future__ import annotations

import os

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

from .. import eem as eem_mod


def _get_cmap(name: str):
    for kwargs in ({}, {"source": "matplotlib"}):
        try:
            return pg.colormap.get(name, **kwargs)
        except Exception:  # noqa: BLE001
            continue
    return None


class MapWindow(QtWidgets.QMainWindow):
    """Generic multi-panel 2-D heatmap viewer (used by 2D-COS and EEM)."""

    def __init__(self, panels: list[dict], title: str = "2D map", parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(460 * len(panels) + 120, 540)
        self.glw = pg.GraphicsLayoutWidget()
        self.setCentralWidget(self.glw)
        self._items: list[dict] = []
        for p in panels:
            self._add_panel(p)
        self.status = self.statusBar()
        self._readout = QtWidgets.QLabel("  ")
        self.status.addPermanentWidget(self._readout)
        self._proxy = pg.SignalProxy(self.glw.scene().sigMouseMoved, rateLimit=50,
                                     slot=self._on_move)
        self._build_toolbar()

    def _add_panel(self, p: dict) -> None:
        plot = self.glw.addPlot(title=p.get("label", ""))
        plot.setLabel("bottom", p.get("xlabel", "x"))
        plot.setLabel("left", p.get("ylabel", "y"))
        x = np.asarray(p["x"], float)
        y = np.asarray(p["y"], float)
        # Z is passed in ImageItem-native orientation: Z[x_index, y_index].
        Zxy = np.asarray(p["Z"], float)
        img = pg.ImageItem(Zxy)
        plot.addItem(img)
        img.setRect(QtCore.QRectF(float(x[0]), float(y[0]),
                                  float(x[-1] - x[0]), float(y[-1] - y[0])))
        cmap = _get_cmap("CET-D1" if p.get("diverging") else p.get("cmap", "viridis"))
        finite = Zxy[np.isfinite(Zxy)]
        if p.get("diverging"):
            a = float(np.nanmax(np.abs(finite))) if finite.size else 1.0
            levels = (-a, a)
        else:
            levels = (float(finite.min()), float(finite.max())) if finite.size else (0.0, 1.0)
        if cmap is not None:
            img.setColorMap(cmap)
        img.setLevels(levels)
        try:
            cb = pg.ColorBarItem(values=levels, colorMap=cmap)
            cb.setImageItem(img, insert_in=plot)
        except Exception:  # noqa: BLE001
            pass
        pen = pg.mkPen(180, 180, 180, style=QtCore.Qt.PenStyle.DashLine)
        vL = pg.InfiniteLine(angle=90, pen=pen)
        hL = pg.InfiniteLine(angle=0, pen=pen)
        plot.addItem(vL, ignoreBounds=True)
        plot.addItem(hL, ignoreBounds=True)
        self._items.append({"plot": plot, "img": img, "x": x, "y": y, "Zxy": Zxy,
                            "vL": vL, "hL": hL, "label": p.get("label", "")})

    def update_panel(self, i: int, x, y, Z) -> None:
        it = self._items[i]
        Zxy = np.asarray(Z, float)        # native orientation Z[x_index, y_index]
        it["Zxy"], it["x"], it["y"] = Zxy, np.asarray(x, float), np.asarray(y, float)
        it["img"].setImage(Zxy)
        finite = Zxy[np.isfinite(Zxy)]
        if finite.size:
            it["img"].setLevels((float(finite.min()), float(finite.max())))

    def _on_move(self, evt) -> None:
        pos = evt[0]
        for it in self._items:
            if it["plot"].sceneBoundingRect().contains(pos):
                mp = it["plot"].getViewBox().mapSceneToView(pos)
                xx, yy = mp.x(), mp.y()
                it["vL"].setPos(xx)
                it["hL"].setPos(yy)
                xi = int(np.clip(np.searchsorted(it["x"], xx), 0, it["x"].size - 1))
                yi = int(np.clip(np.searchsorted(it["y"], yy), 0, it["y"].size - 1))
                zxy = it["Zxy"]
                val = zxy[min(xi, zxy.shape[0] - 1), min(yi, zxy.shape[1] - 1)]
                self._readout.setText(f"{it['label']}   x={xx:.4g}   y={yy:.4g}   z={val:.4g}")
                return

    def _build_toolbar(self) -> None:
        tb = self.addToolBar("Map")
        tb.addAction(QtWidgets.QApplication.style().standardIcon(
            QtWidgets.QStyle.StandardPixmap.SP_DialogSaveButton), "Export image",
            self._export_image)
        tb.addAction("Export data…", self._export_data)

    def _export_image(self) -> None:
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export image",
                                                        "map.png", "PNG (*.png)")
        if path:
            self.glw.grab().save(path)
            self.status.showMessage(f"Saved {os.path.abspath(path)}", 6000)

    def _export_data(self) -> None:
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export matrix",
                                                        "matrix.csv", "CSV (*.csv)")
        if not path:
            return
        it = self._items[0]
        x, y, Zxy = it["x"], it["y"], it["Zxy"]      # Zxy is (nx, ny)
        import csv
        with open(path, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh)
            w.writerow([""] + [f"{v:.6g}" for v in x])
            for j in range(y.size):
                w.writerow([f"{y[j]:.6g}"] + [f"{Zxy[i, j]:.6g}" for i in range(x.size)])
        self.status.showMessage(f"Saved {os.path.abspath(path)}", 6000)


class EEMWindow(MapWindow):
    """EEM heatmap with 3-D surface and scatter-removal controls."""

    def __init__(self, eem, parent=None):
        self.eem_orig = eem
        self.eem = eem
        super().__init__([self._panel(eem)], title=f"EEM — {eem.name}", parent=parent)

    @staticmethod
    def _panel(e) -> dict:
        # ImageItem-native: Z[em_index, ex_index] -> e.Z is (ex, em) so transpose.
        return {"x": e.em, "y": e.ex, "Z": e.Z.T, "label": e.name,
                "xlabel": "Emission (nm)", "ylabel": "Excitation (nm)", "cmap": "viridis"}

    def _build_toolbar(self) -> None:
        super()._build_toolbar()
        tb = self.addToolBar("EEM")
        tb.addAction("3D surface →", self._open_3d)
        self._act_scatter = tb.addAction("Remove scatter")
        self._act_scatter.setCheckable(True)
        self._act_scatter.toggled.connect(self._toggle_scatter)

    def _toggle_scatter(self, on: bool) -> None:
        self.eem = eem_mod.remove_scatter(self.eem_orig) if on else self.eem_orig
        self.update_panel(0, self.eem.em, self.eem.ex, self.eem.Z.T)
        self.status.showMessage("Rayleigh scatter masked." if on else "Scatter restored.", 4000)

    def _open_3d(self) -> None:
        dlg = Surface3DDialog(self.eem, self)
        dlg.show()
        self._dlg3d = dlg


class Surface3DDialog(QtWidgets.QDialog):
    """A rotatable 3-D EEM surface (matplotlib embedded in Qt)."""

    def __init__(self, eem, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"EEM 3D surface — {eem.name}")
        self.resize(740, 580)
        os.environ.setdefault("QT_API", "pyside6")
        import matplotlib
        matplotlib.use("QtAgg", force=False)
        from matplotlib.backends.backend_qtagg import (FigureCanvasQTAgg,
                                                       NavigationToolbar2QT)
        from matplotlib.figure import Figure

        fig = Figure(figsize=(7, 5.2))
        ax = fig.add_subplot(111, projection="3d")
        EM, EX = np.meshgrid(eem.em, eem.ex)       # (n_ex, n_em)
        Z = eem.Z.copy()
        if np.isnan(Z).any():
            Z = np.nan_to_num(Z, nan=float(np.nanmin(Z)))
        ax.plot_surface(EM, EX, Z, cmap="viridis", linewidth=0, antialiased=True)
        ax.set_xlabel("Emission (nm)")
        ax.set_ylabel("Excitation (nm)")
        ax.set_zlabel("Intensity")
        ax.set_title(eem.name)
        canvas = FigureCanvasQTAgg(fig)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(NavigationToolbar2QT(canvas, self))
        layout.addWidget(canvas)
