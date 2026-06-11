"""The central plotting widget, built on pyqtgraph."""
from __future__ import annotations

import os

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

from ..spectrum import Spectrum, SpectrumSet


class PlotView(QtWidgets.QWidget):
    """Displays a :class:`SpectrumSet` with overlay, stacking and a crosshair."""

    #: emitted on mouse move with the cursor position in data coordinates
    cursorMoved = QtCore.Signal(float, float)

    def __init__(self, document: SpectrumSet, parent=None):
        super().__init__(parent)
        self.document = document

        # display state
        self.stack_offset = 0.0      # fraction of global y-span per spectrum
        self.flip_x = False
        self.show_grid = True
        self.log_y = False

        self.plot_widget = pg.PlotWidget()
        self.plot = self.plot_widget.getPlotItem()
        self.plot.showGrid(x=True, y=True, alpha=0.25)
        self.plot.setMenuEnabled(True)
        self.legend = self.plot.addLegend(offset=(-10, 10))

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.plot_widget)

        self._curves: list[pg.PlotDataItem] = []
        self._peak_items: list = []      # scatter + text markers from find_peaks
        self._fit_items: list = []       # fitted sum + component overlays

        # crosshair
        pen = pg.mkPen(color=(120, 120, 120), width=1, style=QtCore.Qt.PenStyle.DashLine)
        self.vline = pg.InfiniteLine(angle=90, movable=False, pen=pen)
        self.hline = pg.InfiniteLine(angle=0, movable=False, pen=pen)
        self.plot.addItem(self.vline, ignoreBounds=True)
        self.plot.addItem(self.hline, ignoreBounds=True)
        self._proxy = pg.SignalProxy(self.plot.scene().sigMouseMoved,
                                     rateLimit=60, slot=self._on_mouse_move)

        self.refresh()

    # ---- drawing ---------------------------------------------------------
    def _global_span(self) -> float:
        spans = [s.yrange[1] - s.yrange[0] for s in self.document.visible_spectra()
                 if s.npoints]
        return max(spans) if spans else 1.0

    def refresh(self) -> None:
        """Redraw all visible spectra according to the current display state."""
        for c in self._curves:
            self.plot.removeItem(c)
        self._curves.clear()
        self.legend.clear()

        span = self._global_span()
        visible = self.document.visible_spectra()
        x_label = y_label = ""
        for i, spec in enumerate(visible):
            if spec.npoints == 0:
                continue
            offset = i * self.stack_offset * span
            pen = pg.mkPen(color=spec.color or "#1f77b4",
                           width=float(spec.meta.get("line_width", 1.4)))
            curve = self.plot.plot(spec.x, spec.y + offset, pen=pen, name=spec.name,
                                   antialias=True)
            self._curves.append(curve)
            x_label, y_label = spec.x_label, spec.y_label

        self.plot.setLabel("bottom", x_label or "x")
        self.plot.setLabel("left", y_label or "y")
        self.plot.showGrid(x=self.show_grid, y=self.show_grid, alpha=0.25)
        self.plot.getViewBox().invertX(self.flip_x)
        self.plot.setLogMode(x=False, y=self.log_y)

    def autoscale(self) -> None:
        self.plot.enableAutoRange()
        self.plot.autoRange()

    # ---- display toggles -------------------------------------------------
    def set_stack_offset(self, frac: float) -> None:
        self.stack_offset = max(0.0, float(frac))
        self.refresh()

    def set_flip_x(self, flip: bool) -> None:
        self.flip_x = bool(flip)
        self.plot.getViewBox().invertX(self.flip_x)

    def set_grid(self, on: bool) -> None:
        self.show_grid = bool(on)
        self.plot.showGrid(x=on, y=on, alpha=0.25)

    def set_log_y(self, on: bool) -> None:
        self.log_y = bool(on)
        self.plot.setLogMode(x=False, y=on)

    def suggest_flip(self) -> None:
        """Flip X automatically for IR/Raman where high→low is conventional."""
        vis = self.document.visible_spectra()
        if vis and vis[0].x_unit in ("cm-1", "raman_cm-1"):
            self.set_flip_x(True)

    # ---- analysis overlays ----------------------------------------------
    def mark_peaks(self, peaks, color: str = "#d62728") -> None:
        """Draw markers + position labels, staggering labels of close peaks.

        Labels of peaks closer than ~4% of the visible x-range are stacked into
        an ascending staircase so their text never overlaps horizontally.
        """
        self.clear_peaks()
        if not peaks:
            return
        xs = [p.center for p in peaks]
        ys = [p.height for p in peaks]
        scatter = pg.ScatterPlotItem(xs, ys, symbol="t1", size=11,
                                     brush=pg.mkBrush(color), pen=pg.mkPen("k", width=0.5))
        self.plot.addItem(scatter)
        self._peak_items.append(scatter)

        (vx0, vx1), (vy0, vy1) = self.plot.getViewBox().viewRange()
        x_span = (vx1 - vx0) or 1.0
        y_span = (vy1 - vy0) or 1.0
        near = 0.04 * x_span      # peaks closer than this share a cluster
        step = 0.05 * y_span      # one vertical stagger tier
        last_x, level = None, 0
        for p in sorted(peaks, key=lambda q: q.center):
            level = level + 1 if (last_x is not None and abs(p.center - last_x) < near) \
                else 0
            last_x = p.center
            txt = pg.TextItem(f"{p.center:.4g}", color=color, anchor=(0.5, 1.0))
            txt.setPos(p.center, p.height + step * (0.3 + level))
            self.plot.addItem(txt)
            self._peak_items.append(txt)

    def show_fit(self, x, total, comp_curves, baseline=None) -> None:
        """Overlay the fitted sum (bold dashed) and each component (thin dashed)."""
        self.clear_fit()
        for c in comp_curves:
            item = self.plot.plot(x, c, pen=pg.mkPen("#7f7f7f", width=1,
                                  style=QtCore.Qt.PenStyle.DashLine))
            self._fit_items.append(item)
        if baseline is not None:
            item = self.plot.plot(x, baseline, pen=pg.mkPen("#bcbd22", width=1,
                                  style=QtCore.Qt.PenStyle.DotLine))
            self._fit_items.append(item)
        total_item = self.plot.plot(x, total, pen=pg.mkPen("#000000", width=2,
                                    style=QtCore.Qt.PenStyle.DashLine))
        self._fit_items.append(total_item)

    def show_region(self, x_lo: float, x_hi: float) -> None:
        """Shade an integrated x-region (re-uses the fit-overlay slot)."""
        self.clear_fit()
        region = pg.LinearRegionItem(values=(x_lo, x_hi), movable=False,
                                     brush=pg.mkBrush(31, 119, 180, 50))
        region.setZValue(-10)
        self.plot.addItem(region)
        self._fit_items.append(region)

    def clear_peaks(self) -> None:
        for it in self._peak_items:
            self.plot.removeItem(it)
        self._peak_items.clear()

    def clear_fit(self) -> None:
        for it in self._fit_items:
            self.plot.removeItem(it)
        self._fit_items.clear()

    def clear_analysis(self) -> None:
        self.clear_peaks()
        self.clear_fit()

    # ---- crosshair -------------------------------------------------------
    def _on_mouse_move(self, evt) -> None:
        pos = evt[0]
        if not self.plot.sceneBoundingRect().contains(pos):
            return
        mp = self.plot.getViewBox().mapSceneToView(pos)
        x, y = mp.x(), mp.y()
        self.vline.setPos(x)
        self.hline.setPos(y)
        self.cursorMoved.emit(x, y)

    # ---- export ----------------------------------------------------------
    def export_image(self, path: str, width: int = 1920) -> None:
        """Export the current view.

        Raster formats (png/jpg/bmp/tif) use pyqtgraph's WYSIWYG exporter;
        vector formats (svg/pdf/eps) are rendered with matplotlib, which is
        far more robust than pyqtgraph's SVG writer and also yields PDF/EPS.
        """
        ext = os.path.splitext(path)[1].lower()
        if ext in (".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff"):
            import pyqtgraph.exporters as exporters
            exp = exporters.ImageExporter(self.plot)
            exp.parameters()["width"] = int(width)
            exp.export(path)
        else:
            self._export_vector(path)

    def _export_vector(self, path: str) -> None:
        import matplotlib
        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
        span = self._global_span()
        visible = self.document.visible_spectra()
        x_label = y_label = ""
        for i, spec in enumerate(visible):
            if spec.npoints == 0:
                continue
            offset = i * self.stack_offset * span
            ax.plot(spec.x, spec.y + offset, color=spec.color or "#1f77b4",
                    lw=float(spec.meta.get("line_width", 1.2)), label=spec.name)
            x_label, y_label = spec.x_label, spec.y_label
        ax.set_xlabel(x_label or "x")
        ax.set_ylabel(y_label or "y")
        if self.flip_x:
            ax.invert_xaxis()
        if self.show_grid:
            ax.grid(alpha=0.3)
        if self.log_y:
            ax.set_yscale("log")
        if visible:
            ax.legend(fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
