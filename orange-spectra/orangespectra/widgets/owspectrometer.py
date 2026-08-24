"""Spectrometer — turn a diffraction-grating spectrum photo into a spectrum.

A Theremino-style spectrometer working on a captured image: rotate the photo so
the dispersion axis is horizontal, take a horizontal strip, read an intensity
profile along the columns, and calibrate pixel -> wavelength from known
reference lines. Outputs an Orange spectrum Table.
"""
import os

import matplotlib
import numpy as np
from AnyQt.QtWidgets import QFileDialog, QSlider

matplotlib.use("Qt5Agg")
# The Qt backend must be imported *after* matplotlib.use(); do not let an
# import sorter hoist the lines below above it.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from Orange.data import Table
from Orange.widgets import gui, settings
from Orange.widgets.widget import Msg, Output, OWWidget

from .. import mplfonts  # noqa: E402, F401  (CJK-capable preview fonts)
from ..spectrometer import (CAL_MODELS, CHANNELS, FLUORESCENT_LINES, ROTATIONS,
                            extract_profile, find_profile_peaks,
                            image_to_spectrum, load_rgb, pixel_to_wavelength,
                            rotate_rgb)
from ..table_io import table_from_spectra
from ._help import add_help


class OWSpectrometer(OWWidget):
    name = "Spectrometer"
    description = ("Turn a diffraction-grating spectrum photo (webcam "
                   "spectrometer) into a calibrated intensity-vs-wavelength "
                   "spectrum.")
    icon = "icons/spectrometer.svg"
    priority = 14
    keywords = ["spectrometer", "webcam", "grating", "theremino", "camera",
                "calibration", "peak", "rotate", "分光", "光柵", "旋轉", "峰"]

    class Outputs:
        spectrum = Output("Spectra", Table)

    class Error(OWWidget.Error):
        load_failed = Msg("{}")
        bad_calibration = Msg("Calibration: {}")

    class Warning(OWWidget.Warning):
        exact_fit = Msg("{n} calibration lines for a degree-{d} fit: R² is 1 by "
                        "construction, not a quality check. Add more lines.")
        folded = Msg("The fit turns over inside the image; {n} pixel(s) past "
                     "the turning point were dropped so the spectrum does not "
                     "fold onto itself.")

    class Information(OWWidget.Information):
        no_image = Msg("Choose a spectrum photo to read.")

    image_path: str = settings.Setting("")
    channel_idx: int = settings.Setting(0)
    rotate_idx: int = settings.Setting(0)              # index into ROTATIONS
    row_center_pct: int = settings.Setting(50)         # % of height
    row_frac_pct: int = settings.Setting(15)           # band height %
    flip: bool = settings.Setting(False)
    model_idx: int = settings.Setting(0)
    cal_text: str = settings.Setting("")               # "pixel=nm" per line/comma
    zoom_cal: bool = settings.Setting(True)            # plot only near the fit
    cursor_px: int = settings.Setting(0)               # peak cursor, in pixels
    peak_prom_pct: int = settings.Setting(5)           # % of full profile range
    peak_dist: int = settings.Setting(8)               # min peak separation, px
    peak_smooth: int = settings.Setting(5)             # smoothing before peak find
    assign_nm: str = settings.Setting("546.1")         # nm to attach to the cursor
    want_main_area = True

    def __init__(self):
        super().__init__()
        self._rgb = None            # the photo as loaded (unrotated)
        self._profile = None        # ROI profile, after rotate + flip
        self._peaks = []            # dicts from find_profile_peaks
        self._coeffs = None         # current pixel -> nm polynomial, or None
        self._cal_span = (None, None)   # pixel range the calibration covers
        self._spec_x = None         # x values actually plotted (nm or pixel)
        self._cursor_artists = []

        add_help(self,
                 "選一張『相機＋繞射光柵』拍到的光譜照片 → 需要時先旋轉，讓色散方向"
                 "變成水平 → 取一條水平帶、沿水平方向讀強度 → 用已知譜線把像素校準成"
                 "波長（例如日光燈汞線 435.8 / 546.1 / 611.6 nm）→ 輸出強度對波長的"
                 "光譜。可用 Peak cursor 游標對準峰、讀出像素與波長，並一鍵寫入校準表。"
                 "\nWebcam-spectrometer image → calibrated spectrum.",
                 "spectrometer")

        fbox = gui.widgetBox(self.controlArea, "Image")
        gui.button(fbox, self, "Choose spectrum photo…", callback=self._choose)
        self.file_label = gui.label(fbox, self, "(no file)")
        gui.comboBox(
            fbox, self, "rotate_idx", label="Rotate 旋轉:",
            items=[label for label, _ in ROTATIONS], callback=self._recompute,
            orientation="horizontal",
            tooltip="手持分光器拍的照片色散方向常是垂直的；轉 90°/270° 後就能沿用"
                    "水平 ROI。旋轉是無損的（只換軸，不做內插）。")

        rbox = gui.widgetBox(self.controlArea, "Strip (ROI) & channel")
        gui.comboBox(rbox, self, "channel_idx", label="Channel:",
                     items=CHANNELS, callback=self._recompute,
                     orientation="horizontal")
        gui.hSlider(rbox, self, "row_center_pct", minValue=0, maxValue=100,
                    label="Strip centre (% h):", callback=self._recompute)
        gui.hSlider(rbox, self, "row_frac_pct", minValue=1, maxValue=100,
                    label="Strip height (% h):", callback=self._recompute)
        gui.checkBox(rbox, self, "flip", "Flip (reverse the axis)",
                     callback=self._recompute)

        pbox = gui.widgetBox(self.controlArea, "Peak cursor 峰值游標")
        self.cursor_ctrl = gui.hSlider(
            pbox, self, "cursor_px", minValue=0, maxValue=1,
            label="Cursor (px):", callback=self._cursor_moved)
        self.cursor_label = gui.label(self, pbox, "")
        jump = gui.hBox(pbox)
        gui.button(jump, self, "◀ 前一峰", callback=lambda: self._jump_peak(-1))
        gui.button(jump, self, "下一峰 ▶", callback=lambda: self._jump_peak(+1))

        dbox = gui.widgetBox(pbox, "Peak detection 自動找峰")
        gui.spin(dbox, self, "peak_prom_pct", 1, 60, label="Min prominence (%):",
                 controlWidth=60, callback=self._recompute)
        gui.spin(dbox, self, "peak_dist", 1, 500, label="Min distance (px):",
                 controlWidth=60, callback=self._recompute)
        gui.spin(dbox, self, "peak_smooth", 1, 51, step=2,
                 label="Smoothing (px):", controlWidth=60,
                 callback=self._recompute)
        self.peak_label = gui.label(self, dbox, "")

        abox = gui.hBox(pbox)
        gui.comboBox(abox, self, "assign_nm", label="λ =",
                     items=[f"{v}" for v in FLUORESCENT_LINES],
                     editable=True, sendSelectedValue=True,
                     orientation="horizontal")
        gui.button(abox, self, "寫進校準表 ↓", callback=self._assign_cursor,
                   tooltip="把游標目前的次像素位置與左邊的波長，追加到下面的 "
                           "pixel=nm pairs")

        cbox = gui.widgetBox(self.controlArea, "Wavelength calibration")
        gui.comboBox(cbox, self, "model_idx", label="Fit:",
                     items=[n for n, _ in CAL_MODELS], callback=self._recompute,
                     orientation="horizontal")
        self.cal_edit = gui.lineEdit(cbox, self, "cal_text",
                                     label="pixel=nm pairs:",
                                     callback=self._recompute)
        gui.checkBox(cbox, self, "zoom_cal",
                     "Zoom plot to the calibrated range 只顯示校準區間",
                     callback=self._recompute,
                     tooltip="多項式在校準點以外會急速發散；勾選後只是把圖縮到"
                             "校準區間附近，輸出的光譜不受影響。")
        gui.label(cbox, self, "e.g. 120=435.8, 410=546.1, 560=611.6\n"
                              "(empty → x axis stays in pixels)")
        self.info_label = gui.label(
            gui.widgetBox(self.controlArea, "Status"), self, "No image.")

        self.figure = Figure(figsize=(7, 6))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.mainArea.layout().addWidget(self.canvas)
        self.ax_img = self.figure.add_subplot(211)
        self.ax_spec = self.figure.add_subplot(212)

        if self.image_path and os.path.exists(self.image_path):
            self._load(self.image_path)
        else:
            self.Information.no_image()

    # ------------------------------------------------------------- helpers
    @property
    def rotate_deg(self) -> int:
        """Clockwise rotation currently selected, in degrees."""
        return ROTATIONS[self.rotate_idx][1]

    @property
    def _slider(self):
        """The QSlider inside the cursor control (gui.hSlider may wrap it)."""
        ctrl = getattr(self, "cursor_ctrl", None)
        if isinstance(ctrl, QSlider):
            return ctrl
        return ctrl.findChild(QSlider) if ctrl is not None else None

    def _choose(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Choose spectrum photo", self.image_path or "",
            "Images (*.png *.jpg *.jpeg *.tif *.tiff *.bmp)")
        if path:
            self._load(path)

    def _load(self, path):
        self.Error.load_failed.clear()
        self.Information.no_image.clear()
        try:
            self._rgb = load_rgb(path)
        except Exception as exc:                       # noqa: BLE001
            self.Error.load_failed(str(exc))
            self._rgb = None
            return
        self.image_path = path
        self.file_label.setText(os.path.basename(path))
        self._recompute()

    def _parse_calibration(self):
        pts = []
        for tok in self.cal_text.replace(";", ",").replace("\n", ",").split(","):
            tok = tok.strip()
            if not tok:
                continue
            if "=" not in tok:
                raise ValueError(f"'{tok}' is not pixel=nm")
            p, w = tok.split("=", 1)
            pts.append((float(p), float(w)))
        return pts

    # ------------------------------------------------------- peak / cursor
    def _find_peaks(self):
        """Detect peaks in the current ROI profile (pixel coordinates)."""
        y = self._profile
        if y is None or y.size < 3:
            self._peaks = []
            self.peak_label.setText("")
            return
        rng = float(np.ptp(y)) or 1.0
        self._peaks = find_profile_peaks(
            y, min_prominence=self.peak_prom_pct / 100.0 * rng,
            min_distance=int(self.peak_dist), smooth=int(self.peak_smooth))
        if not self._peaks:
            self.peak_label.setText("找不到峰（試著調低 prominence）")
            return
        shown = ", ".join(f"{p['position']:.1f}" for p in self._peaks[:8])
        more = " …" if len(self._peaks) > 8 else ""
        self.peak_label.setText(f"{len(self._peaks)} 個峰 (px): {shown}{more}")

    def _cursor_position(self):
        """(sub-pixel position, snapped?) — snaps to a peak within ±4 px."""
        if self._profile is None or self._profile.size == 0:
            return 0.0, False
        i = int(np.clip(self.cursor_px, 0, self._profile.size - 1))
        near = [p for p in self._peaks if abs(p["index"] - i) <= 4]
        if near:
            best = min(near, key=lambda p: abs(p["index"] - i))
            return best["position"], True
        return float(i), False

    def _jump_peak(self, step):
        if not self._peaks:
            return
        pos = [p["index"] for p in self._peaks]
        if step > 0:
            later = [q for q in pos if q > self.cursor_px]
            self.cursor_px = later[0] if later else pos[0]
        else:
            earlier = [q for q in pos if q < self.cursor_px]
            self.cursor_px = earlier[-1] if earlier else pos[-1]
        sl = self._slider
        if sl is not None:
            sl.setValue(int(self.cursor_px))
        self._cursor_moved()

    def _cursor_moved(self):
        """Cheap update: readout text + the cursor line, no full recompute."""
        if self._profile is None or self._profile.size == 0:
            self.cursor_label.setText("")
            return
        sub, snapped = self._cursor_position()
        i = int(np.clip(round(sub), 0, self._profile.size - 1))
        nm = (pixel_to_wavelength(sub, coeffs=self._coeffs)
              if self._coeffs is not None else None)
        txt = f"pixel {sub:8.2f}{'  ← peak' if snapped else '        '}   " \
              f"I {self._profile[i]:7.2f}   "
        if nm is None:
            txt += "λ  (未校準)"
        else:
            txt += f"λ {nm:7.2f} nm"
            lo, hi = self._cal_span
            if lo is not None and not (lo <= sub <= hi):
                txt += "  ⚠ 外插 extrapolated"
        self.cursor_label.setText(txt)
        self._draw_cursor(sub)

    def _assign_cursor(self):
        """Append 'cursor=λ' to the calibration line and re-fit."""
        if self._profile is None:
            return
        try:
            nm = float(str(self.assign_nm).strip())
        except (TypeError, ValueError):
            self.Error.bad_calibration(f"'{self.assign_nm}' is not a number")
            return
        sub, _ = self._cursor_position()
        entry = f"{sub:.2f}={nm:g}"
        current = (self.cal_text or "").strip().rstrip(",")
        self.cal_text = f"{current}, {entry}" if current else entry
        if getattr(self, "cal_edit", None) is not None:
            self.cal_edit.setText(self.cal_text)
        self._recompute()

    def _draw_cursor(self, pos_px):
        """Draw the cursor (and the detected peaks) on top of the spectrum."""
        for artist in self._cursor_artists:
            try:
                artist.remove()
            except (ValueError, NotImplementedError):
                pass
        self._cursor_artists = []
        if self._profile is None:
            return

        def to_x(p):
            if self._spec_x is None:
                return float(p)
            return float(np.interp(p, np.arange(self._profile.size),
                                   self._spec_x))

        for p in self._peaks:
            self._cursor_artists.append(self.ax_spec.axvline(
                to_x(p["position"]), color="#999999", ls=":", lw=0.7, zorder=1))
        self._cursor_artists.append(self.ax_spec.axvline(
            to_x(pos_px), color="#c44e52", ls="-", lw=1.6, zorder=4))
        self.canvas.draw_idle()

    # ---------------------------------------------------------- main logic
    def _recompute(self):
        self.Error.bad_calibration.clear()
        self.Warning.exact_fit.clear()
        self.Warning.folded.clear()
        self.ax_img.clear()
        self.ax_spec.clear()
        self._cursor_artists = []
        if self._rgb is None:
            self._profile, self._peaks, self._coeffs, self._spec_x = \
                None, [], None, None
            self._cal_span = (None, None)
            self.info_label.setText("No image.")
            self.canvas.draw_idle()
            self.Outputs.spectrum.send(None)
            return

        # The ROI is taken from the rotated image, so the preview must show it.
        view = rotate_rgb(self._rgb, self.rotate_deg)
        h = view.shape[0]
        row_center = self.row_center_pct / 100.0 * h
        row_frac = self.row_frac_pct / 100.0

        try:
            calib = self._parse_calibration()
        except ValueError as exc:
            self.Error.bad_calibration(str(exc))
            self.Outputs.spectrum.send(None)
            self.canvas.draw_idle()
            return

        try:
            spec = image_to_spectrum(
                self._rgb, channel=CHANNELS[self.channel_idx],
                row_center=row_center, row_frac=row_frac,
                calibration=calib or None, model=CAL_MODELS[self.model_idx][0],
                name=os.path.splitext(os.path.basename(self.image_path))[0]
                or "spectrum", flip=self.flip, rotate=self.rotate_deg)
        except ValueError as exc:
            self.Error.bad_calibration(str(exc))
            self.Outputs.spectrum.send(None)
            self.canvas.draw_idle()
            return

        # Profile in *calibration* coordinates: rotated + flipped, never the
        # ascending-wavelength reversal image_to_spectrum may apply at the end.
        profile = extract_profile(self._rgb, channel=CHANNELS[self.channel_idx],
                                  row_center=row_center, row_frac=row_frac,
                                  rotate=self.rotate_deg)
        self._profile = profile[::-1] if self.flip else profile
        cal = spec.get("calibration")
        self._coeffs = cal["coeffs"] if cal else None
        if cal:
            ps = [p for p, _ in cal["points"]]
            self._cal_span = (min(ps), max(ps))
            if len(ps) == cal["degree"] + 1:
                self.Warning.exact_fit(n=len(ps), d=cal["degree"])
        else:
            self._cal_span = (None, None)
        self._spec_x = (np.polyval(self._coeffs,
                                   np.arange(self._profile.size, dtype=float))
                        if self._coeffs is not None else None)

        clipped = spec.get("clipped")
        if clipped:
            self.Warning.folded(n=self._profile.size - (clipped[1] - clipped[0]))

        n = int(self._profile.size)
        sl = self._slider
        if sl is not None:
            sl.setMaximum(max(0, n - 1))
        self.cursor_px = int(np.clip(self.cursor_px, 0, max(0, n - 1)))
        self._find_peaks()

        # top: the (rotated) image with the ROI strip drawn
        self.ax_img.imshow(np.clip(view / 255.0, 0, 1), aspect="auto")
        half = max(0.5, row_frac * h / 2.0)
        self.ax_img.add_patch(Rectangle(
            (0, row_center - half), view.shape[1], 2 * half, fill=False,
            edgecolor="#ffd54f", lw=1.5))
        self.ax_img.set_xticks([])
        self.ax_img.set_yticks([])
        title = "ROI strip"
        if self.rotate_deg:
            title += f"  (rotated {self.rotate_deg}° CW)"
        self.ax_img.set_title(title, fontsize=9)

        # bottom: the extracted spectrum
        self.ax_spec.plot(spec["x"], spec["y"], color="#4c72b0", lw=1.0)
        self.ax_spec.set_xlabel(spec["x_label"])
        self.ax_spec.set_ylabel(f"intensity ({CHANNELS[self.channel_idx]})")
        if cal:
            for p, w in cal["points"]:
                self.ax_spec.axvline(w, color="#c44e52", ls="--", lw=0.8)
            if self.zoom_cal and self._spec_x is not None:
                lo_px, hi_px = self._cal_span
                pad = 0.2 * max(hi_px - lo_px, 1.0)
                edges = [np.polyval(self._coeffs, lo_px - pad),
                         np.polyval(self._coeffs, hi_px + pad)]
                self.ax_spec.set_xlim(min(edges), max(edges))
        self.figure.tight_layout()
        self._cursor_moved()          # draws peaks + cursor, then draw_idle()

        out = table_from_spectra([spec])
        out.name = "spectrometer"
        msg = f"{spec['x'].size} points"
        if self.rotate_deg:
            msg += f"  ·  rotated {self.rotate_deg}°"
        if cal:
            msg += (f"  ·  {CAL_MODELS[self.model_idx][0]} calibration, "
                    f"R²={cal['r2']:.4f}, {len(cal['points'])} lines")
        else:
            msg += "  ·  no calibration (x = pixel)"
        if self._peaks:
            msg += f"  ·  {len(self._peaks)} peaks"
        self.info_label.setText(msg)
        self.Outputs.spectrum.send(out)

    def send_report(self):
        self.report_items("Spectrometer", [
            ("Image", os.path.basename(self.image_path) if self.image_path else "—"),
            ("Rotation", ROTATIONS[self.rotate_idx][0]),
            ("Channel", CHANNELS[self.channel_idx]),
            ("Calibration", self.cal_text or "none"),
            ("Peaks (px)", ", ".join(f"{p['position']:.1f}"
                                     for p in self._peaks) or "—")])
        if self._rgb is not None:
            self.report_plot(self.figure)


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWSpectrometer).run()
