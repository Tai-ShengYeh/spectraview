"""Spectrometer — turn a diffraction-grating spectrum photo into a spectrum.

A Theremino-style spectrometer working on a captured image: take a horizontal
strip, read an intensity profile along the columns, and calibrate pixel ->
wavelength from known reference lines. Outputs an Orange spectrum Table.
"""
import os

import numpy as np

from AnyQt.QtWidgets import QFileDialog

import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from Orange.data import Table
from Orange.widgets import gui, settings
from Orange.widgets.widget import Msg, Output, OWWidget

from ..spectrometer import CAL_MODELS, CHANNELS, image_to_spectrum, load_rgb
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
                "calibration", "分光", "光柵"]

    class Outputs:
        spectrum = Output("Spectra", Table)

    class Error(OWWidget.Error):
        load_failed = Msg("{}")
        bad_calibration = Msg("Calibration: {}")

    class Information(OWWidget.Information):
        no_image = Msg("Choose a spectrum photo to read.")

    image_path: str = settings.Setting("")
    channel_idx: int = settings.Setting(0)
    row_center_pct: int = settings.Setting(50)         # % of height
    row_frac_pct: int = settings.Setting(15)           # band height %
    flip: bool = settings.Setting(False)
    model_idx: int = settings.Setting(0)
    cal_text: str = settings.Setting("")               # "pixel=nm" per line/comma
    want_main_area = True

    def __init__(self):
        super().__init__()
        self._rgb = None

        add_help(self,
                 "選一張『相機＋繞射光柵』拍到的光譜照片 → 取一條水平帶、沿水平方向"
                 "讀強度 → 用已知譜線把像素校準成波長（例如日光燈汞線 435.8 / 546.1 / "
                 "611.6 nm）→ 輸出強度對波長的光譜。校準格式：每行或逗號分隔的 "
                 "`像素=波長`。\nWebcam-spectrometer image → calibrated spectrum.",
                 "spectrometer")

        fbox = gui.widgetBox(self.controlArea, "Image")
        gui.button(fbox, self, "Choose spectrum photo…", callback=self._choose)
        self.file_label = gui.label(fbox, self, "(no file)")

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

        cbox = gui.widgetBox(self.controlArea, "Wavelength calibration")
        gui.comboBox(cbox, self, "model_idx", label="Fit:",
                     items=[n for n, _ in CAL_MODELS], callback=self._recompute,
                     orientation="horizontal")
        gui.lineEdit(cbox, self, "cal_text",
                     label="pixel=nm pairs:", callback=self._recompute)
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

    def _recompute(self):
        self.Error.bad_calibration.clear()
        self.ax_img.clear()
        self.ax_spec.clear()
        if self._rgb is None:
            self.info_label.setText("No image.")
            self.canvas.draw_idle()
            self.Outputs.spectrum.send(None)
            return

        h = self._rgb.shape[0]
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
                or "spectrum", flip=self.flip)
        except ValueError as exc:
            self.Error.bad_calibration(str(exc))
            self.Outputs.spectrum.send(None)
            self.canvas.draw_idle()
            return

        # top: the image with the ROI strip drawn
        self.ax_img.imshow(np.clip(self._rgb / 255.0, 0, 1), aspect="auto")
        half = max(0.5, row_frac * h / 2.0)
        self.ax_img.add_patch(Rectangle(
            (0, row_center - half), self._rgb.shape[1], 2 * half, fill=False,
            edgecolor="#ffd54f", lw=1.5))
        self.ax_img.set_xticks([]); self.ax_img.set_yticks([])
        self.ax_img.set_title("ROI strip", fontsize=9)

        # bottom: the extracted spectrum
        self.ax_spec.plot(spec["x"], spec["y"], color="#4c72b0", lw=1.0)
        self.ax_spec.set_xlabel(spec["x_label"])
        self.ax_spec.set_ylabel(f"intensity ({CHANNELS[self.channel_idx]})")
        cal = spec.get("calibration")
        if cal:
            for p, w in cal["points"]:
                self.ax_spec.axvline(w, color="#c44e52", ls="--", lw=0.8)
        self.figure.tight_layout()
        self.canvas.draw_idle()

        out = table_from_spectra([spec])
        out.name = "spectrometer"
        msg = f"{spec['x'].size} points"
        if cal:
            msg += (f"  ·  {CAL_MODELS[self.model_idx][0]} calibration, "
                    f"R²={cal['r2']:.4f}, {len(cal['points'])} lines")
        else:
            msg += "  ·  no calibration (x = pixel)"
        self.info_label.setText(msg)
        self.Outputs.spectrum.send(out)

    def send_report(self):
        self.report_items("Spectrometer", [
            ("Image", os.path.basename(self.image_path) if self.image_path else "—"),
            ("Channel", CHANNELS[self.channel_idx]),
            ("Calibration", self.cal_text or "none")])
        if self._rgb is not None:
            self.report_plot(self.figure)


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWSpectrometer).run()
