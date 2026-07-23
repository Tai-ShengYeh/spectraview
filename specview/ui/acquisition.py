"""Live acquisition dialog: Ocean Optics (seabreeze) + Vernier Go Direct.

Opened from the main window's  儀器 (Device) ▸ 光譜儀擷取…  menu. Non-modal, so
the main plot stays usable while acquiring.

Ocean Optics tab
    connect → live preview (worker thread; the GUI never blocks on the
    integration time) → optional dark / reference frames → raw counts,
    dark-corrected, %T or absorbance → 「加入光譜清單」 hands the spectrum to
    the main window.

Go Direct tab
    Best-effort godirect connection (device info + channel list + live scalar
    values — works for ordinary Go Direct sensors). For the SpectroVis Plus the
    public godirect protocol carries no spectrum, so the tab also offers a
    watch folder: export CSV from Vernier "Spectral Analysis" into that folder
    and each new file is imported automatically.
"""
from __future__ import annotations

import os
import threading
import time
import traceback

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets
import pyqtgraph as pg

from ..devices import DeviceError
from ..devices import godirect_dev, oceanoptics
from ..formats.binary_io import MissingDependency
from ..spectrum import Spectrum


# ============================================================ worker thread
class _AcqWorker(QtCore.QThread):
    """Continuously reads spectra from an open device.

    The read blocks for the integration time, so it lives in its own thread;
    frames arrive on the GUI thread through the ``frame`` signal.
    """

    frame = QtCore.Signal(object)          # np.ndarray (intensities)
    failed = QtCore.Signal(str)

    def __init__(self, device, parent=None):
        super().__init__(parent)
        self._dev = device
        self._stop = threading.Event()
        self._lock = threading.Lock()
        self._pending_ms: float | None = None
        self._averages = 1

    # thread-safe parameter setters (called from the GUI thread) ------------
    def request_integration(self, ms: float) -> None:
        with self._lock:
            self._pending_ms = ms

    def set_averages(self, n: int) -> None:
        with self._lock:
            self._averages = max(1, int(n))

    def stop(self) -> None:
        self._stop.set()

    # -----------------------------------------------------------------------
    def run(self) -> None:
        while not self._stop.is_set():
            try:
                with self._lock:
                    ms, self._pending_ms = self._pending_ms, None
                    avg = self._averages
                if ms is not None:
                    self._dev.set_integration_time_ms(ms)
                acc = None
                for _ in range(avg):
                    if self._stop.is_set():
                        return
                    y = self._dev.read_intensities()
                    acc = y if acc is None else acc + y
                self.frame.emit(acc / avg)
            except Exception as exc:  # noqa: BLE001
                if not self._stop.is_set():
                    self.failed.emit(str(exc))
                return


# ============================================================ Ocean Optics tab
class _OceanTab(QtWidgets.QWidget):
    def __init__(self, add_spectrum, status, parent=None):
        super().__init__(parent)
        self._add_spectrum = add_spectrum     # callable(Spectrum, label)
        self._status = status                 # callable(str)
        self._dev: oceanoptics.OceanOpticsDevice | None = None
        self._worker: _AcqWorker | None = None
        self._x = None                        # wavelength axis
        self._y = None                        # latest raw frame
        self._dark = None
        self._ref = None
        self._n_captured = 0
        self._build()

    # ------------------------------------------------------------------ UI
    def _build(self) -> None:
        lay = QtWidgets.QVBoxLayout(self)

        # --- connection row
        row = QtWidgets.QHBoxLayout()
        self.cmb_dev = QtWidgets.QComboBox()
        self.cmb_dev.setMinimumWidth(240)
        self.btn_refresh = QtWidgets.QPushButton("重新掃描")
        self.btn_conn = QtWidgets.QPushButton("連線")
        row.addWidget(QtWidgets.QLabel("裝置:"))
        row.addWidget(self.cmb_dev, 1)
        row.addWidget(self.btn_refresh)
        row.addWidget(self.btn_conn)
        lay.addLayout(row)

        # --- parameters row
        prm = QtWidgets.QHBoxLayout()
        self.spn_it = QtWidgets.QDoubleSpinBox()
        self.spn_it.setRange(0.001, 65_000.0)
        self.spn_it.setValue(100.0)
        self.spn_it.setDecimals(3)
        self.spn_it.setSuffix(" ms")
        self.spn_avg = QtWidgets.QSpinBox()
        self.spn_avg.setRange(1, 1000)
        self.spn_avg.setValue(1)
        self.cmb_mode = QtWidgets.QComboBox()
        self.cmb_mode.addItems(["Raw counts", "扣暗背景 (S − D)",
                                "穿透率 %T", "吸收度 A"])
        prm.addWidget(QtWidgets.QLabel("積分時間:"))
        prm.addWidget(self.spn_it)
        prm.addWidget(QtWidgets.QLabel("平均次數:"))
        prm.addWidget(self.spn_avg)
        prm.addWidget(QtWidgets.QLabel("顯示模式:"))
        prm.addWidget(self.cmb_mode, 1)
        lay.addLayout(prm)

        # --- live preview
        self.plot = pg.PlotWidget()
        self.plot.setLabel("bottom", "Wavelength (nm)")
        self.plot.setLabel("left", "Counts")
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        self.curve = self.plot.plot(pen=pg.mkPen("#1f77b4", width=1.4))
        lay.addWidget(self.plot, 1)

        # --- acquisition buttons
        btns = QtWidgets.QHBoxLayout()
        self.btn_live = QtWidgets.QPushButton("▶ 開始即時預覽")
        self.btn_live.setCheckable(True)
        self.btn_dark = QtWidgets.QPushButton("存為暗背景 D")
        self.btn_ref = QtWidgets.QPushButton("存為參考 R")
        self.btn_capture = QtWidgets.QPushButton("➕ 加入光譜清單")
        for b in (self.btn_live, self.btn_dark, self.btn_ref, self.btn_capture):
            btns.addWidget(b)
        lay.addLayout(btns)

        self.lbl_info = QtWidgets.QLabel(
            "未連線。 提示: 一次只能有一個程式使用光譜儀 (請先關閉 OceanView)。")
        self.lbl_info.setWordWrap(True)
        self.lbl_info.setStyleSheet("color: palette(mid);")
        lay.addWidget(self.lbl_info)

        # wiring
        self.btn_refresh.clicked.connect(self.refresh_devices)
        self.btn_conn.clicked.connect(self.toggle_connect)
        self.btn_live.toggled.connect(self.toggle_live)
        self.btn_dark.clicked.connect(lambda: self._store_frame("dark"))
        self.btn_ref.clicked.connect(lambda: self._store_frame("ref"))
        self.btn_capture.clicked.connect(self.capture)
        self.spn_it.valueChanged.connect(self._push_params)
        self.spn_avg.valueChanged.connect(self._push_params)
        self.cmb_mode.currentIndexChanged.connect(self._redraw)
        self._set_acq_enabled(False)

    def _set_acq_enabled(self, on: bool) -> None:
        for w in (self.btn_live, self.btn_dark, self.btn_ref,
                  self.btn_capture, self.spn_it, self.spn_avg, self.cmb_mode):
            w.setEnabled(on)

    # ------------------------------------------------------------ connection
    def refresh_devices(self) -> None:
        self.cmb_dev.clear()
        try:
            devs = oceanoptics.list_devices()
        except MissingDependency as exc:
            QtWidgets.QMessageBox.information(self, "需要安裝 seabreeze", str(exc))
            return
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.warning(self, "掃描失敗", str(exc))
            return
        for d in devs:
            self.cmb_dev.addItem(d.label, d.ident)
        self.lbl_info.setText(f"找到 {len(devs)} 台 Ocean Optics 光譜儀。"
                              if devs else
                              "沒有找到光譜儀 — 檢查 USB 線與驅動 (seabreeze_os_setup)。")

    def toggle_connect(self) -> None:
        if self._dev is not None:
            self.disconnect_device()
            return
        serial = self.cmb_dev.currentData()
        dev = oceanoptics.OceanOpticsDevice(serial=serial)
        try:
            dev.open()
        except (MissingDependency, DeviceError) as exc:
            QtWidgets.QMessageBox.warning(self, "連線失敗", str(exc))
            return
        self._dev = dev
        self._x = dev.wavelengths()
        lo, hi = dev.integration_limits_ms
        self.spn_it.setRange(lo, hi)
        self.spn_it.setValue(float(np.clip(self.spn_it.value(), lo, hi)))
        dev.set_integration_time_ms(self.spn_it.value())
        self.btn_conn.setText("中斷連線")
        self._set_acq_enabled(True)
        self.lbl_info.setText(
            f"已連線 {dev.label} — {self._x.size} 像素, "
            f"{self._x[0]:.1f}–{self._x[-1]:.1f} nm, "
            f"積分時間 {lo:.3g}–{hi:.3g} ms")
        self._status(f"Connected: {dev.label}")

    def disconnect_device(self) -> None:
        self._stop_worker()
        if self._dev is not None:
            self._dev.close()
            self._dev = None
        self.btn_conn.setText("連線")
        self.btn_live.setChecked(False)
        self._set_acq_enabled(False)
        self.lbl_info.setText("已中斷連線。")

    # ------------------------------------------------------------ live loop
    def toggle_live(self, on: bool) -> None:
        if on:
            if self._dev is None:
                self.btn_live.setChecked(False)
                return
            self.btn_live.setText("⏸ 停止即時預覽")
            self._worker = _AcqWorker(self._dev)
            self._worker.frame.connect(self._on_frame)
            self._worker.failed.connect(self._on_failed)
            self._push_params()
            self._worker.start()
        else:
            self.btn_live.setText("▶ 開始即時預覽")
            self._stop_worker()

    def _stop_worker(self) -> None:
        if self._worker is not None:
            self._worker.stop()
            self._worker.wait(3000)
            self._worker = None

    def _push_params(self) -> None:
        if self._worker is not None:
            self._worker.request_integration(self.spn_it.value())
            self._worker.set_averages(self.spn_avg.value())
        elif self._dev is not None:
            try:
                self._dev.set_integration_time_ms(self.spn_it.value())
            except DeviceError:
                pass

    def _on_frame(self, y) -> None:
        self._y = np.asarray(y, dtype=float)
        self._redraw()

    def _on_failed(self, msg: str) -> None:
        self.btn_live.setChecked(False)
        QtWidgets.QMessageBox.warning(self, "擷取中斷", msg)

    # ------------------------------------------------------------ processing
    def _store_frame(self, kind: str) -> None:
        y = self._current_raw()
        if y is None:
            return
        if kind == "dark":
            self._dark = y.copy()
            self.lbl_info.setText("已儲存暗背景 D (蓋住入光口再按一次可更新)。")
        else:
            self._ref = y.copy()
            self.lbl_info.setText("已儲存參考 R (100% 透射/空白溶液)。")
        self._redraw()

    def _current_raw(self):
        """Latest raw frame — from live preview, or a one-shot read."""
        if self._y is not None:
            return self._y
        if self._dev is None:
            return None
        try:
            return self._dev.read_intensities()
        except DeviceError as exc:
            QtWidgets.QMessageBox.warning(self, "讀取失敗", str(exc))
            return None

    def _processed(self):
        """Apply the display mode. Returns (y, y_unit) or None."""
        y = self._y
        if y is None:
            return None
        mode = self.cmb_mode.currentIndex()
        if mode == 0:
            return y, "counts"
        if self._dark is None:
            return y, "counts"          # no dark yet -> raw
        s = y - self._dark
        if mode == 1:
            return s, "counts"
        if self._ref is None:
            return s, "counts"          # no reference yet -> dark-corrected
        r = self._ref - self._dark
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(r > 0, s / r, np.nan)
            if mode == 2:
                return t * 100.0, "%T"
            a = -np.log10(np.clip(t, 1e-10, None))
        return a, "absorbance"

    def _redraw(self) -> None:
        out = self._processed()
        if out is None:
            return
        y, unit = out
        self.curve.setData(self._x, y)
        self.plot.setLabel("left", {"counts": "Counts", "%T": "Transmittance (%)",
                                    "absorbance": "Absorbance"}[unit])

    # ------------------------------------------------------------ capture
    def capture(self) -> None:
        if self._y is None:                      # not previewing: take one shot
            y = self._current_raw()
            if y is None:
                return
            self._y = y
        out = self._processed()
        if out is None:
            return
        y, unit = out
        self._n_captured += 1
        name = f"{self._dev.label} #{self._n_captured}"
        spec = Spectrum(self._x.copy(), y, name=name, x_unit="nm", y_unit=unit,
                        meta={"source": "live:oceanoptics",
                              "integration_ms": self.spn_it.value(),
                              "averages": self.spn_avg.value(),
                              "dark_corrected": self._dark is not None,
                              "referenced": self._ref is not None})
        self._add_spectrum(spec, f"Captured {name}")

    # ------------------------------------------------------------ teardown
    def shutdown(self) -> None:
        self.disconnect_device()


# ============================================================ Go Direct tab
class _GoDirectTab(QtWidgets.QWidget):
    def __init__(self, load_paths, status, parent=None):
        super().__init__(parent)
        self._load_paths = load_paths          # callable(list[str])
        self._status = status
        self._dev: godirect_dev.GoDirectDevice | None = None
        self._timer = QtCore.QTimer(self)
        self._timer.setInterval(400)
        self._timer.timeout.connect(self._poll)
        self._watcher: QtCore.QFileSystemWatcher | None = None
        self._seen: set[str] = set()
        self._build()

    def _build(self) -> None:
        lay = QtWidgets.QVBoxLayout(self)

        # --- connection group
        g1 = QtWidgets.QGroupBox("Go Direct 連線 (godirect)")
        v1 = QtWidgets.QVBoxLayout(g1)
        row = QtWidgets.QHBoxLayout()
        self.rad_usb = QtWidgets.QRadioButton("USB")
        self.rad_ble = QtWidgets.QRadioButton("藍牙 BLE")
        self.rad_usb.setChecked(True)
        self.btn_conn = QtWidgets.QPushButton("掃描並連線")
        row.addWidget(self.rad_usb)
        row.addWidget(self.rad_ble)
        row.addStretch(1)
        row.addWidget(self.btn_conn)
        v1.addLayout(row)
        self.txt_sensors = QtWidgets.QPlainTextEdit()
        self.txt_sensors.setReadOnly(True)
        self.txt_sensors.setMaximumHeight(140)
        self.txt_sensors.setPlaceholderText("連線後顯示裝置與感測通道…")
        v1.addWidget(self.txt_sensors)
        lay.addWidget(g1)

        # --- limitation note
        note = QtWidgets.QLabel("⚠ " + godirect_dev.SPECTRO_NOTE)
        note.setWordWrap(True)
        note.setStyleSheet("color: #b5651d;")
        lay.addWidget(note)

        # --- watch folder group
        g2 = QtWidgets.QGroupBox("監看 Spectral Analysis 匯出資料夾 (自動載入新 CSV)")
        v2 = QtWidgets.QVBoxLayout(g2)
        row2 = QtWidgets.QHBoxLayout()
        self.edit_dir = QtWidgets.QLineEdit()
        self.edit_dir.setPlaceholderText("選擇 Spectral Analysis 匯出 CSV 的資料夾…")
        self.btn_dir = QtWidgets.QPushButton("瀏覽…")
        self.chk_watch = QtWidgets.QCheckBox("啟用監看")
        row2.addWidget(self.edit_dir, 1)
        row2.addWidget(self.btn_dir)
        row2.addWidget(self.chk_watch)
        v2.addLayout(row2)
        self.lbl_watch = QtWidgets.QLabel(
            "步驟: Spectral Analysis 取光譜 → File ▸ Export ▸ CSV 存到上面資料夾 "
            "→ SpectraView 立即自動載入。")
        self.lbl_watch.setWordWrap(True)
        self.lbl_watch.setStyleSheet("color: palette(mid);")
        v2.addWidget(self.lbl_watch)
        lay.addWidget(g2)
        lay.addStretch(1)

        self.btn_conn.clicked.connect(self.toggle_connect)
        self.btn_dir.clicked.connect(self.choose_dir)
        self.chk_watch.toggled.connect(self.toggle_watch)

    # ------------------------------------------------------------ godirect
    def toggle_connect(self) -> None:
        if self._dev is not None:
            self.disconnect_device()
            return
        dev = godirect_dev.GoDirectDevice(use_usb=self.rad_usb.isChecked(),
                                          use_ble=self.rad_ble.isChecked())
        try:
            dev.open()
        except (MissingDependency, DeviceError) as exc:
            QtWidgets.QMessageBox.warning(self, "Go Direct 連線失敗", str(exc))
            return
        self._dev = dev
        self.btn_conn.setText("中斷連線")
        lines = [f"已連線: {dev.label}"]
        try:
            names = dev.sensor_names()
            lines += [f"  • {n}" for n in names] or ["  (無公開感測通道)"]
        except DeviceError as exc:
            lines.append(f"  通道讀取失敗: {exc}")
        if dev.is_spectrometer():
            lines.append("")
            lines.append("此裝置是光譜儀 — godirect 協定不含光譜資料,"
                         "請改用下方的匯出資料夾監看。")
        else:
            try:
                dev.start(period_ms=500)
                self._timer.start()
                lines.append("開始讀取即時數值…")
            except DeviceError as exc:
                lines.append(f"啟動讀取失敗: {exc}")
        self.txt_sensors.setPlainText("\n".join(lines))
        self._status(f"Go Direct connected: {dev.label}")

    def _poll(self) -> None:
        if self._dev is None:
            return
        try:
            vals = self._dev.read_values()
        except DeviceError:
            return
        if vals:
            head = self.txt_sensors.toPlainText().splitlines()
            keep = [l for l in head if not l.startswith("→ ")]
            keep += [f"→ {k}: {v:.4g}" for k, v in vals.items()]
            self.txt_sensors.setPlainText("\n".join(keep))

    def disconnect_device(self) -> None:
        self._timer.stop()
        if self._dev is not None:
            self._dev.close()
            self._dev = None
        self.btn_conn.setText("掃描並連線")
        self.txt_sensors.setPlainText("已中斷連線。")

    # ------------------------------------------------------------ watcher
    def choose_dir(self) -> None:
        d = QtWidgets.QFileDialog.getExistingDirectory(self, "選擇匯出資料夾")
        if d:
            self.edit_dir.setText(d)

    def toggle_watch(self, on: bool) -> None:
        if on:
            folder = self.edit_dir.text().strip()
            if not os.path.isdir(folder):
                QtWidgets.QMessageBox.information(self, "資料夾無效",
                                                  "請先選擇存在的資料夾。")
                self.chk_watch.setChecked(False)
                return
            self._seen = {os.path.join(folder, f) for f in os.listdir(folder)}
            self._watcher = QtCore.QFileSystemWatcher([folder], self)
            self._watcher.directoryChanged.connect(self._on_dir_changed)
            self.lbl_watch.setText(f"監看中: {folder}")
        else:
            if self._watcher is not None:
                self._watcher.deleteLater()
                self._watcher = None
            self.lbl_watch.setText("已停止監看。")

    def _on_dir_changed(self, folder: str) -> None:
        try:
            now = {os.path.join(folder, f) for f in os.listdir(folder)}
        except OSError:
            return
        fresh = sorted(p for p in now - self._seen
                       if godirect_dev.is_spectral_analysis_export(p))
        self._seen = now
        if not fresh:
            return
        # Give the writing program a moment to finish the file.
        QtCore.QTimer.singleShot(
            600, lambda paths=fresh: self._import_paths(paths))

    def _import_paths(self, paths: list[str]) -> None:
        ready = [p for p in paths if os.path.isfile(p) and os.path.getsize(p) > 0]
        if ready:
            try:
                self._load_paths(ready)
                self._status(f"自動載入 {len(ready)} 個檔案")
            except Exception:  # noqa: BLE001
                traceback.print_exc()

    def shutdown(self) -> None:
        self.disconnect_device()
        self.chk_watch.setChecked(False)


# ============================================================ InnoSpectra tab
class _IscScanWorker(QtCore.QThread):
    """One InnoSpectra scan (takes seconds) off the GUI thread."""

    done = QtCore.Signal(list)      # list[Spectrum]
    failed = QtCore.Signal(str)

    def __init__(self, device, ref, kinds, name, parent=None):
        super().__init__(parent)
        self._args = (device, ref, kinds, name)

    def run(self) -> None:
        device, ref, kinds, name = self._args
        try:
            self.done.emit(device.scan(ref=ref, kinds=kinds, name=name))
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(str(exc))


class _InnoSpectraTab(QtWidgets.QWidget):
    """InnoSpectra NIRScan via the official Python SDK (32-bit bridge)."""

    def __init__(self, add_spectrum, status, parent=None):
        super().__init__(parent)
        self._add_spectrum = add_spectrum
        self._status = status
        self._dev = None
        self._worker = None
        self._n = 0
        self._settings = QtCore.QSettings("SpectraView", "acquisition")
        self._build()

    def _build(self) -> None:
        from ..devices import innospectra as isc

        lay = QtWidgets.QVBoxLayout(self)

        g0 = QtWidgets.QGroupBox("SDK 設定 (一次設定,自動記住)")
        f0 = QtWidgets.QFormLayout(g0)
        self.edit_sdk = QtWidgets.QLineEdit(
            self._settings.value("isc/sdk_dir", ""))
        self.edit_sdk.setPlaceholderText("ISC_NIRScan_PyQt-master 資料夾 (含 SDK/iscpy.pyd)")
        b1 = QtWidgets.QPushButton("瀏覽…")
        b1.clicked.connect(lambda: self._pick_dir(self.edit_sdk))
        h1 = QtWidgets.QHBoxLayout(); h1.addWidget(self.edit_sdk, 1); h1.addWidget(b1)
        f0.addRow("SDK 資料夾:", h1)
        self.edit_py = QtWidgets.QLineEdit(
            self._settings.value("isc/python32", isc._default_python32() or ""))
        self.edit_py.setPlaceholderText(r"32 位元 Python 3.11 (例 %LOCALAPPDATA%\Programs\Python\Python311-32\python.exe)")
        b2 = QtWidgets.QPushButton("瀏覽…")
        b2.clicked.connect(lambda: self._pick_file(self.edit_py))
        h2 = QtWidgets.QHBoxLayout(); h2.addWidget(self.edit_py, 1); h2.addWidget(b2)
        f0.addRow("32-bit Python:", h2)
        lay.addWidget(g0)

        row = QtWidgets.QHBoxLayout()
        self.btn_conn = QtWidgets.QPushButton("連線")
        self.lbl_info = QtWidgets.QLabel("未連線。掃描前請關閉 ISC NIRScan GUI。")
        self.lbl_info.setWordWrap(True)
        row.addWidget(self.btn_conn)
        row.addWidget(self.lbl_info, 1)
        lay.addLayout(row)

        prm = QtWidgets.QHBoxLayout()
        self.spn_rep = QtWidgets.QSpinBox(); self.spn_rep.setRange(1, 50); self.spn_rep.setValue(6)
        self.cmb_lamp = QtWidgets.QComboBox()
        self.cmb_lamp.addItem("燈:掃描時自動", "auto")
        self.cmb_lamp.addItem("燈:常亮", "on")
        self.cmb_lamp.addItem("燈:關", "off")
        self.cmb_ref = QtWidgets.QComboBox()
        self.cmb_ref.addItem("參考:前次儲存", "previous")
        self.cmb_ref.addItem("參考:出廠內建", "factory")
        self.cmb_ref.addItem("參考:本次掃描當新參考", "new")
        prm.addWidget(QtWidgets.QLabel("平均次數:")); prm.addWidget(self.spn_rep)
        prm.addWidget(self.cmb_lamp); prm.addWidget(self.cmb_ref, 1)
        lay.addLayout(prm)

        kinds = QtWidgets.QHBoxLayout()
        self.chk_abs = QtWidgets.QCheckBox("吸收度"); self.chk_abs.setChecked(True)
        self.chk_int = QtWidgets.QCheckBox("強度")
        self.chk_refl = QtWidgets.QCheckBox("反射率")
        self.chk_cont = QtWidgets.QCheckBox("連續掃描")
        kinds.addWidget(QtWidgets.QLabel("輸出:"))
        for w in (self.chk_abs, self.chk_int, self.chk_refl):
            kinds.addWidget(w)
        kinds.addStretch(1)
        kinds.addWidget(self.chk_cont)
        lay.addLayout(kinds)

        self.plot = pg.PlotWidget()
        self.plot.setLabel("bottom", "Wavelength (nm)")
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        self.curve = self.plot.plot(pen=pg.mkPen("#2ca02c", width=1.4))
        lay.addWidget(self.plot, 1)

        btns = QtWidgets.QHBoxLayout()
        self.btn_scan = QtWidgets.QPushButton("📷 掃描並加入清單")
        self.btn_saveref = QtWidgets.QPushButton("掃描並存為新參考")
        btns.addWidget(self.btn_scan); btns.addWidget(self.btn_saveref)
        lay.addLayout(btns)

        self.btn_conn.clicked.connect(self.toggle_connect)
        self.btn_scan.clicked.connect(self.do_scan)
        self.btn_saveref.clicked.connect(self.do_save_ref)
        self._set_enabled(False)

    # ------------------------------------------------------------- helpers
    def _pick_dir(self, edit) -> None:
        d = QtWidgets.QFileDialog.getExistingDirectory(self, "選擇資料夾")
        if d:
            edit.setText(d)

    def _pick_file(self, edit) -> None:
        f, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "選擇 python.exe", "", "Python (python*.exe);;All files (*.*)")
        if f:
            edit.setText(f)

    def _set_enabled(self, on: bool) -> None:
        for w in (self.btn_scan, self.btn_saveref, self.spn_rep,
                  self.cmb_lamp, self.cmb_ref):
            w.setEnabled(on)

    # ------------------------------------------------------------- connect
    def toggle_connect(self) -> None:
        from ..devices import innospectra as isc
        if self._dev is not None:
            self.disconnect_device()
            return
        sdk = self.edit_sdk.text().strip()
        py = self.edit_py.text().strip() or None
        if not os.path.isdir(sdk):
            QtWidgets.QMessageBox.information(
                self, "請先設定 SDK",
                "請指定解壓後的 ISC_NIRScan_PyQt-master 資料夾。")
            return
        self._settings.setValue("isc/sdk_dir", sdk)
        if py:
            self._settings.setValue("isc/python32", py)
        dev = isc.InnoSpectraDevice(sdk_dir=sdk, python32=py)
        try:
            dev.open()
        except DeviceError as exc:
            QtWidgets.QMessageBox.warning(self, "InnoSpectra 連線失敗", str(exc))
            return
        self._dev = dev
        self.btn_conn.setText("中斷連線")
        rng = ""
        lo, hi = dev.info.get("MinWavelength"), dev.info.get("MaxWavelength")
        if lo and hi:
            rng = f",{lo}–{hi} nm"
        self.lbl_info.setText(f"已連線 {dev.label}{rng}")
        self._set_enabled(True)
        self._status(f"Connected: {dev.label}")

    def disconnect_device(self) -> None:
        self.chk_cont.setChecked(False)
        if self._worker is not None:
            self._worker.wait(1000)
        if self._dev is not None:
            self._dev.close()
            self._dev = None
        self.btn_conn.setText("連線")
        self.lbl_info.setText("已中斷連線。")
        self._set_enabled(False)

    # ------------------------------------------------------------- scanning
    def _kinds(self):
        kinds = []
        if self.chk_abs.isChecked():
            kinds.append("absorbance")
        if self.chk_int.isChecked():
            kinds.append("intensity")
        if self.chk_refl.isChecked():
            kinds.append("reflectance")
        return tuple(kinds) or ("absorbance",)

    def do_scan(self, _=False, save_ref: bool = False) -> None:
        if self._dev is None or self._worker is not None:
            return
        try:
            self._dev.set_params(repeats=self.spn_rep.value(),
                                 lamp=self.cmb_lamp.currentData())
        except DeviceError as exc:
            QtWidgets.QMessageBox.warning(self, "設定失敗", str(exc))
            return
        self._pending_save_ref = save_ref
        ref = "new" if save_ref else self.cmb_ref.currentData()
        self._n += 1
        self.btn_scan.setEnabled(False)
        self.btn_saveref.setEnabled(False)
        self._status("InnoSpectra 掃描中…")
        self._worker = _IscScanWorker(self._dev, ref, self._kinds(),
                                      f"{self._dev.label} #{self._n}")
        self._worker.done.connect(self._on_scanned)
        self._worker.failed.connect(self._on_scan_failed)
        self._worker.finished.connect(self._scan_cleanup)
        self._worker.start()

    def do_save_ref(self) -> None:
        self.do_scan(save_ref=True)

    def _on_scanned(self, specs) -> None:
        if getattr(self, "_pending_save_ref", False) and self._dev is not None:
            try:
                self._dev.save_reference()
                self._status("已儲存為裝置的新參考。")
            except DeviceError as exc:
                QtWidgets.QMessageBox.warning(self, "儲存參考失敗", str(exc))
        for s in specs:
            self._add_spectrum(s, f"Captured {s.name}")
        if specs:
            self.curve.setData(specs[0].x, specs[0].y)
            self.plot.setLabel("left", specs[0].y_label)

    def _on_scan_failed(self, msg: str) -> None:
        self.chk_cont.setChecked(False)
        QtWidgets.QMessageBox.warning(self, "掃描失敗", msg)

    def _scan_cleanup(self) -> None:
        self._worker = None
        if self._dev is not None:
            self.btn_scan.setEnabled(True)
            self.btn_saveref.setEnabled(True)
            if self.chk_cont.isChecked():
                QtCore.QTimer.singleShot(200, self.do_scan)

    def shutdown(self) -> None:
        self.disconnect_device()


# ============================================================ dialog shell
class AcquisitionDialog(QtWidgets.QDialog):
    """Non-modal instrument window; owns both tabs and shuts them down on close."""

    def __init__(self, add_spectrum, load_paths, status, parent=None):
        super().__init__(parent)
        self.setWindowTitle("光譜儀擷取 — Ocean Optics / Go Direct / InnoSpectra")
        self.resize(760, 560)
        self.setWindowFlag(QtCore.Qt.WindowType.Window)   # own taskbar entry

        tabs = QtWidgets.QTabWidget()
        self.ocean = _OceanTab(add_spectrum, status)
        self.godirect = _GoDirectTab(load_paths, status)
        self.innospectra = _InnoSpectraTab(add_spectrum, status)
        tabs.addTab(self.ocean, "Ocean Optics (USB)")
        tabs.addTab(self.godirect, "Vernier Go Direct")
        tabs.addTab(self.innospectra, "InnoSpectra NIR")

        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(tabs)

    def closeEvent(self, ev: QtGui.QCloseEvent) -> None:  # noqa: N802
        self.ocean.shutdown()
        self.godirect.shutdown()
        self.innospectra.shutdown()
        super().closeEvent(ev)
