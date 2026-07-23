"""Vernier Go Direct devices (SpectroVis Plus etc.) via the ``godirect`` package.

IMPORTANT LIMITATION (verified against Vernier's official docs, 2026):
    Vernier's public ``godirect`` Python library explicitly does **not**
    support Go Direct *spectrometers* — the full-spectrum data stream of the
    Go Direct SpectroVis Plus (GDX-SVISPL) uses a private protocol that only
    Vernier's own "Spectral Analysis" app implements.
    https://github.com/VernierST/godirect-examples  ("Go Direct spectrometers
    ... do not work with the godirect library.")

What this module therefore offers:

1. :func:`list_devices` / :class:`GoDirectDevice` — a best-effort godirect
   connection. For ordinary Go Direct sensors this streams real values; for
   the SpectroVis Plus it can connect and identify the unit, but no spectrum
   channels are exposed, and a clear explanation is raised.

2. :func:`is_spectral_analysis_export` + the *watch-folder* workflow in the
   acquisition dialog: run Vernier **Spectral Analysis** alongside SpectraView,
   export (or auto-save) CSV files into a folder, and SpectraView imports each
   new file the moment it appears. This is the reliable way to get SpectroVis
   Plus spectra into SpectraView today.

Install (one-time)::

    pip install godirect
    # USB needs hidapi, BLE needs bleak — both are pulled in by godirect.
"""
from __future__ import annotations

import os

from ..formats.binary_io import MissingDependency
from .base import DeviceError, DeviceInfo, DeviceNotFound

_INSTALL_HINT = (
    "Vernier Go Direct support needs the 'godirect' package.\n\n"
    "    pip install godirect\n\n"
    "USB uses hidapi; Bluetooth uses bleak (installed automatically)."
)

SPECTRO_NOTE = (
    "Vernier 的公開 godirect 程式庫不支援 Go Direct 光譜儀的光譜資料串流 "
    "(官方文件明載 spectrometers 不適用)。\n"
    "建議做法:開啟 Vernier「Spectral Analysis」程式取光譜並匯出 CSV,"
    "SpectraView 會監看匯出資料夾並自動載入新檔案。"
)


def _gd():
    try:
        from godirect import GoDirect
        return GoDirect
    except ImportError as exc:
        raise MissingDependency(_INSTALL_HINT) from exc


def list_devices(use_usb: bool = True, use_ble: bool = False,
                 ble_timeout: int = 8) -> list[DeviceInfo]:
    """Enumerate Go Direct devices over USB and/or Bluetooth LE."""
    GoDirect = _gd()
    gd = GoDirect(use_usb=use_usb, use_ble=use_ble)
    out = []
    try:
        devices = gd.list_devices() or []
        for d in devices:
            name = getattr(d, "name", None) or str(d)
            out.append(DeviceInfo(
                backend="godirect",
                ident=name,
                label=name,
                extra={"type": getattr(d, "type", "?")},
            ))
    finally:
        try:
            gd.quit()
        except Exception:  # noqa: BLE001
            pass
    return out


class GoDirectDevice:
    """Best-effort Go Direct connection (see module docstring for limits).

    Not a :class:`SpectrometerDevice` — Go Direct exposes named scalar
    channels, not a wavelength axis. The dialog shows live channel values.
    """

    backend = "godirect"

    def __init__(self, name: str | None = None,
                 use_usb: bool = True, use_ble: bool = False):
        self._name = name
        self._use_usb = use_usb
        self._use_ble = use_ble
        self._gd = None
        self._dev = None
        self.label = name or "Go Direct"

    # -- lifecycle -----------------------------------------------------------
    def open(self) -> None:
        if self._dev is not None:
            return
        GoDirect = _gd()
        self._gd = GoDirect(use_usb=self._use_usb, use_ble=self._use_ble)
        dev = None
        try:
            devices = self._gd.list_devices() or []
            if self._name:
                for d in devices:
                    if getattr(d, "name", "") == self._name:
                        dev = d
                        break
            elif devices:
                dev = devices[0]
            if dev is None:
                raise DeviceNotFound(
                    "找不到 Go Direct 裝置。\n"
                    "請確認 USB 線或藍牙 (裝置藍燈閃爍表示等待配對),"
                    "並關閉其他佔用裝置的程式 (Spectral Analysis / Graphical Analysis)。")
            if not dev.open(auto_start=False):
                raise DeviceError("Go Direct 裝置開啟失敗 (可能被其他程式佔用)。")
        except (DeviceError, MissingDependency):
            self._quit()
            raise
        except Exception as exc:  # noqa: BLE001
            self._quit()
            raise DeviceError(f"Go Direct 連線失敗: {exc}") from exc
        self._dev = dev
        self.label = getattr(dev, "name", None) or self.label

    def close(self) -> None:
        if self._dev is not None:
            try:
                self._dev.stop()
                self._dev.close()
            except Exception:  # noqa: BLE001
                pass
            self._dev = None
        self._quit()

    def _quit(self) -> None:
        if self._gd is not None:
            try:
                self._gd.quit()
            except Exception:  # noqa: BLE001
                pass
            self._gd = None

    @property
    def is_open(self) -> bool:
        return self._dev is not None

    # -- introspection ----------------------------------------------------------
    def is_spectrometer(self) -> bool:
        """SpectroVis-family units identify themselves by name (GDX-SVISPL...)."""
        n = (self.label or "").upper()
        return "SVIS" in n or "SPECTRO" in n

    def sensor_names(self) -> list[str]:
        """Names of all channels the godirect protocol exposes on this unit."""
        if self._dev is None:
            raise DeviceError("裝置尚未連線。")
        try:
            sensors = self._dev.list_sensors() or {}
            out = []
            for s in sensors.values():
                desc = getattr(s, "sensor_description", None) or str(s)
                units = getattr(s, "sensor_units", "") or ""
                out.append(f"{desc} ({units})" if units else desc)
            return out
        except Exception as exc:  # noqa: BLE001
            raise DeviceError(f"讀取感測器清單失敗: {exc}") from exc

    # -- scalar streaming (works for ordinary Go Direct sensors) ------------------
    def start(self, period_ms: int = 250) -> None:
        if self._dev is None:
            raise DeviceError("裝置尚未連線。")
        self._dev.enable_default_sensors()
        self._dev.start(period=period_ms)

    def read_values(self) -> dict[str, float]:
        """Return {channel description: latest value} for one reading."""
        if self._dev is None:
            raise DeviceError("裝置尚未連線。")
        out: dict[str, float] = {}
        if self._dev.read():
            for s in self._dev.get_enabled_sensors() or []:
                vals = getattr(s, "values", None)
                if vals:
                    desc = getattr(s, "sensor_description", "sensor")
                    out[desc] = float(vals[-1])
                    s.clear()
        return out

    def stop(self) -> None:
        if self._dev is not None:
            try:
                self._dev.stop()
            except Exception:  # noqa: BLE001
                pass


# ---------------------------------------------------------------------------
# Spectral Analysis export helpers (the practical SpectroVis Plus route)
# ---------------------------------------------------------------------------

_SA_EXTS = {".csv", ".txt"}


def is_spectral_analysis_export(path: str) -> bool:
    """Cheap check: does this look like a spectrum file we should auto-import?"""
    return os.path.splitext(path)[1].lower() in _SA_EXTS
