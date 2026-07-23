"""InnoSpectra NIRScan (NIR-S-G1 / NIR-M-T1 ...) via the official Python SDK.

The SDK (``ISC_NIRScan_PyQt``, BSD-licensed by InnoSpectra) ships a compiled
extension ``SDK/iscpy.pyd`` built for **32-bit Python 3.11 on Windows**.
Because SpectraView usually runs on 64-bit Python, this backend talks to the
SDK through a small bridge subprocess (``isc_bridge.py``) executed with a
32-bit interpreter. If the *current* interpreter happens to be 32-bit Python
3.11, the SDK is used in-process instead — same API either way.

Setup (one-time)
    1. Keep the unzipped ``ISC_NIRScan_PyQt-master`` folder somewhere stable.
    2. Install 32-bit Python 3.11 from python.org (choose "Windows installer
       (32-bit)"); no extra packages are needed for the bridge.
    3. In SpectraView: 儀器 ▸ 光譜儀擷取 ▸ InnoSpectra tab, point to both paths.

A scan takes seconds (repeats × exposure); acquisition is therefore
scan-based rather than streaming.
"""
from __future__ import annotations

import json
import os
import struct
import subprocess
import sys

import numpy as np

from ..spectrum import Spectrum
from .base import DeviceError

#: y-data types a scan can produce -> (SDK field, SpectraView y_unit)
RESULT_KINDS = {
    "absorbance": ("absorbance", "absorbance"),
    "intensity": ("intensity", "counts"),
    "reflectance": ("reflectance", "reflectance"),
    "reference": ("reference", "counts"),
}

REF_CHOICES = {"factory": 0, "previous": 1, "new": 2}


def _default_python32() -> str | None:
    """Best-guess locations of a 32-bit Python 3.11 on Windows."""
    for p in (
        os.path.expandvars(r"%LOCALAPPDATA%\Programs\Python\Python311-32\python.exe"),
        r"C:\Python311-32\python.exe",
        r"C:\Program Files (x86)\Python311-32\python.exe",
    ):
        if os.path.isfile(p):
            return p
    return None


class InnoSpectraDevice:
    """One InnoSpectra NIRScan unit, driven through the official SDK.

    ``sdk_dir``  – the unzipped ``ISC_NIRScan_PyQt-master`` folder
    ``python32`` – path to a 32-bit Python 3.11 interpreter (None = try
                   in-process import, which only works if SpectraView itself
                   runs on 32-bit Python 3.11)
    """

    backend = "innospectra"

    def __init__(self, sdk_dir: str, python32: str | None = None):
        self.sdk_dir = sdk_dir
        self.python32 = python32
        self._proc: subprocess.Popen | None = None
        self._session = None          # in-process ISCSession
        self.info: dict = {}
        self.label = "InnoSpectra NIRScan"

    # ------------------------------------------------------------ lifecycle
    def open(self) -> None:
        if self.is_open:
            return
        if self._can_inprocess():
            from . import isc_bridge
            self._session = isc_bridge.ISCSession(self.sdk_dir)
            self.info = self._session.open()
        else:
            py = self.python32 or _default_python32()
            if not py or not os.path.isfile(py):
                raise DeviceError(
                    "找不到 32 位元 Python 3.11。InnoSpectra 的 iscpy.pyd 只能在 "
                    "32-bit Python 3.11 下執行:\n"
                    "1) 到 python.org 下載安裝 Python 3.11 (Windows installer, 32-bit)\n"
                    "2) 在 InnoSpectra 頁籤指定該 python.exe 路徑")
            bridge = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  "isc_bridge.py")
            flags = 0x08000000 if os.name == "nt" else 0   # CREATE_NO_WINDOW
            try:
                self._proc = subprocess.Popen(
                    [py, bridge, self.sdk_dir],
                    stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, text=True, encoding="utf-8",
                    creationflags=flags)
            except OSError as exc:
                raise DeviceError(f"無法啟動橋接程式: {exc}") from exc
            hello = self._rpc({"cmd": "hello"}, timeout=20)
            if hello.get("bits") != 32:
                self.close()
                raise DeviceError(
                    f"指定的 Python 是 {hello.get('bits')} 位元 "
                    f"(v{hello.get('python')}); iscpy.pyd 需要 32 位元 3.11。")
            self.info = self._rpc({"cmd": "open"}, timeout=30).get("info", {})
        model = self.info.get("ModelName") or "NIRScan"
        serial = self.info.get("SerialNumber") or "?"
        self.label = f"InnoSpectra {model} ({serial})"

    def close(self) -> None:
        if self._session is not None:
            try:
                self._session.close()
            except Exception:  # noqa: BLE001
                pass
            self._session = None
        if self._proc is not None:
            try:
                if self._proc.poll() is None:
                    self._send_line({"cmd": "close"})
                    self._proc.wait(timeout=5)
            except Exception:  # noqa: BLE001
                self._proc.kill()
            self._proc = None

    @property
    def is_open(self) -> bool:
        return self._session is not None or (
            self._proc is not None and self._proc.poll() is None)

    # ------------------------------------------------------------ commands
    def set_params(self, repeats: int | None = None, lamp: str | None = None,
                   pga: int | str | None = None) -> None:
        """repeats: scans to average; lamp: off/on/auto; pga: gain or 'auto'."""
        if self._session is not None:
            self._session.set_params(repeats=repeats, lamp=lamp, pga=pga)
        else:
            self._rpc({"cmd": "set", "repeats": repeats,
                       "lamp": lamp, "pga": pga}, timeout=15)

    def scan(self, ref: str = "previous",
             kinds: tuple[str, ...] = ("absorbance",),
             name: str | None = None) -> list[Spectrum]:
        """Perform one scan; return one Spectrum per requested kind.

        ref: 'factory' (出廠參考) / 'previous' (前次參考) / 'new' (本次當新參考)
        """
        refcode = REF_CHOICES.get(ref)
        if refcode is None:
            raise DeviceError(f"ref 必須是 {list(REF_CHOICES)},收到 {ref!r}")
        if self._session is not None:
            res = self._session.perform_scan(refcode)
        else:
            res = self._rpc({"cmd": "scan", "ref": refcode}, timeout=300)
        x = np.asarray(res.get("wavelength") or [], dtype=float)
        if x.size == 0:
            raise DeviceError("掃描沒有回傳波長資料。")
        out = []
        base = name or self.label
        for kind in kinds:
            field, y_unit = RESULT_KINDS.get(kind, (None, None))
            y = np.asarray(res.get(field) or [], dtype=float) if field else None
            if y is None or y.size != x.size:
                continue
            out.append(Spectrum(
                x.copy(), y, name=f"{base} [{kind}]", x_unit="nm", y_unit=y_unit,
                meta={"source": "live:innospectra", "device": self.label,
                      "reference": ref, **{k: v for k, v in self.info.items()
                                           if k in ("ModelName", "SerialNumber")}}))
        if not out:
            raise DeviceError("掃描結果中沒有所要求的資料類型。")
        return out

    def save_reference(self) -> None:
        """After scan(ref='new'), persist it as the device's local reference."""
        if self._session is not None:
            self._session.save_reference()
        else:
            self._rpc({"cmd": "save_ref"}, timeout=15)

    # ------------------------------------------------------------ plumbing
    def _can_inprocess(self) -> bool:
        return (os.name == "nt" and struct.calcsize("P") * 8 == 32
                and sys.version_info[:2] == (3, 11))

    def _send_line(self, obj: dict) -> None:
        assert self._proc and self._proc.stdin
        self._proc.stdin.write(json.dumps(obj) + "\n")
        self._proc.stdin.flush()

    def _rpc(self, obj: dict, timeout: float = 60) -> dict:
        if self._proc is None or self._proc.poll() is not None:
            raise DeviceError("橋接程式未啟動或已結束。")
        import threading
        self._send_line(obj)
        result: dict = {}

        def _read():
            line = self._proc.stdout.readline()
            if line:
                try:
                    result.update(json.loads(line))
                except json.JSONDecodeError:
                    result.update({"ok": False, "error": f"壞回應: {line[:200]}"})

        t = threading.Thread(target=_read, daemon=True)
        t.start()
        t.join(timeout)
        if t.is_alive():
            self._proc.kill()
            self._proc = None
            raise DeviceError(f"橋接程式逾時 (> {timeout:g}s),已中止。")
        if not result:
            err = ""
            try:
                err = self._proc.stderr.read() or ""
            except Exception:  # noqa: BLE001
                pass
            self._proc = None
            raise DeviceError(f"橋接程式意外結束。{err[:400]}")
        if not result.get("ok"):
            raise DeviceError(result.get("error", "未知錯誤"))
        return result

    # ------------------------------------------------------------ sugar
    def __enter__(self):
        self.open()
        return self

    def __exit__(self, *exc):
        self.close()
        return False
