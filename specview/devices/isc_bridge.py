#!/usr/bin/env python
"""Bridge process for the InnoSpectra NIRScan Python SDK (32-bit only).

InnoSpectra's official Python SDK (ISC_NIRScan_PyQt: ``SDK/iscpy.pyd``) is a
CPython extension compiled for **32-bit Python 3.11 on Windows**. SpectraView
normally runs on 64-bit Python, so this script is executed as a *separate*
32-bit interpreter and speaks JSON-lines over stdin/stdout:

    python32 isc_bridge.py <path-to-ISC_NIRScan_PyQt-master>

Request  -> one JSON object per line:  {"cmd": "...", ...}
Response -> one JSON object per line:  {"ok": true, ...} | {"ok": false, "error": "..."}

Commands
    hello                          -> {"ok", "python", "bits"}
    open                           -> device info
    set   {repeats, lamp, pga}     -> apply scan parameters
    scan  {ref}                    -> full scan result arrays
                                      ref: 0 factory / 1 previous / 2 new reference
    save_ref                       -> store last ref-scan (after scan with ref=2)
    close                          -> release device

This file is also imported by ``specview.devices.innospectra`` for the rare
case where SpectraView itself already runs on 32-bit Python 3.11 (then the
SDK is used in-process and no subprocess is needed). It therefore has NO
dependencies on specview or numpy — standard library only.
"""
from __future__ import annotations

import json
import os
import struct
import sys


class ISCSession:
    """Thin wrapper around the InnoSpectra SDK Device/Scan objects."""

    LAMP = {"off": 0, "on": 1, "auto": 2}

    def __init__(self, sdk_root: str):
        self.sdk_root = os.path.abspath(sdk_root)
        self.device = None
        self.scan = None
        self._opened = False

    # ------------------------------------------------------------- helpers
    def _import_sdk(self):
        if self.device is not None:
            return
        if not os.path.isdir(os.path.join(self.sdk_root, "SDK")):
            raise RuntimeError(
                "SDK folder not found: %s (expected the ISC_NIRScan_PyQt-master "
                "directory that contains SDK/iscpy.pyd)" % self.sdk_root)
        if self.sdk_root not in sys.path:
            sys.path.insert(0, self.sdk_root)
        try:
            from SDK import device as _device, scan as _scan
        except ImportError as exc:
            bits = struct.calcsize("P") * 8
            raise RuntimeError(
                "Could not import InnoSpectra SDK (iscpy): %s\n"
                "This interpreter is Python %s (%d-bit); the SDK needs "
                "32-bit Python 3.11 on Windows." %
                (exc, sys.version.split()[0], bits)) from exc
        self.device = _device.Device()
        self.scan = _scan.Scan()

    # ------------------------------------------------------------- commands
    def open(self) -> dict:
        self._import_sdk()
        self.device.Init()
        ret = self.device.Open()
        if ret < 0:
            raise RuntimeError(
                "Device.Open() failed (ret=%d). Check the USB cable and close "
                "the ISC NIRScan GUI (only one program may hold the device)." % ret)
        self._opened = True
        self.device.Information()
        info = dict(self.device.DevInfo)
        # keep JSON clean: only simple fields
        return {k: v for k, v in info.items()
                if isinstance(v, (str, int, float))}

    def set_params(self, repeats=None, lamp=None, pga=None) -> None:
        self._require_open()
        if repeats is not None:
            if self.scan.SetScanNumRepeats(int(repeats)) < 0:
                raise RuntimeError("SetScanNumRepeats(%r) failed" % repeats)
        if lamp is not None:
            code = self.LAMP.get(str(lamp))
            if code is None:
                raise RuntimeError("lamp must be one of %s" % list(self.LAMP))
            if self.scan.SetLamp(code) < 0:
                raise RuntimeError("SetLamp(%r) failed" % lamp)
        if pga is not None:
            if pga == "auto":
                self.scan.SetFixedPGAGain(False, 0)
            else:
                self.scan.SetFixedPGAGain(True, int(pga))

    def perform_scan(self, ref: int = 1) -> dict:
        self._require_open()
        ret = self.scan.PerformScan(int(ref))
        if ret < 0:
            raise RuntimeError("PerformScan failed (ret=%d)" % ret)
        ret = self.scan.GetScanResult()
        if ret < 0:
            raise RuntimeError("GetScanResult failed (ret=%d)" % ret)
        est = None
        try:
            est = self.scan.GetEstimatedScanTime()
        except Exception:
            pass
        return {
            "wavelength": list(self.scan.WaveLength),
            "intensity": list(self.scan.Intensity),
            "absorbance": list(self.scan.Absorbance),
            "reflectance": list(self.scan.Reflectance),
            "reference": list(self.scan.ReferenceIntensity),
            "estimated_s": est,
        }

    def save_reference(self) -> None:
        self._require_open()
        if self.scan.SaveReferenceScan() < 0:
            raise RuntimeError("SaveReferenceScan failed")

    def close(self) -> None:
        if self.device is not None and self._opened:
            try:
                self.device.Close()
            except Exception:
                pass
            self._opened = False

    def _require_open(self) -> None:
        if not self._opened:
            raise RuntimeError("device not open (send {\"cmd\": \"open\"} first)")


# ------------------------------------------------------------------ stdio loop
def _serve(sdk_root: str) -> None:
    session = ISCSession(sdk_root)
    out = sys.stdout
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            req = json.loads(line)
            cmd = req.get("cmd")
            if cmd == "hello":
                resp = {"ok": True, "python": sys.version.split()[0],
                        "bits": struct.calcsize("P") * 8}
            elif cmd == "open":
                resp = {"ok": True, "info": session.open()}
            elif cmd == "set":
                session.set_params(repeats=req.get("repeats"),
                                   lamp=req.get("lamp"), pga=req.get("pga"))
                resp = {"ok": True}
            elif cmd == "scan":
                resp = {"ok": True, **session.perform_scan(req.get("ref", 1))}
            elif cmd == "save_ref":
                session.save_reference()
                resp = {"ok": True}
            elif cmd == "close":
                session.close()
                resp = {"ok": True, "bye": True}
            else:
                resp = {"ok": False, "error": "unknown cmd: %r" % cmd}
        except Exception as exc:  # noqa: BLE001 — everything goes to the client
            resp = {"ok": False, "error": str(exc)}
        out.write(json.dumps(resp) + "\n")
        out.flush()
        if resp.get("bye"):
            break
    session.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("usage: python isc_bridge.py <ISC_NIRScan_PyQt-master dir>\n")
        sys.exit(2)
    _serve(sys.argv[1])
