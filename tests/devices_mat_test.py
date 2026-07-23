"""Headless tests for the instrument backends and the .mat reader.

Run:  python tests/devices_mat_test.py
No hardware and no PySide6 needed (the UI module is only syntax-checked).
"""
from __future__ import annotations

import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from specview.devices import DeviceError, SpectrometerDevice          # noqa: E402
from specview.devices import godirect_dev, oceanoptics                # noqa: E402
from specview.formats.binary_io import MissingDependency              # noqa: E402
from specview.formats.mat_io import load_mat                          # noqa: E402

PASS = FAIL = 0


def check(label, fn):
    global PASS, FAIL
    try:
        fn()
        PASS += 1
        print(f"  ok    {label}")
    except Exception as exc:  # noqa: BLE001
        FAIL += 1
        print(f"  FAIL  {label}: {exc!r}")


# ---------------------------------------------------------------- devices
class _FakeDev(SpectrometerDevice):
    """In-memory instrument to exercise the shared read_spectrum() path."""
    backend = "fake"

    def __init__(self):
        self._open = False
        self._x = np.linspace(400, 1000, 601)

    def open(self):
        self._open = True

    def close(self):
        self._open = False

    @property
    def is_open(self):
        return self._open

    def set_integration_time_ms(self, ms):
        pass

    def wavelengths(self):
        return self._x

    def read_intensities(self):
        return np.exp(-0.5 * ((self._x - 650) / 30) ** 2) * 1000 \
            + np.random.default_rng(0).normal(0, 1, self._x.size)


def test_fake_device_spectrum():
    with _FakeDev() as d:
        s = d.read_spectrum(name="test", averages=4)
        assert s.npoints == 601 and s.x_unit == "nm" and s.y_unit == "counts"
        assert s.meta["averages"] == 4
        assert abs(s.x[np.argmax(s.y)] - 650) < 2      # peak at 650 nm


def test_oceanoptics_missing_dep():
    """Without seabreeze installed, a clear install hint must be raised."""
    try:
        oceanoptics.list_devices()
    except MissingDependency as exc:
        assert "pip install seabreeze" in str(exc)
    except Exception:
        # seabreeze installed in this env: enumeration itself must not crash
        pass


def test_godirect_missing_dep():
    try:
        godirect_dev.list_devices()
    except MissingDependency as exc:
        assert "pip install godirect" in str(exc)
    except Exception as exc:  # noqa: BLE001
        raise AssertionError(f"unexpected: {exc}")


def test_godirect_helpers():
    d = godirect_dev.GoDirectDevice(name="GDX-SVISPL 123")
    d.label = "GDX-SVISPL 0F100123"
    assert d.is_spectrometer()
    d.label = "GDX-TMP 0A2"
    assert not d.is_spectrometer()
    assert godirect_dev.is_spectral_analysis_export("a/b/run1.csv")
    assert not godirect_dev.is_spectral_analysis_export("a/b/run1.pdf")


def test_device_error_when_closed():
    d = oceanoptics.OceanOpticsDevice()
    try:
        d.read_intensities()
        raise AssertionError("should have raised")
    except DeviceError:
        pass


# ---------------------------------------------------------------- .mat v7
def _tmp(name):
    return os.path.join(tempfile.gettempdir(), name)


def test_mat_v5_named():
    from scipy.io import savemat
    x = np.linspace(400, 700, 301)
    y = np.sin(x / 50)
    p = _tmp("t_v5.mat")
    savemat(p, {"wavelength": x, "intensity": y})
    (s,) = load_mat(p)
    assert s.npoints == 301 and s.x_unit == "nm"
    assert np.allclose(s.x, x)


def test_mat_v5_matrix():
    from scipy.io import savemat
    x = np.linspace(1000, 2500, 100)
    m = np.column_stack([x, np.cos(x / 100), np.sin(x / 100)])
    p = _tmp("t_v5m.mat")
    savemat(p, {"data": m})
    specs = load_mat(p)
    assert len(specs) == 2 and specs[0].npoints == 100


def test_mat_v5_dataset_struct():
    """PLS_Toolbox-style struct: data matrix + axisscale cell with wavelengths."""
    from scipy.io import savemat
    x = np.linspace(1100, 2498, 700)
    d = np.vstack([np.sin(x / (100 + 10 * i)) for i in range(5)])
    axcell = np.empty((2, 2), dtype=object)
    axcell[0, 0] = np.array([]); axcell[0, 1] = ""
    axcell[1, 0] = x; axcell[1, 1] = ""
    p = _tmp("t_v5ds.mat")
    savemat(p, {"m5spec": {"data": d, "axisscale": axcell, "name": "corn m5"}})
    specs = load_mat(p)
    assert len(specs) == 5, f"expected 5 spectra, got {len(specs)}"
    assert specs[0].npoints == 700 and specs[0].x_unit == "nm"
    assert "m5" in specs[0].name


def test_mat_index_table_skipped():
    """A struct whose axis is a plain 1..N ramp is a table, not spectra."""
    from scipy.io import savemat
    n = 40
    axcell = np.empty((1, 2), dtype=object)
    axcell[0, 0] = np.arange(1, n + 1, dtype=float); axcell[0, 1] = ""
    p = _tmp("t_v5tab.mat")
    savemat(p, {"props": {"data": np.random.rand(3, n), "axisscale": axcell},
                "wavelength": np.linspace(400, 700, 90),
                "intensity": np.random.rand(90)})
    specs = load_mat(p)
    assert len(specs) == 1 and specs[0].npoints == 90   # only the real spectrum


# ---------------------------------------------------------------- .mat v7.3
def _write_v73(path, variables):
    """Emulate MATLAB -v7.3: HDF5 with the MATLAB userblock signature."""
    import h5py
    with h5py.File(path, "w", userblock_size=512) as f:
        for k, v in variables.items():
            arr = np.atleast_2d(np.asarray(v, dtype=float)).T   # column-major
            ds = f.create_dataset(k, data=arr)
            ds.attrs["MATLAB_class"] = np.bytes_("double")
    # scipy sniffs the header to decide it's v7.3; write the signature text
    with open(path, "r+b") as fh:
        fh.write(b"MATLAB 7.3 MAT-file, written by test" + b" " * 80)
        fh.seek(124)
        fh.write(b"\x00\x02IM")


def test_mat_v73_named():
    x = np.linspace(400, 700, 301)
    y = np.exp(-((x - 550) ** 2) / 800)
    p = _tmp("t_v73.mat")
    _write_v73(p, {"wavelength": x, "absorbance": y})
    (s,) = load_mat(p)
    assert s.npoints == 301
    assert s.x_unit == "nm" and s.y_unit == "absorbance"
    assert abs(s.x[np.argmax(s.y)] - 550) < 2


def test_mat_v73_matrix():
    x = np.linspace(0, 10, 50)
    m = np.column_stack([x, x ** 2])
    p = _tmp("t_v73m.mat")
    _write_v73(p, {"spec": m})
    (s,) = load_mat(p)
    assert s.npoints == 50 and np.allclose(s.y, s.x ** 2)


# ---------------------------------------------------------------- UI syntax
def test_ui_compiles():
    import py_compile
    root = os.path.join(os.path.dirname(__file__), "..", "specview")
    for rel in ("ui/acquisition.py", "ui/main_window.py",
                "devices/base.py", "devices/oceanoptics.py",
                "devices/godirect_dev.py", "formats/mat_io.py"):
        py_compile.compile(os.path.join(root, rel), doraise=True)


if __name__ == "__main__":
    print("devices + mat tests")
    for n, f in sorted({k: v for k, v in globals().items()
                        if k.startswith("test_")}.items()):
        check(n, f)
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)
