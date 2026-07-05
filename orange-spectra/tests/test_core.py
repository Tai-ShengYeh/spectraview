"""Headless tests for orangespectra.core (no Orange / Qt needed).

Run:  python tests/test_core.py
"""
import json
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from orangespectra import core  # noqa: E402

PASS = FAIL = 0


def check(name, cond):
    global PASS, FAIL
    PASS, FAIL = PASS + cond, FAIL + (not cond)
    print(f"  {'ok  ' if cond else 'FAIL'} {name}")


print("== source resolution ==")
check("bare id -> IRUG URL",
      core.resolve_source(4119) == "http://www.irug.org/jcamp-details?id=4119")
check("URL passthrough", core.resolve_source("https://x.org/a.jdx") == "https://x.org/a.jdx")
try:
    core.resolve_source("nope")
    check("garbage raises", False)
except ValueError:
    check("garbage raises", True)

print("== IRUG jqPlot page ==")
# Real IRUG format: jqPlotData.series['Submitter'] = {"<wavenumber>":<intensity>,...}
# (wavenumber is the quoted key). Config options are word:number and unquoted,
# so they never match the quoted-key pattern.
_irug = ("<html><body><script>var jqPlotData={};jqPlotData.series=[];"
         "jqPlotData.series['Submitter'] = "
         '{"1900.0":0.10,"1899.0":0.12,"1898.0":0.20,"1897.0":0.15,"1896.0":0.05};'
         "$.jqplot('c',[jqPlotData.series],{axes:{xaxis:{min:100,max:1900}}});"
         "</script></body></html>")
s = core.parse_irug_jqplot(_irug, source="http://www.irug.org/jcamp-details?id=4119")
check("5 points, config not counted", s is not None and s["x"].size == 5)
check("x sorted ascending", np.all(np.diff(s["x"]) > 0))
check("pairing kept (y@1900=0.10)",
      abs(np.interp(1900, s["x"], s["y"]) - 0.10) < 1e-9)
check("named from id", s["name"] == "IRUG 4119")
check("cm-1 label", "cm-1" in s["x_label"])
# Multiple series (sample + a reference match): take the submitted one.
_irug_multi = ("<script>jqPlotData.series['ref A'] = {\"10.0\":9,\"11.0\":9,\"12.0\":9};"
               "jqPlotData.series['Submitter'] = {\"10.0\":1,\"11.0\":2,\"12.0\":3};</script>")
sm = core.parse_irug_jqplot(_irug_multi, source="id=1")
check("multi-series: picks Submitter", sm["x"].size == 3 and abs(np.interp(12, sm["x"], sm["y"]) - 3) < 1e-9)

print("== SOPRANO Dygraph page ==")
_sop = ("<h1>PR1</h1><h2>Raman 785nm</h2><script>g=new Dygraph(el,"
        "[[101.0,0],[102.0,4.5],[103.0,-1.25]],"
        "{labels:[\"Raman shift\",\"PR1\"],xlabel:\"Raman shift\"});</script>")
s2 = core.parse_soprano(_sop, source="https://soprano.kikirpa.be/?id=PR1")
check("3 points", s2 is not None and s2["x"].size == 3)
check("y values", abs(s2["y"][1] - 4.5) < 1e-9 and abs(s2["y"][2] + 1.25) < 1e-9)
check("named from labels", s2["name"] == "PR1")

print("== JCAMP (AFFN) ==")
_jdx = ("##TITLE=t\n##XUNITS=1/CM\n##FIRSTX=4000\n##LASTX=3996\n##NPOINTS=5\n"
        "##XFACTOR=1\n##YFACTOR=0.01\n##XYPOINTS=(XY..XY)\n"
        "4000,10\n3999,20\n3998,30\n3997,25\n3996,15\n##END=\n")
s3 = core.parse_jcamp(_jdx, source="https://x.org/t.jdx")
check("pairs parsed, YFACTOR applied",
      s3["x"].size == 5 and abs(max(s3["y"]) - 0.30) < 1e-9)
_jdx2 = ("##TITLE=t2\n##XUNITS=1/CM\n##FIRSTX=100\n##LASTX=104\n##NPOINTS=5\n"
         "##XYDATA=(X++(Y..Y))\n100 1 2 3\n103 4 5\n##END=\n")
s3b = core.parse_jcamp(_jdx2)
check("AFFN XYDATA parsed (5 y values on linspace)",
      s3b["x"].size == 5 and abs(s3b["x"][-1] - 104) < 1e-9 and s3b["y"][-1] == 5)

print("== CSV ==")
s4 = core.parse_csv("wavenumber,intensity\n400,1\n401,2\n402,4\n",
                    source="https://x.org/pb15.csv")
check("csv parsed with header", s4["x"].size == 3 and s4["x_label"] == "wavenumber")
check("csv named from file", s4["name"] == "pb15")

print("== load_spectrum_url dispatch (injected fetch) ==")
s5 = core.load_spectrum_url(4119, fetch=lambda u: (_irug.encode(), "text/html"))
check("IRUG page via URL", s5["name"] == "IRUG 4119")
s6 = core.load_spectrum_url("https://soprano.kikirpa.be/?id=PR1",
                            fetch=lambda u: (_sop.encode(), "text/html"))
check("SOPRANO page via URL", s6["name"] == "PR1")
try:
    core.load_spectrum_url("https://x.org/none",
                           fetch=lambda u: (b"<html>no data</html>", "text/html"))
    check("no-data raises", False)
except ValueError:
    check("no-data raises", True)

print("== similarity ==")
x = np.linspace(0, 100, 400)
g = lambda c, w: np.exp(-((x - c) ** 2) / (2 * w * w))  # noqa: E731
ya = g(30, 4) + 0.6 * g(70, 6)
yb = 2.5 * ya                       # scaled copy -> perfectly similar
yc = g(50, 5)                       # different spectrum
sc_same = core.similarity_scores(x, ya, x, yb)
sc_diff = core.similarity_scores(x, ya, x, yc)
check("scaled copy: correlation ~1", abs(sc_same["correlation"] - 1) < 1e-9)
check("scaled copy: cosine ~1, sam ~0",
      abs(sc_same["cosine"] - 1) < 1e-9 and sc_same["sam"] < 1e-4)
check("different spectrum scores lower",
      sc_diff["correlation"] < 0.5 and sc_diff["euclid"] > sc_same["euclid"])
check("no overlap -> zero/inf",
      core.similarity_scores([0, 1], [1, 1], [5, 6], [1, 1])["correlation"] == 0.0)

print("== library (SpectraView .speclib round-trip) ==")
lib = [core.make_spectrum(x, ya, name="A", x_label="wavenumber (cm-1)"),
       core.make_spectrum(x, yc, name="C", x_label="wavenumber (cm-1)")]
p = os.path.join(tempfile.mkdtemp(), "t.speclib")
core.save_library(lib, p, name="test")
obj = json.load(open(p, encoding="utf-8"))
check("speclib JSON schema (SpectraView-compatible)",
      "library" in obj and obj["library"][0]["x_unit"] == "cm-1"
      and "x" in obj["library"][0] and "y" in obj["library"][0])
back = core.load_library(p)
check("round-trip 2 entries", len(back) == 2 and back[0]["name"] == "A")
check("round-trip values", np.allclose(back[1]["y"], yc))
hits = core.search_library(core.make_spectrum(x, ya + 0.01, name="q"), back)
check("search ranks the true match first", hits[0]["name"] == "A")
check("sam ranking also finds A",
      core.search_library(core.make_spectrum(x, ya, name="q"), back,
                          rank_by="sam")[0]["name"] == "A")

print("== mixture NNLS ==")
refs = [core.make_spectrum(x, g(30, 4), name="R1"),
        core.make_spectrum(x, g(70, 6), name="R2"),
        core.make_spectrum(x, g(50, 5), name="R3")]
mix = core.make_spectrum(x, 3.0 * g(30, 4) + 1.0 * g(70, 6) + 0.05, name="mix")
res = core.mixture_nnls(mix, refs, fit_offset=True)
check("recovers coefficients 3:1",
      abs(res["coeffs"][0] - 3.0) < 0.05 and abs(res["coeffs"][1] - 1.0) < 0.05)
check("absent component ~0", res["coeffs"][2] < 0.02)
check("fractions sum to 1", abs(res["fractions"].sum() - 1.0) < 1e-9)
check("offset recovered", abs(res["offset"] - 0.05) < 0.02)
check("R^2 ~ 1", res["r_squared"] > 0.999)
try:
    core.mixture_nnls(mix, [])
    check("empty refs raises", False)
except ValueError:
    check("empty refs raises", True)

print("== merge ==")
gx, ys = core.merge_spectra([core.make_spectrum([0, 1, 2], [1, 2, 3]),
                             core.make_spectrum([0.5, 1.5, 2.5], [5, 6, 7])])
check("merge overlap grid", gx[0] >= 0.5 and gx[-1] <= 2.0 and len(ys) == 2)

print("== aquagram (aquaphotomics) ==")
check("12 standard WAMACs", len(core.WAMACS) == 12)
# NIR spectra over the water region; two "groups" differing at one water band.
wl = np.linspace(1300, 1600, 300)
base = 0.5 + 0.1 * np.exp(-((wl - 1440) ** 2) / (2 * 30 ** 2))
specs_a = [core.make_spectrum(wl, base + 0.02 * i + rng_noise, name=f"A{i}",
                              x_label="wavelength (nm)")
           for i, rng_noise in enumerate([0.0, 0.001, -0.001])]
bumped = base + 0.05 * np.exp(-((wl - 1492) ** 2) / (2 * 8 ** 2))
specs_b = [core.make_spectrum(wl, bumped + 0.02 * i, name=f"B{i}",
                              x_label="wavelength (nm)") for i in range(3)]
allspec = specs_a + specs_b

raw = core.aquagram_coordinates(allspec, normalization="raw")
check("raw: n×12 matrix", raw["values"].shape == (6, 12))
check("raw: values are the absorbance (~0.5)", 0.4 < raw["values"][0].mean() < 0.8)
check("raw: WAMACs covered by 1300-1600 range", raw["covered"])

aq = core.aquagram_coordinates(allspec, normalization="aquagram")
check("aquagram: each band standardized (col mean ~0)",
      np.allclose(aq["values"].mean(axis=0), 0, atol=1e-9))
check("aquagram: each band unit std",
      np.allclose(aq["values"].std(axis=0), 1, atol=1e-9))
# 1492 nm band separates group B (bumped) from A -> its column splits by sign
i1492 = list(core.WAMACS).index(1492.0)
colA = aq["values"][:3, i1492].mean()
colB = aq["values"][3:, i1492].mean()
check("aquagram: 1492nm band separates the two groups", colB - colA > 1.0)

snv = core.aquagram_coordinates(allspec, normalization="snv")
check("snv: n×12 matrix, not standardized across set",
      snv["values"].shape == (6, 12)
      and not np.allclose(snv["values"].mean(axis=0), 0, atol=1e-6))

check("custom bands honored",
      core.aquagram_coordinates(allspec, wamacs=[1400, 1450, 1500])
      ["values"].shape == (6, 3))
try:
    core.aquagram_coordinates(allspec, normalization="nope")
    check("bad normalization raises", False)
except ValueError:
    check("bad normalization raises", True)
# WAMACs outside the measured range -> covered flag False
oor = core.aquagram_coordinates([core.make_spectrum(np.linspace(1400, 1450, 50),
                                 np.ones(50), name="narrow")])
check("out-of-range WAMACs flagged (covered=False)", not oor["covered"])


print("== peak finding ==")
_px = np.linspace(0, 100, 2001)
_py = (np.exp(-((_px - 20) / 2.0) ** 2) + 0.6 * np.exp(-((_px - 55) / 1.5) ** 2)
       + 0.3 * np.exp(-((_px - 80) / 3.0) ** 2))
_pk = core.find_spectrum_peaks(_px, _py, min_height_frac=0.1,
                               min_prominence_frac=0.05)
check("finds the 3 peaks", len(_pk) == 3)
check("centers correct", all(abs(p["center"] - c) < 0.2
                             for p, c in zip(_pk, (20, 55, 80))))
check("gaussian FWHM correct (2*sqrt(ln2)*sigma)",
      abs(_pk[0]["fwhm"] - 2 * np.sqrt(np.log(2)) * 2.0) < 0.15)
check("sorted by center", _pk == sorted(_pk, key=lambda p: p["center"]))
check("height threshold filters",
      len(core.find_spectrum_peaks(_px, _py, min_height_frac=0.7)) == 1)
check("min_distance merges",
      len(core.find_spectrum_peaks(_px, _py, min_height_frac=0.1,
                                   min_prominence_frac=0.05,
                                   min_distance=50)) < 3)
check("smoothing tolerated",
      len(core.find_spectrum_peaks(_px, _py, min_height_frac=0.1,
                                   min_prominence_frac=0.05,
                                   smooth_window=7)) == 3)
check("no peaks in flat line",
      core.find_spectrum_peaks(_px, np.ones_like(_px)) == [])

print("== XRF identification ==")
from orangespectra import xrf  # noqa: E402
check("53 elements Na-U", len(xrf.ELEMENTS) == 53)
_m = xrf.identify_energy(6.40, tol=0.05)
check("6.40 keV -> Fe Ka1", _m and _m[0]["symbol"] == "Fe" and _m[0]["line"] == "Ka1")
_m = xrf.identify_energy(10.55, tol=0.05)
check("10.55 keV -> Pb La1 (closest first)", _m[0]["symbol"] == "Pb")
check("As Ka1 also within tol at 10.54",
      any(x["symbol"] == "As" for x in xrf.identify_energy(10.55, tol=0.05)))
check("line filter works",
      all(x["line"] in ("Ka1", "Kb1")
          for x in xrf.identify_energy(10.55, tol=0.2, line_filter={"Ka1", "Kb1"})))
_r = xrf.identify_peaks([6.404, 8.048, 99.0], tol=0.05)
check("identify_peaks best per energy",
      _r[0]["best"]["symbol"] == "Fe" and _r[1]["best"]["symbol"] == "Cu"
      and _r[2]["best"] is None)

print("== PLS-DA ==")
_rs = np.random.RandomState(1)
_grid = np.linspace(0, 1, 60)
_X, _labs = [], []
for _ci, _c in enumerate((0.25, 0.5, 0.75)):
    for _ in range(10):
        _X.append(np.exp(-((_grid - _c) / 0.05) ** 2) * (1 + 0.1 * _rs.randn())
                  + 0.02 * _rs.randn(60))
        _labs.append(f"class{_ci}")
_res = core.plsda_fit(np.array(_X), _labs, n_components=3)
check("separable classes -> 100% training accuracy", _res["accuracy"] == 1.0)
check("scores shape (n, A)", _res["scores"].shape == (30, 3))
check("loadings shape (p, A)", _res["loadings"].shape == (60, 3))
check("VIP length p, mean(VIP^2)=1",
      _res["vip"].shape == (60,) and abs((_res["vip"] ** 2).mean() - 1.0) < 1e-6)
check("confusion diagonal", _res["confusion"].trace() == 30)
check("explained X variance <= 1 and decreasing-ish",
      _res["explained_x_variance"].sum() <= 1.0 + 1e-9)
_sub = [0, 1, 10, 11, 20, 21]     # 2 samples from each of the 3 classes
check("components clipped to <= n-1",
      core.plsda_fit(np.array(_X)[_sub], [_labs[i] for i in _sub],
                     n_components=99)["n_components"] <= 5)
try:
    core.plsda_fit(np.array(_X), ["same"] * 30, n_components=2)
    check("single class raises", False)
except ValueError:
    check("single class raises", True)


print("== embedded JS chart parsing ==")
_plotly = ('<script>var t={x:[400,401,402,403,404],y:[0.1,0.4,0.9,0.5,0.2]};'
           'Plotly.newPlot("d",[t]);</script>')
_sp = core.load_spectrum_url("https://infra-art.eu/sample/PSE320",
                             fetch=lambda u: (_plotly.encode(), "text/html"))
check("Plotly x/y arrays parsed", _sp is not None and _sp["x"].size == 5)
check("named from URL last segment", _sp["name"] == "PSE320")
check("y pairing kept", abs(np.interp(402, _sp["x"], _sp["y"]) - 0.9) < 1e-9)
_hc = ('<script>Highcharts.chart("c",{series:[{data:[[4000,0.02],[3999,0.03],'
       '[3998,0.08],[3997,0.05],[3996,0.01]]}]});</script>')
_sp2 = core.load_spectrum_url("https://x.org/s/AB12",
                              fetch=lambda u: (_hc.encode(), "text/html"))
check("Highcharts pair array parsed", _sp2["x"].size == 5)
check("sorted ascending", _sp2["x"][0] == 3996.0 and np.all(np.diff(_sp2["x"]) > 0))
_cj = ('<script>new Chart(x,{data:{labels:[1000,1001,1002,1003,1004,1005],'
       'datasets:[{data:[5,6,9,7,6,5]}]}});</script>')
_sp3 = core.load_spectrum_url("https://x.org/c/CJ",
                              fetch=lambda u: (_cj.encode(), "text/html"))
check("Chart.js labels+data parsed", _sp3["x"].size == 6)
check("embedded parser does not hijack IRUG",
      core.load_spectrum_url(4119, fetch=lambda u: (
          '<script>jqPlotData.series[\'Submitter\'] = '
          '{"1900.0":0.1,"1899.0":0.2,"1898.0":0.3};</script>'.encode(),
          "text/html"))["name"] == "IRUG 4119")
try:
    core.load_spectrum_url("https://x.org/none",
                           fetch=lambda u: (b"<html>nothing here</html>", "text/html"))
    check("no-chart page raises", False)
except ValueError:
    check("no-chart page raises", True)


print("== files: folder / zip / matrix / NetCDF ==")
from orangespectra import files as _files  # noqa: E402
_fd = tempfile.mkdtemp()
open(os.path.join(_fd, "a.csv"), "w").write("wl,i\n400,1\n401,2\n402,4\n")
open(os.path.join(_fd, "b.dx"), "w").write(
    "##TITLE=B\n##XUNITS=1/CM\n##FIRSTX=100\n##LASTX=104\n##NPOINTS=5\n"
    "##XYDATA=(X++(Y..Y))\n100 1 2 3\n103 4 5\n##END=\n")
open(os.path.join(_fd, "wide.csv"), "w").write(
    "name,1000,1001,1002,1003\nS1,1,2,3,4\nS2,4,3,2,1\n")
open(os.path.join(_fd, "long.csv"), "w").write(
    "wl,alpha,beta\n500,1,9\n501,2,8\n502,3,7\n503,4,6\n")
open(os.path.join(_fd, "junk.md"), "w").write("not a spectrum")
from scipy.io import netcdf_file as _ncf
_nc = _ncf(os.path.join(_fd, "apple.cdf"), "w")
_nc.createDimension("sample", 3); _nc.createDimension("wavelength", 12)
_wl = _nc.createVariable("wavelength", "d", ("wavelength",))
_wl[:] = np.linspace(400, 411, 12)
_ab = _nc.createVariable("absorbance", "d", ("sample", "wavelength"))
_ab[:] = np.random.RandomState(0).rand(3, 12)
_nc.close()

_wide = _files.load_spectra_path(os.path.join(_fd, "wide.csv"))
check("wide matrix: 2 named spectra",
      len(_wide) == 2 and _wide[0]["name"] == "S1"
      and _wide[0]["x"][0] == 1000.0)
_long = _files.load_spectra_path(os.path.join(_fd, "long.csv"))
check("long matrix: 2 column spectra",
      len(_long) == 2 and {_s["name"] for _s in _long} == {"alpha", "beta"})
_ncs = _files.load_netcdf(os.path.join(_fd, "apple.cdf"))
check("NetCDF: 3 rows, wavelength axis",
      len(_ncs) == 3 and _ncs[0]["x_label"] == "wavelength"
      and _ncs[0]["x"][0] == 400.0)
_all = _files.load_spectra_folder(_fd)
check("folder skips junk, loads all spectra", len(_all) == 1 + 1 + 2 + 2 + 3)
import zipfile as _zf
_zp = os.path.join(_fd, "pack.zip")
with _zf.ZipFile(_zp, "w") as _z:
    for _fn in ("a.csv", "wide.csv", "apple.cdf"):
        _z.write(os.path.join(_fd, _fn), _fn)
check("zip archive loads without extraction",
      len(_files.load_spectra_zip(_zp)) == 1 + 2 + 3)
_nc2 = _ncf(os.path.join(_fd, "pair.cdf"), "w")
_nc2.createDimension("n", 15)
_tv = _nc2.createVariable("time", "d", ("n",)); _tv[:] = np.arange(15.0)
_yv = _nc2.createVariable("tic", "d", ("n",))
_yv[:] = np.random.RandomState(1).rand(15)
_nc2.close()
_p = _files.load_netcdf(os.path.join(_fd, "pair.cdf"))
check("NetCDF 1-D pair fallback", len(_p) == 1 and _p[0]["x_label"] == "time")


print("== spectrometer (image -> spectrum) ==")
from orangespectra import spectrometer as _spm  # noqa: E402
_H, _W = 40, 640
_img = np.zeros((_H, _W, 3))
for _col, _rgb in [(100, (0.2, 0.4, 1.0)), (500, (1.0, 0.3, 0.1))]:
    _band = np.exp(-((np.arange(_W) - _col) / 4.0) ** 2)
    for _c, _ch in enumerate(_rgb):
        _img[:, :, _c] += _band * _ch * 255
_img = np.clip(_img, 0, 255)
_prof = _spm.extract_profile(_img, channel="luminance")
check("profile length = width", _prof.size == _W)
check("emission lines are profile maxima",
      _prof[100] == _prof[95:106].max() and _prof[500] == _prof[495:506].max())
check("red channel favours the red line, blue the blue line",
      _spm.extract_profile(_img, "red")[500] > _spm.extract_profile(_img, "blue")[500]
      and _spm.extract_profile(_img, "blue")[100] > _spm.extract_profile(_img, "red")[100])
_c, _r2, _deg = _spm.fit_calibration([100, 500], [435.8, 611.6], "linear")
check("linear calibration R^2 ~ 1 (degree 1)", _r2 > 0.999 and _deg == 1)
_sp = _spm.image_to_spectrum(_img, calibration=[(100, 435.8), (500, 611.6)])
check("calibrated -> wavelength axis, ascending",
      _sp["x_label"] == "wavelength (nm)" and np.all(np.diff(_sp["x"]) > 0))
check("line pixel maps to its wavelength",
      abs(_sp["x"][np.argmin(np.abs(_sp["x"] - 435.8))] - 435.8) < 1.0)
_raw = _spm.image_to_spectrum(_img, calibration=None)
check("no calibration -> pixel axis",
      _raw["x_label"] == "pixel" and _raw["x"][0] == 0 and _raw["x"][-1] == _W - 1)
_q = _spm.image_to_spectrum(_img, calibration=[(50, 400), (300, 530), (600, 680)],
                            model="quadratic")
check("quadratic calibration uses degree 2", _q["calibration"]["degree"] == 2)
_flip = _spm.image_to_spectrum(_img, calibration=None, flip=True)
check("flip reverses the profile", np.allclose(_flip["y"][::-1], _raw["y"]))
try:
    _spm.fit_calibration([1], [2])
    check("too few calibration points raises", False)
except ValueError:
    check("too few calibration points raises", True)

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)
