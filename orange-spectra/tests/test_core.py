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

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)
