"""Headless self-test for the non-GUI core (no Qt needed).

Run:  python tests/smoke_test.py
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from specview import axes, processing  # noqa: E402
from specview.demo import load_demo_set  # noqa: E402
from specview.formats import load_any, save_csv, save_jcamp  # noqa: E402
from specview.formats.jcamp import load_jcamp  # noqa: E402
from specview.spectrum import Spectrum  # noqa: E402

PASS, FAIL = 0, 0


def check(name, cond):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  ok   {name}")
    else:
        FAIL += 1
        print(f"  FAIL {name}")


print("== model & demo ==")
specs = load_demo_set()
check("demo set has 5 spectra", len(specs) == 5)
check("spectrum x sorted ascending", all(np.all(np.diff(s.x) > 0) for s in specs))
ftir = specs[0]
check("ftir is cm-1 / absorbance", ftir.x_unit == "cm-1" and ftir.y_unit == "absorbance")

print("== axis conversions ==")
nm = axes.convert_x(ftir, "nm")
back = axes.convert_x(nm, "cm-1")
check("cm-1 -> nm -> cm-1 round trip", np.allclose(np.sort(back.x), np.sort(ftir.x), rtol=1e-6))
check("nm conversion flips/keeps monotonic", np.all(np.diff(nm.x) > 0))
ab = specs[2]  # uv-vis absorbance
t = axes.convert_y(ab, "transmittance")
ab2 = axes.convert_y(t, "absorbance")
check("A -> T -> A round trip", np.allclose(ab2.y, ab.y, atol=1e-6))
check("transmittance in (0,1]", np.all(t.y > 0) and np.all(t.y <= 1.0 + 1e-9))
raman = specs[1]
nm_r = axes.convert_x(raman, "nm", laser_nm=785.0)
back_r = axes.convert_x(nm_r, "raman_cm-1", laser_nm=785.0)
check("raman shift -> nm -> raman round trip",
      np.allclose(np.sort(back_r.x), np.sort(raman.x), atol=1e-3))

print("== processing ==")
sm = processing.savitzky_golay(ftir, 21, 3)
check("SG smoothing keeps length", sm.npoints == ftir.npoints)
check("SG reduces noise std", np.std(np.diff(sm.y)) < np.std(np.diff(ftir.y)))
d1 = processing.derivative(ftir, 1, 21, 3)
check("1st derivative length", d1.npoints == ftir.npoints)
rb = processing.baseline_rubberband(ftir)
check("rubberband baseline -> min near 0", rb.y.min() >= -1e-6)
als = processing.baseline_als(ftir, 1e5, 0.01)
check("ALS baseline runs", np.isfinite(als.y).all())
poly = processing.baseline_polynomial(ftir, 3)
check("poly baseline runs", np.isfinite(poly.y).all())
nmax = processing.normalize(ftir, "max")
check("normalize max -> peak 1", abs(np.max(np.abs(nmax.y)) - 1.0) < 1e-9)
narea = processing.normalize(ftir, "area")
_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
check("normalize area -> integral 1", abs(_trapz(np.abs(narea.y), narea.x) - 1.0) < 1e-6)
snv = processing.snv(ftir)
check("SNV -> mean~0 std~1", abs(snv.y.mean()) < 1e-9 and abs(snv.y.std() - 1) < 1e-9)
det = processing.detrend(ftir, 1)
check("detrend runs", np.isfinite(det.y).all())
cutk = processing.cut(ftir, 1500, 1800, keep=True)
check("cut keeps subrange", cutk.x.min() >= 1500 - 1e-9 and cutk.x.max() <= 1800 + 1e-9)
ipl = processing.interpolate(ftir, n=512)
check("interpolate to 512 pts", ipl.npoints == 512)
msc = processing.msc(specs[3], specs[3].copy())
check("msc runs", np.isfinite(msc.y).all())

print("== arithmetic ==")
a = Spectrum(np.linspace(0, 10, 101), np.ones(101) * 2, name="A")
b = Spectrum(np.linspace(0, 10, 101), np.ones(101) * 3, name="B")
check("A + B = 5", np.allclose(processing.combine(a, b, "add").y, 5))
check("A - B = -1", np.allclose(processing.combine(a, b, "sub").y, -1))
check("A * B = 6", np.allclose(processing.combine(a, b, "mul").y, 6))
avg = processing.average([a, b], "mean")
check("mean(A,B) = 2.5", np.allclose(avg.y, 2.5))

print("== SpectrumSet identity with mixed lengths (regression) ==")
# Spectrum must use identity equality: dataclass __eq__ would compare numpy
# arrays during `in`/index/remove and crash on differing lengths.
from specview.spectrum import SpectrumSet  # noqa: E402

doc = SpectrumSet()
s_short = Spectrum(np.linspace(0, 10, 200), np.ones(200), name="short", y_unit="%R")
s_long = Spectrum(np.linspace(0, 10, 257), np.full(257, 5.0), name="long", y_unit="%R")
doc.add(s_short)            # index 0, length 200
doc.add(s_long)             # index 1, length 257
check("Spectrum uses identity equality", (s_long == s_long) and not (s_long == s_short))
try:
    # replacing the 2nd spectrum scans the list -> compares against the 200-pt one
    doc.replace(s_long, axes.convert_y(s_long, "absorbance"))
    doc.remove(doc[0])
    check("replace/remove across mismatched lengths (no broadcast crash)", len(doc) == 1)
except Exception as exc:  # noqa: BLE001
    check(f"mixed-length membership crashed: {exc}", False)

print("== analysis: peaks / fit / integrate ==")
from specview import analysis as A  # noqa: E402

xx = np.linspace(400, 800, 2000)
def _g(c, amp, s):  # noqa: E306
    return amp * np.exp(-(xx - c) ** 2 / (2 * s ** 2))
syn = Spectrum(xx, _g(500, 1.0, 8) + _g(650, 0.6, 12), name="syn",
               x_unit="nm", y_unit="absorbance")
pks = A.find_peaks(syn, 0.1, 0.05)
check("find_peaks finds 2 peaks", len(pks) == 2)
check("peak centers ~500/650",
      abs(pks[0].center - 500) < 1 and abs(pks[1].center - 650) < 1)
check("peak FWHM ~18.8/28.3",
      abs(pks[0].fwhm - 18.8) < 1.5 and abs(pks[1].fwhm - 28.3) < 2)
fit = A.fit_peaks(syn, "gaussian", 2)
check("gaussian fit R2 > 0.999", fit.r_squared > 0.999)
check("fit recovers centers", abs(fit.components[0].center - 500) < 0.5
      and abs(fit.components[1].center - 650) < 0.5)
check("fit recovers areas (amp·σ·√2π)", abs(fit.components[0].area - 20.05) < 0.5
      and abs(fit.components[1].area - 18.05) < 0.5)
check("components sorted by center",
      fit.components[0].center < fit.components[1].center)
fl = A.fit_peaks(syn, "lorentzian", 2)
check("lorentzian fit runs", len(fl.components) == 2 and np.isfinite(fl.r_squared))
fpv = A.fit_peaks(syn, "pseudovoigt", 2)
check("pseudovoigt fit: eta in [0,1]", all(0 <= c.eta <= 1 for c in fpv.components))
comp = A.component_to_spectrum(syn, fit.x, fit.comp_curves[0], "c1")
check("component_to_spectrum -> Spectrum", comp.npoints == fit.x.size and comp.name == "c1")
intg = A.integrate(syn, 460, 560, "none")
check("integrate area ~ first gaussian", abs(intg["area"] - 20.05) < 0.3)
check("integrate centroid ~500", abs(intg["centroid"] - 500) < 1)
check("integrate linear baseline runs",
      np.isfinite(A.integrate(syn, 400, 800, "linear")["area"]))
vpks = A.find_peaks(Spectrum(xx, 1.0 - _g(500, 1.0, 8), x_unit="nm"),
                    0.1, 0.05, valleys=True)
check("valley detection finds the dip", len(vpks) == 1 and abs(vpks[0].center - 500) < 1)

print("== IO round trips ==")
tmp = tempfile.mkdtemp()
csv_path = os.path.join(tmp, "s.csv")
save_csv(ftir, csv_path)
r = load_any(csv_path)[0]
check("CSV round trip values", np.allclose(r.y, ftir.y, rtol=1e-4))
dx_path = os.path.join(tmp, "s.dx")
save_jcamp(ftir, dx_path)
r2 = load_any(dx_path)[0]
check("JCAMP round trip length", r2.npoints == ftir.npoints)
check("JCAMP round trip values", np.allclose(np.sort(r2.y), np.sort(ftir.y), atol=1e-3))

print("== combined export ==")
from specview.formats import save_combined_csv  # noqa: E402

comb_c = os.path.join(tmp, "comb_cols.csv")
info_c = save_combined_csv([a, b], comb_c, "columns")
check("combined columns: 2 spectra, not resampled",
      info_c == {"n_spectra": 2, "n_points": 101, "resampled": False})
back = load_any(comb_c)
check("combined columns re-imports 2 named spectra",
      len(back) == 2 and back[0].name == "A" and back[1].name == "B")
check("combined columns values preserved",
      np.allclose(back[0].y, 2.0) and np.allclose(back[1].y, 3.0))
comb_r = os.path.join(tmp, "comb_rows.csv")
save_combined_csv([a, b], comb_r, "rows")
rlines = open(comb_r, encoding="utf-8").read().splitlines()
check("combined rows: 1 header + 2 spectrum rows", len(rlines) == 3
      and rlines[1].split(",")[0] == "A" and len(rlines[1].split(",")) == 102)
c_short = Spectrum(np.linspace(0, 10, 57), np.full(57, 7.0), name="c")
info_m = save_combined_csv([a, c_short], os.path.join(tmp, "comb_mixed.csv"), "columns")
check("combined mixed-length interpolates (resampled=True)", info_m["resampled"] is True)

print("== JCAMP ASDF decoder ==")
# DIF/SQZ encoded: y = [10,12,15,15,15,13,10]
dif_text = (
    "##TITLE=DIF test\n##JCAMP-DX=4.24\n##DATA TYPE=INFRARED SPECTRUM\n"
    "##XUNITS=1/CM\n##YUNITS=ARBITRARY UNITS\n##XFACTOR=1.0\n##YFACTOR=1.0\n"
    "##FIRSTX=1\n##LASTX=7\n##NPOINTS=7\n##XYDATA=(X++(Y..Y))\n"
    "1A0KL%%kl\n##END=\n"
)
p = os.path.join(tmp, "dif.dx")
open(p, "w").write(dif_text)
s = load_jcamp(p)[0]
check("DIF/SQZ decode values", np.array_equal(s.y, [10, 12, 15, 15, 15, 13, 10]))

# DUP encoded: y = [5,5,5,5]
dup_text = (
    "##TITLE=DUP test\n##JCAMP-DX=4.24\n##XUNITS=1/CM\n##YUNITS=ARBITRARY UNITS\n"
    "##XFACTOR=1.0\n##YFACTOR=1.0\n##FIRSTX=1\n##LASTX=4\n##NPOINTS=4\n"
    "##XYDATA=(X++(Y..Y))\n1E%U\n##END=\n"
)
p2 = os.path.join(tmp, "dup.dx")
open(p2, "w").write(dup_text)
s2 = load_jcamp(p2)[0]
check("DUP decode values", np.array_equal(s2.y, [5, 5, 5, 5]))

# AFFN XYPOINTS pairs
aff_text = (
    "##TITLE=AFFN pairs\n##XUNITS=NANOMETERS\n##YUNITS=ABSORBANCE\n"
    "##XFACTOR=1.0\n##YFACTOR=1.0\n##XYPOINTS=(XY..XY)\n"
    "200, 0.1\n201, 0.2\n202, 0.35\n##END=\n"
)
p3 = os.path.join(tmp, "aff.dx")
open(p3, "w").write(aff_text)
s3 = load_jcamp(p3)[0]
check("AFFN pairs decode", np.allclose(s3.x, [200, 201, 202]) and np.allclose(s3.y, [0.1, 0.2, 0.35]))
check("AFFN pairs units", s3.x_unit == "nm" and s3.y_unit == "absorbance")

print("== airPLS baseline ==")
xb = np.linspace(0, 1000, 1000)
bg = 100.0 * np.exp(-(xb - 500) ** 2 / (2 * 300 ** 2)) + 50.0   # broad background
pk = 200.0 * np.exp(-(xb - 700) ** 2 / (2 * 10 ** 2))           # sharp peak
sb = Spectrum(xb, bg + pk, name="airpls")
corr = processing.baseline_airpls(sb, lam=1e5, porder=2)
check("airPLS removes background (median ~0)", abs(np.median(corr.y)) < 10)
check("airPLS preserves the peak", corr.y.max() > 150)
check("airPLS flattens peak-free region",
      np.abs(corr.y[xb < 600]).max() < 0.1 * corr.y.max())

print("== JSON IO ==")
import json as _json  # noqa: E402
from specview.formats import save_json  # noqa: E402

jdir = tempfile.mkdtemp()
def _wj(name, obj):  # noqa: E306
    p = os.path.join(jdir, name)
    _json.dump(obj, open(p, "w"))
    return p
sj = load_any(_wj("a.json", {"x": [1, 2, 3], "y": [4, 5, 6]}))[0]
check("JSON x/y arrays", sj.npoints == 3 and np.allclose(sj.y, [4, 5, 6]))
sj = load_any(_wj("b.json", {"name": "S", "x_unit": "nm", "y_unit": "absorbance",
                             "wavelength": [1, 2], "absorbance": [7, 8]}))[0]
check("JSON aliased keys + units",
      sj.name == "S" and sj.x_unit == "nm" and sj.y_unit == "absorbance")
sj = load_any(_wj("c.json", {"data": [[1, 9], [2, 8]]}))[0]
check("JSON data pairs", np.allclose(sj.x, [1, 2]) and np.allclose(sj.y, [9, 8]))
check("JSON spectra list -> N spectra",
      len(load_any(_wj("d.json", {"spectra": [{"x": [1, 2], "y": [1, 2]},
                                              {"x": [1, 2], "y": [3, 4]}]}))) == 2)
check("JSON top-level pairs list",
      load_any(_wj("e.json", [[1, 1], [2, 2], [3, 3]]))[0].npoints == 3)
rt = Spectrum([400, 500], [0.1, 0.2], name="rt", x_unit="nm", y_unit="absorbance")
rp = os.path.join(jdir, "rt.json")
save_json(rt, rp)
bj = load_any(rp)[0]
check("JSON round trip", np.allclose(bj.y, rt.y) and bj.x_unit == "nm" and bj.name == "rt")

print("== MATLAB .mat IO ==")
from scipy.io import savemat  # noqa: E402

mp = os.path.join(jdir, "named.mat")
savemat(mp, {"wavenumber": np.linspace(400, 4000, 50), "intensity": np.arange(50.0)})
sm = load_any(mp)[0]
check("MAT named variables", sm.npoints == 50 and sm.x_unit == "cm-1")
mp2 = os.path.join(jdir, "matrix.mat")
savemat(mp2, {"M": np.column_stack([np.arange(30.0), np.ones(30)])})
check("MAT single 2-D matrix", load_any(mp2)[0].npoints == 30)

print("== library search ==")
from specview.library import SpectralLibrary, similarity_scores  # noqa: E402

xl = np.linspace(400, 1800, 800)
def _bnd(c, a, w):  # noqa: E306
    return a * np.exp(-(xl - c) ** 2 / (2 * w ** 2))
refA = Spectrum(xl, _bnd(700, 1.0, 30) + _bnd(1200, 0.5, 40), name="A", x_unit="cm-1")
refB = Spectrum(xl, _bnd(900, 0.8, 25) + _bnd(1500, 0.6, 35), name="B", x_unit="cm-1")
lib = SpectralLibrary("t")
lib.add(refA)
lib.add(refB)
lp = os.path.join(tempfile.mkdtemp(), "t.speclib")
lib.save(lp)
lib2 = SpectralLibrary.load(lp)
check("library save/load preserves entries", len(lib2) == 2 and lib2.entries[0].name == "A")
hits = lib2.search(refA.copy(), top_n=2)
check("library search ranks the match first",
      hits[0]["name"] == "A" and hits[0]["scores"]["correlation"] > 0.99)
check("similarity of identical spectra ~ 1",
      abs(similarity_scores(refA, refA.copy())["correlation"] - 1) < 1e-9)

print("== mixture NNLS ==")
mix = Spectrum(xl, 0.6 * refA.y + 0.4 * refB.y, name="mix", x_unit="cm-1")
mres = A.mixture_nnls(mix, [refA, refB], fit_offset=True)
check("mixture fractions recovered (0.6 / 0.4)",
      abs(mres.fractions[0] - 0.6) < 0.02 and abs(mres.fractions[1] - 0.4) < 0.02)
check("mixture R^2 ~ 1", mres.r_squared > 0.999)

print("== XRF element identification ==")
from specview import xrf  # noqa: E402

ids = xrf.identify_peaks([3.692, 6.404, 8.048], tol=0.05)
check("XRF identifies Ca / Fe / Cu",
      [r["best"]["symbol"] for r in ids] == ["Ca", "Fe", "Cu"])
check("XRF line label is Greek", ids[1]["best"]["line_label"] == "Kα1")
check("XRF tolerance rejects far energies", xrf.identify_peaks([50.0], 0.05)[0]["best"] is None)

print("== 2D correlation (2D-COS / 2T2D) ==")
from specview import cos2d  # noqa: E402
from specview.demo import demo_cos_series  # noqa: E402

series = demo_cos_series(14)
xc, sync, asyn = cos2d.correlation_from_spectra(series, ref="mean")
check("2D-COS synchronous is symmetric", np.allclose(sync, sync.T))
check("2D-COS asynchronous is antisymmetric", np.allclose(asyn, -asyn.T))
i12, i15 = np.argmin(np.abs(xc - 1200)), np.argmin(np.abs(xc - 1500))
check("sync cross-peak (1200,1500) < 0 (anti-correlated)", sync[i12, i15] < 0)
check("async cross-peak (1200,1500) nonzero (sequential)", abs(asyn[i12, i15]) > 1e-3)
H = cos2d.hilbert_noda_matrix(5)
check("Hilbert-Noda zero diagonal & antisymmetric",
      np.allclose(np.diag(H), 0) and np.allclose(H, -H.T))
gx, phi, psi = cos2d.two_trace_from_spectra(series[0], series[-1])
check("2T2D Φ symmetric, Ψ antisymmetric",
      np.allclose(phi, phi.T) and np.allclose(psi, -psi.T))
check("2T2D Ψ nonzero", np.abs(psi).max() > 1e-6)

print("== fluorescence EEM ==")
from specview import eem as eemmod  # noqa: E402
from specview.demo import demo_eem  # noqa: E402

E = demo_eem()
check("EEM Z shape is ex×em", E.Z.shape == (E.ex.size, E.em.size))
check("EEM scatter removal masks Rayleigh", np.isnan(eemmod.remove_scatter(E, width=10).Z).sum() > 0)
ems = [Spectrum(np.linspace(300, 560, 180),
                np.exp(-((np.linspace(300, 560, 180) - (360 + i * 4)) ** 2) / (2 * 18 ** 2)),
                name=f"ex{260 + i * 10}") for i in range(7)]
Ef = eemmod.eem_from_spectra(ems)
check("EEM from 7 emission spectra", Ef.Z.shape == (7, 180) and Ef.ex[0] == 260)
ep = os.path.join(tempfile.mkdtemp(), "e.csv")
with open(ep, "w", encoding="utf-8") as fh:
    fh.write("corner," + ",".join(f"{v:.0f}" for v in E.ex) + "\n")
    for j, emv in enumerate(E.em):
        fh.write(f"{emv:.0f}," + ",".join(f"{E.Z[i, j]:.2f}" for i in range(E.ex.size)) + "\n")
E2 = eemmod.read_eem_matrix(ep, ex_in_columns=True)
check("EEM matrix file round-trip", np.allclose(E2.Z, E.Z, atol=1e-1))

print("== EEM PARAFAC ==")
from specview.demo import demo_eem_stack  # noqa: E402

stack = demo_eem_stack(7)
pres = eemmod.parafac_from_eems(stack, rank=3, nonneg=True)
check("PARAFAC fit > 0.95", pres.fit > 0.95)
check("PARAFAC factor shapes",
      pres.scores.shape == (7, 3) and pres.ex_load.shape == (stack[0].ex.size, 3))
em_peaks = sorted(pres.em[np.argmax(pres.em_load[:, f])] for f in range(3))
check("PARAFAC recovers emission peaks ~330/400/440",
      all(any(abs(p - t) < 12 for p in em_peaks) for t in (330, 400, 440)))
check("PARAFAC loadings are non-negative",
      (pres.ex_load >= 0).all() and (pres.em_load >= 0).all())
check("component_eem is a rank-1 EEM",
      pres.component_eem(0).Z.shape == (stack[0].ex.size, stack[0].em.size))

print("== hetero 2D-COS ==")
sA, sB = demo_cos_series(10), demo_cos_series(10)
gx1, gx2, hsync, hasyn = cos2d.hetero_from_spectra(sA, sB, "mean")
check("hetero map shape matches the two grids", hsync.shape == (gx1.size, gx2.size))
_, homo_sync, _ = cos2d.correlation_from_spectra(sA, "mean")
_, _, hsame, _ = cos2d.hetero_from_spectra(sA, sA, "mean")
check("hetero(X,X) sync == generalized sync", np.allclose(hsame, homo_sync, atol=1e-9))

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)
