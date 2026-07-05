"""Real NNLS mixture example on the Eigenvector CGL NIR dataset.

CGL = casein / glucose / lactate — a full three-component NIR mixture design
(231 samples, 117 wavelengths 1104-2496 nm, targets casein/glucose/lactate/
moisture wt%; Naes & Isaksson). Ideal for the Mixture Analysis (NNLS) widget:

  1. estimate each component's PURE NIR spectrum from the calibration set by
     classical least squares (CLS):  Spure = pinv(C) . X   (C = wt% of the 3
     solids, X = spectra),
  2. take a TEST mixture and unmix it with non-negative least squares
     (orangespectra.core.mixture_nnls) against those 3 pure spectra,
  3. compare the recovered proportions with the reference wt% values.

It also writes example files so you can reproduce it live in Orange:
  examples/cgl_components.speclib   (3 pure-component reference spectra)
  examples/cgl_mixture.csv          (one test mixture spectrum)
  docs/demo_mixture.png             (the figure used in the tutorial)

CGL_nir.mat (~a few hundred KB) is looked up locally, else downloaded from the
Eigenvector archive. Run:  python examples/cgl_mixture_nnls.py [path/to/CGL_nir.mat]
"""
import pathlib
import sys

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "orange-spectra"))
from orangespectra.core import make_spectrum, mixture_nnls, save_library

COMPONENTS = ["casein", "glucose", "lactate"]      # the 3 solids (ignore moisture)

# ------------------------------------------------------------ locate CGL_nir.mat
CANDIDATES = [ROOT / "examples" / "CGL_nir.mat", pathlib.Path.cwd() / "CGL_nir.mat"]
if len(sys.argv) > 1:
    CANDIDATES.insert(0, pathlib.Path(sys.argv[1]))
MIRRORS = ["https://eigenvector.com/data/CGL/CGL_nir.mat",
           "https://www.eigenvector.com/data/CGL/CGL_nir.mat"]

mat_path = next((p for p in CANDIDATES if p.exists()), None)
if mat_path is None:
    import urllib.request
    dst = ROOT / "examples" / "CGL_nir.mat"
    for url in MIRRORS:
        try:
            print(f"Downloading CGL_nir.mat from {url} ...")
            dst.write_bytes(urllib.request.urlopen(url, timeout=60).read())
            mat_path = dst
            break
        except Exception as exc:                       # noqa: BLE001
            print(f"  failed: {exc}")
    if mat_path is None:
        sys.exit("Could not obtain CGL_nir.mat — download it from "
                 "https://eigenvector.com/resources/data-sets/ and pass its path.")


def dset(mat, name):
    """Return (data, axisscale) from an Eigenvector DataSet stored in a .mat."""
    obj = mat[name]
    data = np.asarray(getattr(obj, "data", obj), dtype=float)
    axis = None
    try:                                                # axisscale{2} = columns
        ax = obj.axisscale
        ax = ax[1] if isinstance(ax, (list, np.ndarray)) and len(ax) > 1 else ax
        axis = np.asarray(getattr(ax, "values", ax), dtype=float).ravel()
    except Exception:                                   # noqa: BLE001
        pass
    return data, axis


mat = sio.loadmat(str(mat_path), squeeze_me=True, struct_as_record=False)
have = [k for k in mat if not k.startswith("__")]
print("variables:", have)

Xcal, wl = dset(mat, "Xcal")
Ycal, _ = dset(mat, "Ycal")
Xtest, _ = dset(mat, "Xtest")
Ytest, _ = dset(mat, "Ytest")
if wl is None or wl.size != Xcal.shape[1]:
    wl = np.linspace(1104, 2496, Xcal.shape[1])         # documented axis fallback
C = Ycal[:, :3]                                         # casein, glucose, lactate wt%

# ------------------------------------------------------------ CLS pure spectra
Spure = np.linalg.pinv(C) @ Xcal                        # 3 x n_wl
refs = [make_spectrum(wl, Spure[i], name=COMPONENTS[i],
                      x_label="wavelength (nm)") for i in range(3)]

# ------------------------------------------------------------ NNLS on a test mix
# pick the most balanced 3-component test mixture for a clear demo
frac_test = Ytest[:, :3] / Ytest[:, :3].sum(1, keepdims=True)
k = int(np.argmin(np.abs(frac_test - 1 / 3).sum(1)))
mix = make_spectrum(wl, Xtest[k], name="CGL test mixture", x_label="wavelength (nm)")
res = mixture_nnls(mix, refs, fit_offset=True)

true_frac = frac_test[k]
print("\nTest mixture #%d  (casein/glucose/lactate wt%%): %s"
      % (k, np.round(Ytest[k, :3], 1)))
for name, rec, tru in zip(COMPONENTS, res["fractions"], true_frac):
    print(f"  {name:8s}: NNLS {rec*100:5.1f}%   reference {tru*100:5.1f}%")
print("NNLS fit R^2 = %.4f" % res["r_squared"])

# ------------------------------------------------------------ figure
fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 5),
                               gridspec_kw=dict(width_ratios=[2, 1]))
axL.plot(wl, mix["y"], color="#222", lw=1.8, label="test mixture")
axL.plot(res["x"], res["fit"], color="#c44e52", lw=1.4, ls="--", label="NNLS fit")
cols = ["#4c72b0", "#dd8452", "#55a868"]
for i, name in enumerate(COMPONENTS):
    axL.plot(wl, res["coeffs"][i] * Spure[i], color=cols[i], lw=1.0, alpha=0.9,
             label=f"{name} × {res['coeffs'][i]:.1f}")
axL.set_xlabel("wavelength (nm)"); axL.set_ylabel("absorbance")
axL.set_title(f"CGL NIR mixture unmixed by NNLS  (R²={res['r_squared']:.4f})",
              fontsize=11)
axL.legend(fontsize=8)

x = np.arange(3); w = 0.38
axR.bar(x - w / 2, true_frac * 100, w, color="#9ec9d6", label="reference wt%")
axR.bar(x + w / 2, res["fractions"] * 100, w, color="#0E7C7B", label="NNLS %")
axR.set_xticks(x); axR.set_xticklabels(COMPONENTS)
axR.set_ylabel("proportion (%)"); axR.set_title("recovered vs reference", fontsize=11)
axR.legend(fontsize=8)
for xi, (t, r) in enumerate(zip(true_frac, res["fractions"])):
    axR.text(xi - w / 2, t * 100 + 1, f"{t*100:.0f}", ha="center", fontsize=8)
    axR.text(xi + w / 2, r * 100 + 1, f"{r*100:.0f}", ha="center", fontsize=8)
fig.suptitle("Mixture Analysis (NNLS) on real CGL NIR — casein / glucose / lactate",
             fontsize=13, color="#119C9A")
fig.tight_layout()
fig.savefig(ROOT / "docs" / "demo_mixture.png", dpi=140, facecolor="white",
            bbox_inches="tight")
print("saved", ROOT / "docs" / "demo_mixture.png")

# ------------------------------------------------------------ example files
save_library(refs, str(ROOT / "examples" / "cgl_components.speclib"),
             name="CGL pure components (CLS)")
with open(ROOT / "examples" / "cgl_mixture.csv", "w") as fh:
    fh.write("wavelength,absorbance\n")
    for a, b in zip(wl, Xtest[k]):
        fh.write(f"{a:.1f},{b:.6f}\n")
print("wrote examples/cgl_components.speclib and examples/cgl_mixture.csv")
