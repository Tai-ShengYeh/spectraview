"""Regenerate docs/demo_aquagram.png from the REAL Eigenvector corn NIR dataset.

Corn `m5spec` (80 samples, 1100-2498 nm) fully covers the 12 WAMACs water bands
(1342-1516 nm). Samples are split by reference moisture into low vs high groups;
the classic normalized aquagram (left) separates them at the water bands, next
to SNV-only (right).

The dataset (corn.mat, ~1.4 MB) is looked up locally, else downloaded from a
public mirror. Run:  python docs/regen_aquagram_demo.py [path/to/corn.mat]
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
from orangespectra.core import WAMACS, aquagram_coordinates

# corn.mat lookup: CLI arg, common local spots, else download from a mirror.
CANDIDATES = [
    ROOT / "orange-spectra" / "validation" / "corn.mat",
    ROOT / "docs" / "corn.mat",
    pathlib.Path.cwd() / "corn.mat",
]
if len(sys.argv) > 1:
    CANDIDATES.insert(0, pathlib.Path(sys.argv[1]))
MIRRORS = [
    "https://raw.githubusercontent.com/ryuzakyl/data-bloodhound/master/datasets/nir_corn/data/corn.mat",
    "https://www.eigenvector.com/data/Corn/corn.mat",
    "https://eigenvector.com/wp-content/uploads/2019/06/corn.mat_.zip",
]

mat_path = next((p for p in CANDIDATES if p.exists()), None)
if mat_path is None:
    import io
    import urllib.request
    dst = ROOT / "orange-spectra" / "validation" / "corn.mat"
    for url in MIRRORS:
        try:
            print(f"Downloading corn.mat from {url} ...")
            raw = urllib.request.urlopen(url, timeout=60).read()
            if url.endswith(".zip"):
                import zipfile
                with zipfile.ZipFile(io.BytesIO(raw)) as z:
                    raw = z.read("corn.mat")
            dst.write_bytes(raw)
            mat_path = dst
            break
        except Exception as exc:                       # noqa: BLE001
            print(f"  failed: {exc}")
    if mat_path is None:
        sys.exit("Could not obtain corn.mat — pass a local path as the first argument.")

mat = sio.loadmat(str(mat_path), squeeze_me=True, struct_as_record=False)
wl = np.arange(1100, 2500, 2).astype(float)
X = mat["m5spec"].data.astype(float)                    # 80 x 700 absorbance
moisture = mat["propvals"].data.astype(float)[:, 0]     # column 0 = moisture

order = np.argsort(moisture)
low_idx, high_idx = order[:20], order[-20:]
print(f"moisture: low mean={moisture[low_idx].mean():.2f}%  "
      f"high mean={moisture[high_idx].mean():.2f}%")
specs = ([{"x": wl, "y": X[i], "name": f"L{i}"} for i in low_idx] +
         [{"x": wl, "y": X[i], "name": f"H{i}"} for i in high_idx])
n_low = len(low_idx)
V = aquagram_coordinates(specs, normalization="aquagram")["values"]
snv = aquagram_coordinates(specs, normalization="snv")["values"]


def radar(ax, coords, title, rlim, n_low):
    m = len(WAMACS)
    ang = np.linspace(0, 2 * np.pi, m, endpoint=False)
    angc = np.concatenate([ang, ang[:1]])
    for i, row in enumerate(coords):
        col = "#1f77b4" if i < n_low else "#E36414"
        ax.plot(angc, np.concatenate([row, row[:1]]), lw=0.7, color=col,
                alpha=0.32, zorder=2)
    for grp, col, lab in ((coords[:n_low], "#1f77b4", "Low moisture (mean)"),
                          (coords[n_low:], "#E36414", "High moisture (mean)")):
        mm = np.concatenate([grp.mean(0), grp.mean(0)[:1]])
        ax.plot(angc, mm, lw=2.6, color=col, marker="o", ms=3, label=lab, zorder=4)
        ax.fill(angc, mm, color=col, alpha=0.08, zorder=1)
    ax.plot(angc, np.zeros_like(angc), color="#666", lw=1.0, ls="--", zorder=1)
    ax.set_xticks(ang)
    ax.set_xticklabels([str(int(l)) for l in WAMACS], fontsize=8.5)
    ax.tick_params(axis="x", pad=6)
    ax.set_ylim(rlim)
    lo, hi = int(np.floor(rlim[0])), int(np.ceil(rlim[1]))
    ax.set_yticks(list(range(lo, hi + 1)))
    ax.set_rlabel_position(15)                          # off the 1342 spoke
    ax.tick_params(axis="y", labelsize=7, colors="#8a8a8a")
    ax.set_title(title, color="#0b5c9e", fontsize=11, pad=18)
    ax.grid(alpha=.35)


rmax = float(np.max(np.abs(np.concatenate([V.ravel(), snv.ravel()])))) * 1.10
rlim = (-rmax, rmax)
fig, axes = plt.subplots(1, 2, figsize=(12, 6), subplot_kw=dict(polar=True))
radar(axes[0], V, "Aquagram (SNV + across-sample std)\nreal corn NIR, low vs high moisture",
      rlim, n_low)
radar(axes[1], snv, "SNV-only sampled at WAMACs\n(no across-sample standardization)",
      rlim, n_low)
axes[0].legend(loc="lower center", bbox_to_anchor=(0.5, -0.24), ncol=2,
               fontsize=9, frameon=False)
fig.suptitle("Aquagram on real corn NIR spectra (Eigenvector m5, low vs high moisture)",
             fontsize=13, color="#119C9A", y=1.03)
fig.subplots_adjust(wspace=0.4)
out = ROOT / "docs" / "demo_aquagram.png"
fig.savefig(out, dpi=140, facecolor="white", bbox_inches="tight")
print("saved", out)
