"""Regenerate docs/demo_aquagram.png from a REAL NIR dataset (offline-runnable).

Uses the public chemotools 'fermentation' NIR set (21 spectra, 428-1833 nm,
covers the 12 WAMACs water bands 1342-1516 nm), split by reference glucose into
low vs high groups. Shows the classic normalized aquagram (left) vs SNV-only
(right). Run:  pip install chemotools ; python docs/regen_aquagram_demo.py
"""
import pathlib, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "orange-spectra"))
from orangespectra.core import aquagram_coordinates, WAMACS

import chemotools.datasets as d
X, y = d.load_fermentation_train()
wl = np.array([float(c) for c in X.columns])
Xv = X.values.astype(float)
gluc = y.values.astype(float).ravel()

order = np.argsort(gluc)
n = len(order); half = n // 2
low_idx, high_idx = order[:half], order[half:]      # low vs high glucose
specs = ([{"x": wl, "y": Xv[i], "name": f"low{i}"} for i in low_idx] +
         [{"x": wl, "y": Xv[i], "name": f"high{i}"} for i in high_idx])
n_low = len(low_idx)
V = aquagram_coordinates(specs, normalization="aquagram")["values"]
snv = aquagram_coordinates(specs, normalization="snv")["values"]
print(f"{n} spectra; low glucose mean={gluc[low_idx].mean():.1f}, "
      f"high={gluc[high_idx].mean():.1f} g/L")


def radar(ax, coords, title, rlim, n_low):
    m = len(WAMACS)
    ang = np.linspace(0, 2 * np.pi, m, endpoint=False)
    angc = np.concatenate([ang, ang[:1]])
    for i, row in enumerate(coords):
        col = "#1f77b4" if i < n_low else "#E36414"
        ax.plot(angc, np.concatenate([row, row[:1]]), lw=0.7, color=col,
                alpha=0.35, zorder=2)
    for grp, col, lab in ((coords[:n_low], "#1f77b4", "Low glucose (mean)"),
                          (coords[n_low:], "#E36414", "High glucose (mean)")):
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
    ax.set_rlabel_position(15)                        # off the 1342 spoke
    ax.tick_params(axis="y", labelsize=7, colors="#8a8a8a")
    ax.set_title(title, color="#0b5c9e", fontsize=11, pad=18)
    ax.grid(alpha=.35)


rmax = float(np.max(np.abs(np.concatenate([V.ravel(), snv.ravel()])))) * 1.10
rlim = (-rmax, rmax)
fig, axes = plt.subplots(1, 2, figsize=(12, 6), subplot_kw=dict(polar=True))
radar(axes[0], V, "Aquagram (SNV + across-sample std)\nreal NIR, low vs high glucose",
      rlim, n_low)
radar(axes[1], snv, "SNV-only sampled at WAMACs\n(no across-sample standardization)",
      rlim, n_low)
axes[0].legend(loc="lower center", bbox_to_anchor=(0.5, -0.24), ncol=2,
               fontsize=9, frameon=False)
fig.suptitle("Aquagram on real NIR spectra (fermentation set, low vs high glucose)",
             fontsize=13, color="#119C9A", y=1.03)
fig.subplots_adjust(wspace=0.4)
out = ROOT / "docs" / "demo_aquagram.png"
fig.savefig(out, dpi=140, facecolor="white", bbox_inches="tight")
print("saved", out)
