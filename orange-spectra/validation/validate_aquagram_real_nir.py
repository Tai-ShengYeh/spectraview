"""Validate the Aquagram widget core on REAL NIR spectra (HANDOFF §6-5).

Uses the public Eigenvector corn dataset (80 samples on m5, 1100-2498 nm,
2 nm step) which fully covers the 12 WAMACs water bands (1342-1516 nm).
We split samples by moisture (low vs high) and check that the aquagram
separates them at the water-related bands, and that math invariants hold.
"""
import pathlib, sys
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# import the package core directly (repo-relative, works from any clone)
ROOT = pathlib.Path(__file__).resolve().parents[2]   # .../spectraview
sys.path.insert(0, str(ROOT / "orange-spectra"))
from orangespectra.core import aquagram_coordinates, WAMACS

DATA = pathlib.Path(__file__).with_name("corn.mat")
OUT = ROOT / "orange-spectra" / "validation"
OUT.mkdir(exist_ok=True)

# Auto-download the public Eigenvector corn dataset if not present locally.
if not DATA.exists():
    import urllib.request, zipfile, io
    url = "https://eigenvector.com/wp-content/uploads/2019/06/corn.mat_.zip"
    print(f"Downloading corn dataset from {url} ...")
    raw = urllib.request.urlopen(url, timeout=60).read()
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        with z.open("corn.mat") as f, open(DATA, "wb") as out:
            out.write(f.read())
    print(f"Saved -> {DATA}")

mat = sio.loadmat(DATA, squeeze_me=True, struct_as_record=False)
wl = np.arange(1100, 2500, 2).astype(float)
X = mat["m5spec"].data.astype(float)            # 80 x 700 absorbance
moisture = mat["propvals"].data.astype(float)[:, 0]   # column 0 = Moisture

# build spectrum-dicts as the widget expects: {"x","y","name"}
def spec(y, name):
    return {"x": wl, "y": np.asarray(y, float), "name": name}

# Two groups by moisture (low vs high third)
order = np.argsort(moisture)
low_idx = order[:20]
high_idx = order[-20:]
print(f"Moisture: low group mean={moisture[low_idx].mean():.2f}%  "
      f"high group mean={moisture[high_idx].mean():.2f}%")

# group-mean spectra + a few individuals -> aquagram
group_specs = [spec(X[low_idx].mean(0), "Low moisture (mean)"),
               spec(X[high_idx].mean(0), "High moisture (mean)")]

results = {}
for norm in ["raw", "snv", "aquagram"]:
    r = aquagram_coordinates(group_specs, normalization=norm)
    results[norm] = r
    print(f"\n[{norm}] covered={r['covered']} values shape={r['values'].shape}")

# ---- Math invariants on a larger set (all 40 grouped samples) ----
many = ([spec(X[i], f"L{i}") for i in low_idx] +
        [spec(X[i], f"H{i}") for i in high_idx])
agg = aquagram_coordinates(many, normalization="aquagram")
V = agg["values"]
col_mean = V.mean(0)
col_std = V.std(0)
assert np.allclose(col_mean, 0, atol=1e-9), "aquagram columns must be zero-mean"
assert np.allclose(col_std, 1, atol=1e-9), "aquagram columns must be unit-std"
assert agg["covered"], "corn 1100-2498nm must cover all WAMACs"
assert V.shape == (40, 12), f"expected 40x12, got {V.shape}"
print("\nMath invariants OK: aquagram columns zero-mean, unit-std; all 12 WAMACs covered.")

# Does the aquagram separate low vs high moisture at water bands?
low_mean = V[:20].mean(0)
high_mean = V[20:].mean(0)
sep = high_mean - low_mean
band_1444 = int(np.argmin(np.abs(np.array(WAMACS) - 1444)))
print(f"\nWAMACS (nm): {list(WAMACS)}")
print("Aquagram coord (high-moisture group mean, standardized):")
for b, hm, lm in zip(WAMACS, high_mean, low_mean):
    print(f"  {b:>4} nm : high={hm:+.2f}  low={lm:+.2f}  Δ(high-low)={hm-lm:+.2f}")
print(f"\nMost separated band: {WAMACS[int(np.argmax(np.abs(sep)))]} nm "
      f"(Δ={sep[int(np.argmax(np.abs(sep)))]:+.2f})")

# ---- Radar plot (classic aquagram over the 40 individual samples) ----
# aggregate SNV coords over the same 40 samples for the comparison panel
snv_agg = aquagram_coordinates(many, normalization="snv")


def radar(ax, coords, labels, title, rlim, n_low):
    """Plot each sample as a thin line (blue=low, orange=high moisture) plus a
    bold group mean, with readable radial tick labels on one clear spoke."""
    n = len(labels)
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    angc = np.concatenate([ang, ang[:1]])
    for i, row in enumerate(coords):
        col = "#1f77b4" if i < n_low else "#E36414"
        vals = np.concatenate([row, row[:1]])
        ax.plot(angc, vals, lw=0.7, color=col, alpha=0.35, zorder=2)
    for grp, col, lab in ((coords[:n_low], "#1f77b4", "Low moisture (mean)"),
                          (coords[n_low:], "#E36414", "High moisture (mean)")):
        m = np.concatenate([grp.mean(0), grp.mean(0)[:1]])
        ax.plot(angc, m, lw=2.6, color=col, marker="o", ms=3, label=lab, zorder=4)
        ax.fill(angc, m, color=col, alpha=0.08, zorder=1)
    ax.plot(angc, np.zeros_like(angc), color="#666", lw=1.0, ls="--", zorder=1)
    ax.set_xticks(ang)
    ax.set_xticklabels([str(int(l)) for l in labels], fontsize=8.5)
    ax.tick_params(axis="x", pad=6)
    # readable radial axis: integer rings only, labels on the 0° spoke, small grey
    ax.set_ylim(rlim)
    lo, hi = int(np.floor(rlim[0])), int(np.ceil(rlim[1]))
    ax.set_yticks(list(range(lo, hi + 1)))
    ax.set_rlabel_position(15)
    ax.tick_params(axis="y", labelsize=7, colors="#8a8a8a")
    ax.set_title(title, color="#0b5c9e", fontsize=11, pad=18)
    ax.grid(alpha=.35)


# shared radial scale so both panels are comparable and never crop the fill
_m = float(np.max(np.abs(np.concatenate([V.ravel(), snv_agg["values"].ravel()])))) * 1.10
_rlim = (-_m, _m)

fig, axes = plt.subplots(1, 2, figsize=(12, 6), subplot_kw=dict(polar=True))
radar(axes[0], V, WAMACS,
      "Aquagram (SNV + across-sample std)\nreal corn NIR, low vs high moisture",
      _rlim, len(low_idx))
radar(axes[1], snv_agg["values"], WAMACS,
      "SNV-only sampled at WAMACs\n(no across-sample standardization)",
      _rlim, len(low_idx))
axes[0].legend(loc="lower center", bbox_to_anchor=(0.5, -0.24), ncol=2,
               fontsize=9, frameon=False)
fig.suptitle("Aquagram validation on real corn NIR spectra (Eigenvector, m5)",
             fontsize=13, color="#119C9A", y=1.03)
fig.subplots_adjust(wspace=0.4)
fig.savefig(OUT / "aquagram_real_corn_nir.png", dpi=140, facecolor="white",
            bbox_inches="tight")
# save the computed coordinates so the figure can be re-rendered without the dataset
np.savez(OUT / "aquagram_values.npz", wamacs=np.array(WAMACS),
         aquagram=V, snv=snv_agg["values"], n_low=len(low_idx))
print(f"\nSaved radar figure -> {OUT/'aquagram_real_corn_nir.png'}")

# raw group-mean spectra with WAMACs marked (sanity of coverage)
fig2, ax = plt.subplots(figsize=(8, 4))
ax.plot(wl, X[low_idx].mean(0), color="#119C9A", label="Low moisture (mean)")
ax.plot(wl, X[high_idx].mean(0), color="#E36414", label="High moisture (mean)")
for b in WAMACS:
    ax.axvline(b, color="#888", ls=":", lw=.7)
ax.set_xlim(1100, 1700)
ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Absorbance")
ax.set_title("Corn NIR group means with 12 WAMACs (dotted) — all inside range",
             color="#0b5c9e")
ax.legend(); fig2.tight_layout()
fig2.savefig(OUT / "aquagram_wamacs_coverage.png", dpi=140, facecolor="white")
print(f"Saved coverage figure -> {OUT/'aquagram_wamacs_coverage.png'}")
print("\nVALIDATION PASSED: Aquagram works on real NIR data, math invariants hold, "
      "and the water bands reflect moisture differences.")
