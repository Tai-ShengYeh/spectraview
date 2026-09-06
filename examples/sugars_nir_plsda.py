"""Reproduce the food-analysis example of the orange-spectra paper.

Data: ``examples/sugars_nir_replicates.csv`` -- 45 NIR absorbance spectra
(9 substances x 5 replicate scans; InnoSpectra NIR-S-R14, TI DLP Hadamard
scan, 1600-2400 nm, 200 points). Column names are wavelengths in nm; the
``name`` / ``class`` columns give the substance and sugar-vs-additive label,
exactly the layout Orange's File widget reads.

What it computes (all with ``orangespectra.core``, i.e. the same functions the
widgets call):

1. leave-one-replicate-out library identification (Spectral Library):
   the library is rebuilt from the other four scans of every substance and the
   held-out scan is ranked by Pearson correlation;
2. NNLS unmixing (Mixture Analysis) of every held-out scan against that library;
3. PLS-DA sugar-vs-additive and 9-class, leave-one-out accuracy vs the number
   of components, and the VIP > 1 bands of the 3-component two-class model;
4. the two-panel figure ``paper/fig_sugars_nir.png``.

Run from the repository root:  python examples/sugars_nir_plsda.py
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "orange-spectra"))
from orangespectra import core  # noqa: E402

SUGARS = {"fructose", "glucose", "lactose", "maltose", "sucrose"}


def load_replicates(path):
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    head, body = rows[0], rows[1:]
    x = np.array([float(v) for v in head[3:]])
    X = np.array([[float(v) for v in r[3:]] for r in body])
    names = [r[0] for r in body]
    return x, X, names


def loo_library(x, X, names):
    """Leave-one-replicate-out ranked search; returns (n_correct, margins)."""
    subs = sorted(set(names))
    correct, margins, top1, nnls_frac, nnls_r2 = 0, [], [], [], []
    for i in range(len(names)):
        lib = []
        for s in subs:
            reps = [X[j] for j in range(len(names)) if names[j] == s and j != i]
            lib.append(core.make_spectrum(x, np.mean(reps, axis=0), name=s))
        q = core.make_spectrum(x, X[i], name=names[i])
        hits = core.search_library(q, lib, rank_by="correlation")
        c1, c2 = (h["scores"]["correlation"] for h in hits[:2])
        correct += hits[0]["name"] == names[i]
        margins.append(c1 - c2)
        top1.append(c1)
        r = core.mixture_nnls(q, lib)
        frac = dict(zip(r["names"], r["fractions"]))
        nnls_frac.append(frac[names[i]])
        nnls_r2.append(r["r_squared"])
    return (correct, np.array(margins), np.array(top1),
            np.array(nnls_frac), np.array(nnls_r2))


def loo_plsda(X, labels, n_components):
    ok = 0
    for i in range(len(labels)):
        m = np.arange(len(labels)) != i
        r = core.plsda_fit(X[m], [labels[j] for j in np.nonzero(m)[0]],
                           n_components)
        y_hat = (X[i] - r["x_mean"]) @ r["coef"] + r["y_mean"]
        ok += r["classes"][int(np.argmax(y_hat))] == labels[i]
    return ok


def vip_bands(x, vip, threshold=1.0):
    bands, start = [], None
    above = vip > threshold
    for i, a in enumerate(above):
        if a and start is None:
            start = i
        if start is not None and (not a or i == len(above) - 1):
            end = i if not a else i + 1
            bands.append((float(x[start]), float(x[end - 1]),
                          float(vip[start:end].max())))
            start = None
    return bands


def main():
    x, X, names = load_replicates(os.path.join(HERE, "sugars_nir_replicates.csv"))
    two = ["sugar" if n in SUGARS else "additive" for n in names]
    n = len(names)
    print(f"{n} spectra, {x.size} points, {x[0]:.0f}-{x[-1]:.0f} nm")

    correct, margins, top1, nnls_frac, nnls_r2 = loo_library(x, X, names)
    print(f"library LOO: {correct}/{n} top-1 correct; top-1 correlation >= "
          f"{top1.min():.4f}; margin to runner-up min {margins.min():.4f}, "
          f"median {np.median(margins):.4f}")
    print(f"NNLS LOO: fraction assigned to the right substance min "
          f"{100 * nnls_frac.min():.1f} %, median {100 * np.median(nnls_frac):.1f} %; "
          f"R^2 >= {nnls_r2.min():.4f}")

    for a in (1, 2, 3):
        print(f"PLS-DA sugar vs additive, {a} comp.: LOO {loo_plsda(X, two, a)}/{n}")
    for a in (4, 6, 8):
        print(f"PLS-DA 9-class, {a} comp.: LOO {loo_plsda(X, names, a)}/{n}")

    r = core.plsda_fit(X, two, 3)
    print("explained X variance per component:",
          np.round(r["explained_x_variance"], 3))
    bands = vip_bands(x, r["vip"])
    for lo, hi, mx in bands:
        sel = (x >= lo) & (x <= hi)
        s_mean = X[[i for i, t in enumerate(two) if t == "sugar"]][:, sel].mean()
        a_mean = X[[i for i, t in enumerate(two) if t == "additive"]][:, sel].mean()
        sd = np.mean([X[[i for i, nm in enumerate(names) if nm == s]][:, sel]
                      .std(axis=0).mean() for s in set(names)])
        print(f"VIP>1 band {lo:.0f}-{hi:.0f} nm: max VIP {mx:.2f}; mean A "
              f"sugars {s_mean:.3f} vs additives {a_mean:.3f}; replicate SD {sd:.4f}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; figure skipped")
        return
    fig, ax = plt.subplots(2, 1, figsize=(6.5, 6.2), sharex=True)
    for s in sorted(set(names)):
        ax[0].plot(x, X[[i for i, nm in enumerate(names) if nm == s]].mean(axis=0),
                   lw=1.1, label=s, ls="-" if s in SUGARS else "--")
    ax[0].set_ylabel("absorbance")
    ax[0].legend(fontsize=7, ncol=3, frameon=False)
    ax[0].set_title("(a) five sugars (solid) and four additives (dashed), "
                    "mean of 5 replicates", fontsize=9)
    ax[1].plot(x, r["vip"], color="#c44e52", lw=1.2)
    ax[1].axhline(1, color="#888", ls=":", lw=0.8)
    ax[1].set_ylabel("VIP")
    ax[1].set_xlabel("wavelength (nm)")
    ax[1].set_title("(b) PLS-DA sugar vs additive (3 components): VIP scores",
                    fontsize=9)
    for lo, hi, _ in bands:
        if hi - lo >= 20:
            ax[1].axvspan(lo, hi, color="#c44e52", alpha=0.08)
    plt.tight_layout()
    out = os.path.join(os.path.dirname(HERE), "paper", "fig_sugars_nir.png")
    plt.savefig(out, dpi=200)
    print("figure ->", out)


if __name__ == "__main__":
    main()
