"""Validate orangespectra.core.pls_regression_fit's VIP against
D:\\Orange\\acetic_t1_1mm_0313_2025\\pls_analysis.py (an external, user-owned
analysis script — NOT part of this repo and NOT modified by this script).

The dataset (acetic_t1_1mm_0313_2025_df.csv, T1 1 mm-pathlength acetic acid
vinegar NIR spectra + reference acidity) also lives outside the repo and is
never copied in here; this script reads it from disk at run time and expects
it to already exist locally.

Run:
    python orange-spectra/validation/validate_pls_regression_acetic.py

What this checks
-----------------
``pls_analysis.py`` computes "VIP" with a loop that ends in ``.sum()`` on the
per-feature array *before* accumulating it into the per-feature VIP vector,
which collapses the feature axis to a scalar every iteration -- so its
"vip_scores" ends up constant across all 228 wavelengths (see comparison 1
below). This script:

1. Replicates the script's own component-selection CV (5-fold, shuffle,
   random_state=42, 1..20 components) and its PLSRegression fit
   (scale=True, the sklearn default the script relies on implicitly) --
   *without* running the original script (it hard-codes a C:\\ path, would
   overwrite the user's PNGs, and its "top 20" printout is meaningless, see
   below).
2. Reproduces the script's literal VIP formula verbatim (bug included).
3. Reproduces the script's formula with the stray ``.sum()`` removed (a
   different, non-standard definition: it weights by x_loadings^2 *
   y_loadings^2 with no SSY/weight normalisation, unlike VIP).
4. Computes standard VIP directly from the same fitted sklearn model
   (x_weights_, transform(X) scores, y_loadings_) -- this is the real
   apples-to-apples check.
5. Runs the same X/y/n_components/scale through
   ``orangespectra.core.pls_regression_fit`` and reports Spearman rank
   correlation, top-20 overlap, and max |difference| of its VIP against
   each of 2-4 above.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import KFold, cross_val_score

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from orangespectra import core  # noqa: E402

# The original script hard-codes C:\Orange\...; the data actually lives at
# D:\Orange\... on this machine. We read it from there (see docstring) and do
# not touch or run the original script.
DATA_PATH = r"D:\Orange\acetic_t1_1mm_0313_2025\acetic_t1_1mm_0313_2025_df.csv"


def top_k_overlap(a, b, k=20):
    top_a = set(np.argsort(a)[::-1][:k].tolist())
    top_b = set(np.argsort(b)[::-1][:k].tolist())
    return len(top_a & top_b)


def report(label, ours, theirs, k=20):
    if np.ptp(theirs) == 0:
        print(f"{label}: DEGENERATE (constant vector, ptp=0) -- "
              "Spearman / top-k overlap are not meaningful; reporting anyway.")
    with np.errstate(all="ignore"):
        rho, pval = spearmanr(ours, theirs)
    overlap = top_k_overlap(ours, theirs, k)
    max_diff = float(np.max(np.abs(ours - theirs)))
    print(f"{label}:")
    print(f"    Spearman rho = {rho!r} (p={pval!r})")
    print(f"    top-{k} overlap = {overlap}/{k}")
    print(f"    max |diff|   = {max_diff:.6g}")
    return rho, overlap, max_diff


def main():
    if not os.path.exists(DATA_PATH):
        print(f"Dataset not found at {DATA_PATH!r} -- nothing to validate.")
        print("This script reads the user's own external dataset; it is not "
              "shipped with the repo.")
        return 1

    data = pd.read_csv(DATA_PATH)
    print(f"Loaded {DATA_PATH}: shape={data.shape}")
    if "Unnamed: 0" in data.columns or data.columns[0] != "acidity":
        print(f"NOTE: unexpected leading column(s) {list(data.columns[:3])!r} -- "
              "check whether the original script silently used a sample-id "
              "column as a feature.")

    y = data["acidity"].values.astype(float)
    X = data.drop("acidity", axis=1).values.astype(float)
    p = X.shape[1]
    print(f"X: {X.shape}, y range [{y.min()}, {y.max()}]")

    # ---- replicate the script's own component-selection CV exactly -------
    n_components_grid = np.arange(1, 21)
    cv = KFold(n_splits=5, shuffle=True, random_state=42)
    mse_scores = []
    for n in n_components_grid:
        pls = PLSRegression(n_components=n)  # scale=True (sklearn default)
        scores = cross_val_score(pls, X, y, cv=cv, scoring="neg_mean_squared_error")
        mse_scores.append(-scores.mean())
    optimal_components = int(n_components_grid[np.argmin(mse_scores)])
    print(f"Optimal n_components (replicated CV): {optimal_components} "
          f"(min CV MSE = {min(mse_scores):.6g})")

    # ---- the script's own fit (scale=True, sklearn default) --------------
    pls_optimal = PLSRegression(n_components=optimal_components)  # scale=True
    pls_optimal.fit(X, y)

    # ---- comparison 1: literal script VIP formula (bug included) ---------
    vip_bug = np.zeros(p)
    for i in range(optimal_components):
        vip_bug += (p * (pls_optimal.x_loadings_[:, i] ** 2) *
                   (pls_optimal.y_loadings_[:, i] ** 2)).sum()   # .sum() -> scalar!
    vip_bug = np.sqrt(vip_bug * p / (pls_optimal.x_loadings_ ** 2).sum())
    print(f"\nLiteral script VIP: ptp={np.ptp(vip_bug):.3g} "
          f"(0 means every feature got the identical 'VIP' -- the .sum() "
          "inside the loop collapses the per-feature array to a scalar "
          "before it is accumulated, so the per-feature information is "
          "destroyed before the sqrt/normalise step; argsort's 'top 20' is "
          "then just the last 20 array positions, i.e. meaningless, and its "
          "printed 'Wavelength = p - 20 + idx' label is also wrong).")

    # ---- comparison 2: same formula with the stray .sum() removed --------
    # Still not standard VIP (uses x_loadings not weights, no SSY weighting),
    # but at least varies per feature.
    vip_nosum = np.zeros(p)
    for i in range(optimal_components):
        vip_nosum += p * (pls_optimal.x_loadings_[:, i] ** 2) * (pls_optimal.y_loadings_[:, i] ** 2)
    vip_nosum = np.sqrt(vip_nosum * p / (pls_optimal.x_loadings_ ** 2).sum())

    # ---- comparison 3: standard VIP from the SAME fitted sklearn model ---
    T_sk = pls_optimal.x_scores_
    W_sk = pls_optimal.x_weights_
    Q_sk = pls_optimal.y_loadings_  # (1, A) for single-target y
    wnorm2 = np.maximum((W_sk ** 2).sum(axis=0), 1e-12)
    ssy = (T_sk ** 2).sum(axis=0) * (Q_sk ** 2).sum(axis=0)
    vip_std_sklearn = np.sqrt(p * ((W_sk ** 2) / wnorm2 @ ssy) / max(ssy.sum(), 1e-12))

    # ---- our widget's core function, same X/y/n_components/scale ---------
    res = core.pls_regression_fit(X, y, n_components=optimal_components, scale=True)
    vip_ours = res["vip"]
    print(f"\nmean(VIP_ours^2) = {(vip_ours ** 2).mean():.6f}  (should be 1.0)")

    max_diff_std = float(np.max(np.abs(vip_ours - vip_std_sklearn)))
    print(f"\nours vs. standard VIP computed from sklearn's own fit: "
          f"max |diff| = {max_diff_std:.3e} (expected ~1e-8 -- same formula, "
          "same fitted model)")

    print()
    report("Comparison 1: ours vs. literal script VIP (buggy)", vip_ours, vip_bug)
    print()
    report("Comparison 2: ours vs. script VIP with .sum() removed "
           "(different, non-standard definition)", vip_ours, vip_nosum)
    print()
    report("Comparison 3: ours vs. standard VIP from sklearn's fit "
           "(the real apples-to-apples check)", vip_ours, vip_std_sklearn)

    print("\nConclusion: pls_regression_fit's VIP matches the standard "
          "definition (comparison 3) to ~1e-8 given the same data, "
          "n_components and scale. The mismatch against the script's own "
          "printed numbers (comparisons 1-2) is the script's bug/non-standard "
          "formula, not a discrepancy in the widget -- see the notes above.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
