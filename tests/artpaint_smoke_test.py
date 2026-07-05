# -*- coding: utf-8 -*-
"""Smoke test for the Art Paint hyperspectral teaching pipeline.

Runs the same chemometrics workflow that docs/artpaint.html and
examples/artpaint_pls.ipynb teach: load DSO .mat -> log(1/R) -> per-sample
mean spectra -> PLS LOO-CV -> pixel-level prediction of cube C.

Skips (exit 0 with a SKIP message) when the dataset is not present, so CI
without the restricted .mat files still passes.  Output is ASCII only.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "ArtImageDataA", "ArtImageDataA")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print("[%s] %s%s" % (status, name, (" - " + detail) if detail else ""))
    if not cond:
        FAILURES.append(name)


def main():
    need = ["PaintCubeB.mat", "PaintCubeC.mat", "ArtImageY.mat"]
    if not all(os.path.exists(os.path.join(DATA, f)) for f in need):
        print("[SKIP] ArtImageDataA dataset not found - nothing to test")
        return 0

    import scipy.io as sio
    from sklearn.cross_decomposition import PLSRegression
    from sklearn.model_selection import LeaveOneOut

    m = sio.loadmat(os.path.join(DATA, "PaintCubeB.mat"))
    dso = m["PaintCubeB"][0, 0]
    X = dso["data"].astype(float) / 65536.0
    wl = np.ravel(dso["axisscale"][1, 0])
    mask = m["PaintMask"].astype(int)

    check("cube B shape", X.shape == (57600, 207), str(X.shape))
    check("wavelength range", abs(wl[0] - 988.9) < 0.1 and abs(wl[-1] - 1674.7) < 0.1,
          "%.1f-%.1f nm" % (wl[0], wl[-1]))
    check("mask ids 1..24", mask.min() == 1 and mask.max() == 24)
    check("reflectance in [0,1]", 0 <= X.min() and X.max() <= 1.0,
          "min=%.4f max=%.4f" % (X.min(), X.max()))

    A = np.log10(1.0 / np.clip(X, 1e-6, None))
    mask_f = mask.reshape(-1, order="F")
    S = np.array([A[mask_f == k].mean(axis=0) for k in range(1, 25)])

    ym = sio.loadmat(os.path.join(DATA, "ArtImageY.mat"))
    yTot = ym["yTotalPercent"][0, 0]["data"].astype(float)
    y24 = np.array([yTot[mask_f == k].mean(axis=0) for k in range(1, 25)])
    check("y sums to 100", np.allclose(y24.sum(axis=1), 100, atol=0.1))
    check("sample 1 Prussian ~32.4", abs(y24[0, 0] - 32.44) < 0.1,
          "%.2f" % y24[0, 0])

    yP = y24[:, 0]
    pred = np.zeros(24)
    for tr, te in LeaveOneOut().split(S):
        pls = PLSRegression(n_components=5).fit(S[tr], yP[tr])
        pred[te] = pls.predict(S[te]).ravel()
    rmsecv = float(np.sqrt(np.mean((pred - yP) ** 2)))
    check("PLS LOO RMSECV < 8%", rmsecv < 8.0, "RMSECV=%.2f" % rmsecv)

    mC = sio.loadmat(os.path.join(DATA, "PaintCubeC.mat"))
    XC = np.log10(65536.0 / np.clip(mC["PaintCubeC"][0, 0]["data"].astype(float), 1, None))
    maskC_f = mC["PaintMask"].astype(int).reshape(-1, order="F")

    rng = np.random.default_rng(0)
    idx = rng.choice(A.shape[0], 6000, replace=False)
    pls_pix = PLSRegression(n_components=5).fit(A[idx], yTot[idx, 0])
    predC = pls_pix.predict(XC).ravel()
    predC_mean = np.array([predC[maskC_f == k].mean() for k in range(1, 25)])
    rmsep = float(np.sqrt(np.mean((predC_mean - y24[:, 0]) ** 2)))
    check("cube C per-sample RMSEP < 6%", rmsep < 6.0, "RMSEP=%.2f" % rmsep)

    print("-" * 50)
    if FAILURES:
        print("FAILED: %d check(s): %s" % (len(FAILURES), ", ".join(FAILURES)))
        return 1
    print("All artpaint pipeline checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
