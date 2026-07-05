# -*- coding: utf-8 -*-
"""Convert Art Paint Pigment hypercube .mat files to Orange .tab tables.

The Art Paint dataset (IASIM-10 workshop, Universitat de Barcelona) stores
NIR hyperspectral cubes as PLS_Toolbox DSO structs inside .mat files, which
Orange cannot read directly.  This script flattens a cube into an Orange
native .tab file:

  - one row per pixel, one column per wavelength (207 columns)
  - meta columns: map_x, map_y (pixel coordinates, required by the
    HyperSpectra widget of the Orange-Spectroscopy add-on) and sample (1-24)
  - if ArtImageY.mat sits next to the cube, y columns are added too and
    Prussian is marked as the learning target (class column)

Usage (run inside the ArtImageDataA folder):

  python artpaint_to_orange.py PaintCubeB.mat --stride 2 -o artpaint_B.tab
  python artpaint_to_orange.py PaintCubeC.mat --stride 2 -o artpaint_C.tab

--stride N keeps every Nth pixel in both directions (stride 2 -> 120x120 =
14400 rows, loads fast in Orange; stride 1 keeps all 57600 pixels).
"""
import argparse
import os
import sys

import numpy as np
import scipy.io as sio

Y_LABELS = ["Prussian", "Heliogen", "Ultramarine", "Oil"]


def load_cube(path):
    m = sio.loadmat(path)
    var = next(k for k in m if not k.startswith("__") and k != "PaintMask")
    dso = m[var][0, 0]
    refl = dso["data"].astype(np.float64) / 65536.0  # V = R * 65536
    wl = np.ravel(dso["axisscale"][1, 0]).astype(np.float64)
    mask = m["PaintMask"].astype(int) if "PaintMask" in m else None
    ny, nx = (int(v) for v in np.ravel(dso["imagesize"]))
    return refl, wl, mask, ny, nx, var


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mat", help="PaintCubeB.mat or PaintCubeC.mat")
    ap.add_argument("-o", "--out", default=None, help="output .tab path")
    ap.add_argument("--stride", type=int, default=2,
                    help="keep every Nth pixel (default 2)")
    args = ap.parse_args(argv)

    refl, wl, mask, ny, nx, var = load_cube(args.mat)
    out = args.out or os.path.splitext(args.mat)[0] + ".tab"

    # DSO rows follow MATLAB column-major linear indexing
    xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))  # yy=row, xx=col
    map_x = xx.reshape(-1, order="F")
    map_y = yy.reshape(-1, order="F")
    keep = (map_x % args.stride == 0) & (map_y % args.stride == 0)

    sample = (mask.reshape(-1, order="F")
              if mask is not None else np.zeros(refl.shape[0], int))

    ydir = os.path.dirname(os.path.abspath(args.mat))
    ypath = os.path.join(ydir, "ArtImageY.mat")
    yTot = None
    if os.path.exists(ypath) and refl.shape[0] == ny * nx == 57600:
        yTot = sio.loadmat(ypath)["yTotalPercent"][0, 0]["data"].astype(float)

    idx = np.where(keep)[0]
    wl_names = ["%.2f" % v for v in wl]
    metas = ["map_x", "map_y", "sample"]
    ycols = Y_LABELS if yTot is not None else []

    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(wl_names + ycols + metas) + "\n")
        f.write("\t".join(["continuous"] * (len(wl_names) + len(ycols))
                          + ["continuous", "continuous", "discrete"]) + "\n")
        flags = [""] * len(wl_names)
        if ycols:
            flags += ["class"] + ["meta"] * (len(ycols) - 1)
        flags += ["meta", "meta", "meta"]
        f.write("\t".join(flags) + "\n")
        for i in idx:
            row = ["%.5f" % v for v in refl[i]]
            if yTot is not None:
                row += ["%.3f" % v for v in yTot[i]]
            row += [str(map_x[i]), str(map_y[i]), str(sample[i])]
            f.write("\t".join(row) + "\n")

    print("wrote %s: %d pixels x %d wavelengths (from %s, stride %d)"
          % (out, len(idx), len(wl_names), var, args.stride))
    if yTot is None:
        print("note: ArtImageY.mat not found next to the cube, "
              "no concentration columns added")
    return 0


if __name__ == "__main__":
    sys.exit(main())
