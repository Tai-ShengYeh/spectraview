# -*- coding: utf-8 -*-
"""Smoke test for docs/artpaint.html and its companion example files.

Checks the teaching page is self-contained and well-formed:
  - embedded data block parses as JSON with the expected keys/shapes
  - all 8 sections and the 10-question quiz are present
  - companion files (notebook, Orange converter) exist and are valid
  - the Orange converter produces a well-formed 3-header .tab (runs only
    when the dataset is present; loads it into Orange when installed)

ASCII output only.
"""
import io
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
PAGE = os.path.join(ROOT, "docs", "artpaint.html")
DATA = os.path.join(ROOT, "ArtImageDataA", "ArtImageDataA")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print("[%s] %s%s" % (status, name, (" - " + detail) if detail else ""))
    if not cond:
        FAILURES.append(name)


def main():
    check("artpaint.html exists", os.path.exists(PAGE))
    if not os.path.exists(PAGE):
        return finish()

    html = io.open(PAGE, encoding="utf-8").read()

    m = re.search(r"const D = (\{.*?\});\n", html, re.S)
    check("embedded data block found", m is not None)
    if m:
        d = json.loads(m.group(1))
        check("207 wavelengths", len(d.get("wl", [])) == 207)
        check("24 spectra x 207", len(d.get("spectraB", [])) == 24
              and all(len(r) == 207 for r in d["spectraB"]))
        check("24 y rows x 4", len(d.get("y", [])) == 24
              and all(len(r) == 4 for r in d["y"]))
        check("pca scores present", len(d.get("pca", {}).get("pc1", [])) == 24)
        check("pls cv results present", len(d.get("pls", {}).get("pred", [])) == 24)
        # maskImgB is tiny (24 flat-color blocks compress extremely well)
        check("three base64 images", all(
            len(d.get(k, "")) > 500 for k in ("meanImgB", "maskImgB", "predMapC")))
        check("plausible RMSECV", 0 < d["pls"]["rmsecv"] < 8,
              str(d["pls"]["rmsecv"]))

    for sec in ["sec-data", "sec-explore", "sec-chemo", "sec-python",
                "sec-colab", "sec-r", "sec-orange", "sec-quiz"]:
        check("section " + sec, ('id="%s"' % sec) in html)

    quiz_qs = len(re.findall(r"\{cat:\"(?:chem|tool)\"", html))
    check("10 quiz questions", quiz_qs == 10, "found %d" % quiz_qs)
    check("quiz grading js", "btnGrade" in html and "btnReset" in html)
    check("HyperSpectra mentioned", "HyperSpectra" in html)
    check("dataset credit present", "IASIM-10" in html and "Garc" in html)
    check("no CDN scripts", "<script src=" not in html)

    nb_path = os.path.join(ROOT, "examples", "artpaint_pls.ipynb")
    check("notebook exists", os.path.exists(nb_path))
    if os.path.exists(nb_path):
        nb = json.load(io.open(nb_path, encoding="utf-8"))
        check("notebook has cells", len(nb.get("cells", [])) >= 10,
              "%d cells" % len(nb.get("cells", [])))

    conv = os.path.join(ROOT, "examples", "artpaint_to_orange.py")
    check("orange converter exists", os.path.exists(conv))

    cube = os.path.join(DATA, "PaintCubeB.mat")
    if os.path.exists(conv) and os.path.exists(cube):
        import tempfile
        out = os.path.join(tempfile.gettempdir(), "artpaint_tab_smoke.tab")
        r = subprocess.run([sys.executable, conv, cube, "--stride", "8",
                            "-o", out], capture_output=True, text=True)
        check("converter runs", r.returncode == 0, r.stderr.strip()[:120])
        if r.returncode == 0:
            lines = io.open(out, encoding="utf-8").read().splitlines()
            ncol = len(lines[0].split("\t"))
            check("tab header 207+4+3 cols", ncol == 214, "%d cols" % ncol)
            check("tab has 3 header rows + data",
                  len(lines) == 3 + 30 * 30, "%d lines" % len(lines))
            check("map_x/map_y meta flagged",
                  lines[2].split("\t")[-3:] == ["meta", "meta", "meta"])
            try:
                import Orange
                t = Orange.data.Table(out)
                check("Orange loads tab", len(t) == 900
                      and len(t.domain.attributes) == 207,
                      "rows=%d" % len(t))
            except ImportError:
                print("[SKIP] Orange not installed - tab load check skipped")
        if os.path.exists(out):
            os.remove(out)
    else:
        print("[SKIP] dataset not found - converter run skipped")

    return finish()


def finish():
    print("-" * 50)
    if FAILURES:
        print("FAILED: %d check(s): %s" % (len(FAILURES), ", ".join(FAILURES)))
        return 1
    print("All artpaint page checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
