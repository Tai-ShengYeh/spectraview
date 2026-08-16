"""Headless test for the CJK font setup used by the matplotlib previews.

Run:  python tests/test_mplfonts.py
Skips gracefully if the machine has none of the known CJK fonts.
"""
import os
import sys
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib  # noqa: E402

matplotlib.use("Agg")

from orangespectra.mplfonts import enable_cjk_fonts  # noqa: E402

PASS = FAIL = 0


def check(name, cond):
    global PASS, FAIL
    PASS, FAIL = PASS + cond, FAIL + (not cond)
    print(f"  {'ok  ' if cond else 'FAIL'} {name}")


fam = enable_cjk_fonts()
check("enable_cjk_fonts() runs without error", True)
check("idempotent (same result twice)", enable_cjk_fonts() == fam)

from matplotlib import rcParams  # noqa: E402

check("unicode minus disabled", rcParams["axes.unicode_minus"] is False)

if fam is None:
    print("  (no CJK font on this machine - glyph test skipped)")
else:
    check("chosen family is first in fallback list",
          rcParams["font.sans-serif"][0] == fam)
    import matplotlib.pyplot as plt
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 2, 1], label="硃砂_5")
        ax.legend()
        ax.set_title("石膏 gypsum")
        fig.canvas.draw()
        plt.close(fig)
        missing = [w for w in wlist if "missing from font" in str(w.message)]
    check("CJK labels render without missing-glyph warnings", not missing)

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)
