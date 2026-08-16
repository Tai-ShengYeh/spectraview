"""CJK-capable font setup for the matplotlib previews.

matplotlib's default font (DejaVu Sans) has no CJK glyphs, so Chinese /
Japanese / Korean spectrum names (e.g. a file called 硃砂_5.txt) render as
boxes in the preview plots even though the Qt UI shows them fine. Calling
:func:`enable_cjk_fonts` prepends the first CJK-capable font found on the
system to matplotlib's sans-serif fallback list - Latin text keeps its
usual look, CJK glyphs stop being tofu.
"""

# Ordered by how likely they are to exist per platform: zh-TW / zh-CN /
# ja / ko Windows fonts first, then Noto (Linux), then macOS families.
CJK_CANDIDATES = [
    "Microsoft JhengHei",       # Windows, Traditional Chinese
    "Microsoft YaHei",          # Windows, Simplified Chinese
    "Yu Gothic UI", "Meiryo",   # Windows, Japanese
    "Malgun Gothic",            # Windows, Korean
    "Noto Sans CJK TC", "Noto Sans TC",
    "Noto Sans CJK JP", "Noto Sans CJK SC", "Noto Sans CJK KR",
    "PingFang TC", "PingFang SC", "Hiragino Sans",   # macOS
    "Arial Unicode MS", "SimHei", "PMingLiU",
]

_selected = None
_done = False


def enable_cjk_fonts():
    """Idempotent; returns the chosen family name, or None if the system
    has none of the known CJK fonts (plots then keep the default font)."""
    global _selected, _done
    if _done:
        return _selected
    _done = True
    try:
        from matplotlib import font_manager, rcParams
        available = {f.name for f in font_manager.fontManager.ttflist}
        _selected = next((c for c in CJK_CANDIDATES if c in available), None)
        if _selected:
            fallback = list(rcParams["font.sans-serif"])
            if _selected not in fallback:
                rcParams["font.sans-serif"] = [_selected] + fallback
        rcParams["axes.unicode_minus"] = False   # proper minus sign too
    except Exception:  # noqa: BLE001 - fonts are cosmetic, never fatal
        _selected = None
    return _selected


# Importing this module applies the setup (like matplotlib.use()).
enable_cjk_fonts()
