"""Pure spectral logic for the Orange widgets (no Orange / Qt imports here).

Everything in this module works on plain numpy arrays and simple dicts, so it
can be unit-tested headlessly. The logic mirrors SpectraView's proven modules:

  * URL import  — IRUG jqPlot pages, SOPRANO Dygraph pages, JCAMP-DX (AFFN),
                  and two-column CSV/TSV, matching ``specview.formats``.
  * similarity  — correlation / cosine / spectral angle / Euclidean, matching
                  ``specview.library.similarity_scores``.
  * library     — JSON persistence **compatible with SpectraView .speclib**.
  * mixture     — non-negative least squares, matching
                  ``specview.analysis.mixture_nnls``.

A "spectrum" here is a dict: {name, x, y, x_label, source}.
"""
from __future__ import annotations

import ast
import html as _html
import json
import os
import re
import urllib.parse
import urllib.request

import numpy as np

_EPS = 1e-12
USER_AGENT = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/124.0 Safari/537.36"
)
IRUG_DETAIL_URL = "http://www.irug.org/jcamp-details?id={id}"


def make_spectrum(x, y, name="spectrum", x_label="x", source=""):
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    n = min(x.size, y.size)
    x, y = x[:n], y[:n]
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    order = np.argsort(x, kind="stable")
    x, y = x[order], y[order]
    if x.size:
        uniq = np.concatenate(([True], np.diff(x) != 0))
        x, y = x[uniq], y[uniq]
    return {"name": name, "x": x, "y": y, "x_label": x_label, "source": source}


# ===================================================================== fetch
def resolve_source(id_or_url) -> str:
    """A bare number (or ``irug:N``) becomes the IRUG detail-page URL."""
    s = str(id_or_url).strip()
    if not s:
        raise ValueError("Empty source. Enter an IRUG id or a URL.")
    m = re.fullmatch(r"(?:irug:)?(\d+)", s, re.IGNORECASE)
    if m:
        return IRUG_DETAIL_URL.format(id=m.group(1))
    if not re.match(r"https?://", s, re.IGNORECASE):
        raise ValueError(f"Not an IRUG id or URL: {s!r}")
    return s


def default_fetch(url: str) -> tuple[bytes, str]:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=60) as resp:  # noqa: S310
        return resp.read(), resp.headers.get("Content-Type", "") or ""


# --------------------------------------------------------------- IRUG jqPlot
def parse_irug_jqplot(text: str, source: str = ""):
    """Parse the spectrum IRUG embeds in its interactive (jqPlot) viewer.

    Each spectrum is a JSON object assigned to ``jqPlotData.series['Name']``::

        jqPlotData.series['Submitter'] = {"1900.0":0.384, "1899.1":0.413, ...}

    i.e. the **wavenumber is the quoted key** and the **intensity the value**.
    A page may hold several series (submitted sample + reference matches); we
    take the submitted one when labelled, otherwise the first.
    """
    series = re.findall(
        r"jqPlotData\.series\s*\[\s*['\"]([^'\"]*)['\"]\s*\]\s*=\s*\{([^{}]*)\}",
        text, re.IGNORECASE)
    if series:
        block = next((b for name, b in series if "submit" in name.lower()),
                     series[0][1])
    else:
        scripts = re.findall(r"<script\b[^>]*>(.*?)</script>", text,
                             re.DOTALL | re.IGNORECASE)
        block = next((s for s in scripts if "jqplotdata" in s.lower()), None)
        if block is None:
            block = text if "jqplotdata" in text.lower() else None
        if block is None:
            return None

    num = r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"
    pairs = re.findall(rf'"({num})"\s*:\s*({num})', block)   # real IRUG format
    if len(pairs) < 2:                                       # fallback: bare n:n
        pairs = re.findall(rf'({num})\s*:\s*({num})', block)
    if len(pairs) < 2:
        return None
    m = re.search(r"id=(\d+)", source)
    name = f"IRUG {m.group(1)}" if m else "IRUG spectrum"
    return make_spectrum([float(a) for a, _ in pairs], [float(b) for _, b in pairs],
                         name=name, x_label="wavenumber (cm-1)", source=source)


# ------------------------------------------------------------ SOPRANO Dygraph
def _find_balanced(text: str, start: int, opener="[", closer="]") -> str:
    depth, in_str, escape = 0, None, False
    for i in range(start, len(text)):
        ch = text[i]
        if in_str:
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == in_str:
                in_str = None
            continue
        if ch in ("'", '"'):
            in_str = ch
        elif ch == opener:
            depth += 1
        elif ch == closer:
            depth -= 1
            if depth == 0:
                return text[start:i + 1]
    raise ValueError("Unbalanced SOPRANO data array.")


def _strip_tags(text: str) -> str:
    return _html.unescape(re.sub(r"<[^>]+>", "", text)).strip()


def parse_soprano(text: str, source: str = ""):
    """SOPRANO pages embed the spectrum as a Dygraph [[x, y], ...] array."""
    m = re.search(r"new\s+Dygraph\s*\(", text, re.IGNORECASE)
    if not m:
        return None
    start = text.find("[[", m.end())
    if start < 0:
        return None
    arr = np.asarray(ast.literal_eval(_find_balanced(text, start)), dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 2 or arr.shape[0] < 2:
        return None
    name = ""
    lm = re.search(r"labels\s*:\s*(\[[^\]]+\])", text, re.DOTALL)
    if lm:
        try:
            labels = [_strip_tags(str(v)) for v in ast.literal_eval(lm.group(1))]
            if len(labels) > 1:
                name = labels[1]
        except Exception:
            pass
    if not name:
        qs = urllib.parse.parse_qs(urllib.parse.urlparse(source).query)
        name = qs.get("id", ["SOPRANO spectrum"])[0]
    x_label = ("Raman shift (cm-1)" if "raman" in text.lower()
               else "wavenumber (cm-1)" if "cm" in text.lower() else "x")
    return make_spectrum(arr[:, 0], arr[:, 1], name=name,
                         x_label=x_label, source=source)


# ------------------------------------------------------------ JCAMP-DX (AFFN)
def parse_jcamp(text: str, source: str = ""):
    """Plain-AFFN JCAMP-DX: ##XYPOINTS/##PEAK TABLE pairs, or ##XYDATA with
    numeric lines (compressed SQZ/DIF files need SpectraView's full reader)."""
    if "##TITLE" not in text.upper():
        return None
    ldr = {}
    data_lines, capturing, kind = [], False, None
    for raw in text.splitlines():
        line = raw.split("$$", 1)[0].rstrip()
        if line.strip().startswith("##"):
            mm = re.match(r"##\s*([^=]+)=(.*)", line.strip())
            if not mm:
                continue
            key = mm.group(1).strip().upper()
            if key in ("XYPOINTS", "PEAK TABLE", "PEAKTABLE", "XYDATA"):
                capturing, kind = True, ("xydata" if key == "XYDATA" else "pairs")
                continue
            capturing = False
            ldr[key] = mm.group(2).strip()
        elif capturing and line.strip():
            data_lines.append(line)
    if not data_lines:
        return None

    def num(key, default):
        try:
            return float(ldr.get(key, default))
        except (TypeError, ValueError):
            return default

    xf, yf = num("XFACTOR", 1.0), num("YFACTOR", 1.0)
    name = ldr.get("TITLE", "") or "JCAMP spectrum"
    nums = [[float(t) for t in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", ln)]
            for ln in data_lines]
    if kind == "pairs":
        flat = [v for row in nums for v in row]
        x = np.asarray(flat[0::2]) * xf
        y = np.asarray(flat[1::2]) * yf
    else:  # XYDATA (X++(Y..Y)), AFFN only: first number per line is the X checkpoint
        ys = [v for row in nums for v in row[1:]]
        y = np.asarray(ys) * yf
        firstx, lastx = num("FIRSTX", 0.0), num("LASTX", max(len(ys) - 1, 1))
        x = np.linspace(firstx, lastx, len(ys)) if len(ys) > 1 else np.asarray([firstx])
    xunit = ldr.get("XUNITS", "").upper()
    x_label = "wavenumber (cm-1)" if "CM" in xunit else (
        "wavelength (nm)" if "NANO" in xunit else "x")
    return make_spectrum(x, y, name=name.strip(), x_label=x_label, source=source)


# ----------------------------------------------------------------- CSV / TSV
_DELIM_RE = re.compile(r"[,\t;|]+")
_WS_RE = re.compile(r"\s+")
_COMMENT_RE = re.compile(r"^\s*[#%!']")


def split_fields(line: str) -> list:
    """Split one data line into fields.

    Uses comma / tab / semicolon / pipe when any is present; otherwise falls
    back to runs of whitespace, so space-delimited instrument exports (Raman
    ``.txt``, ``.dat``, ``.asc``, ``.xy`` …) parse like CSV. Returns ``[]`` for
    blank and comment lines (``#``, ``%``, ``!``, ``'``).
    """
    s = line.strip()
    if not s or _COMMENT_RE.match(s):
        return []
    if _DELIM_RE.search(s):
        return [p.strip() for p in _DELIM_RE.split(s) if p.strip() != ""]
    return _WS_RE.split(s)


def parse_csv(text: str, source: str = ""):
    rows = []
    header = None
    for line in text.splitlines():
        parts = split_fields(line)
        if len(parts) < 2:
            continue
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except ValueError:
            if header is None:
                header = parts
    if len(rows) < 2:
        return None
    x_label = header[0] if header else "x"
    name = os.path.splitext(os.path.basename(
        urllib.parse.urlparse(source).path))[0] or "CSV spectrum"
    return make_spectrum([r[0] for r in rows], [r[1] for r in rows],
                         name=name, x_label=x_label, source=source)


# ----------------------------------------------- generic embedded JS charts
_NUM = r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"


def _array_after(text, key_pattern):
    """Return the JS array literal (as a string) that follows ``key_pattern``
    (e.g. ``x:``), or None. The '[' must immediately follow the key."""
    for m in re.finditer(key_pattern, text, re.IGNORECASE):
        b = text.find("[", m.end())
        if b < 0 or (b - m.end()) > 3:
            continue
        try:
            return _find_balanced(text, b)
        except ValueError:
            continue
    return None


def _mk_from_pairs(pairs, source, name=None):
    xy = sorted(((float(a), float(b)) for a, b in pairs), key=lambda p: p[0])
    if name is None:
        seg = urllib.parse.urlparse(source).path.rstrip("/").split("/")[-1]
        name = seg or "spectrum"
    return make_spectrum([p[0] for p in xy], [p[1] for p in xy],
                         name=name, x_label="x", source=source)


def parse_embedded_arrays(text: str, source: str = ""):
    """Last-resort parser for spectra embedded in an interactive JS chart:
    Plotly (``x:[…], y:[…]``), Highcharts / ECharts (``data:[[x,y],…]``),
    Chart.js (``labels:[…]`` + ``data:[…]``), or any ``[[x,y],…]`` array.

    Runs only after the specific parsers, and needs several numeric points, so
    it rarely mis-fires. Helps databases that only render a chart (no file)."""
    # (1) paired arrays: data:[[x,y],…] / series data (Highcharts, ECharts, Plotly)
    for key in (r"['\"]?data['\"]?\s*:", r"['\"]?series['\"]?\s*:"):
        s = _array_after(text, key)
        if s and s.lstrip().startswith("[["):
            pairs = re.findall(rf"\[\s*({_NUM})\s*,\s*({_NUM})", s)
            if len(pairs) >= 5:
                return _mk_from_pairs(pairs, source)
    # any [[x,y],…] block anywhere (longest wins)
    best = None
    for m in re.finditer(r"\[\s*\[", text):
        try:
            s = _find_balanced(text, m.start())
        except ValueError:
            continue
        pairs = re.findall(rf"\[\s*({_NUM})\s*,\s*({_NUM})", s)
        if len(pairs) >= 8 and (best is None or len(pairs) > len(best)):
            best = pairs
    if best:
        return _mk_from_pairs(best, source)
    # (2) separate x / y arrays (Plotly-style)
    xs, ys = _array_after(text, r"['\"]?x['\"]?\s*:"), _array_after(
        text, r"['\"]?y['\"]?\s*:")
    if xs and ys:
        xv = [float(v) for v in re.findall(_NUM, xs)]
        yv = [float(v) for v in re.findall(_NUM, ys)]
        if len(xv) >= 5 and len(xv) == len(yv):
            return _mk_from_pairs(zip(xv, yv), source)
    # (3) Chart.js: labels (x) + data (y)
    lab, dat = _array_after(text, r"labels\s*:"), _array_after(text, r"data\s*:")
    if lab and dat:
        xv = [float(v) for v in re.findall(_NUM, lab)]
        yv = [float(v) for v in re.findall(_NUM, dat)]
        if len(xv) >= 5 and len(xv) == len(yv):
            return _mk_from_pairs(zip(xv, yv), source)
    return None


def load_spectrum_url(id_or_url, fetch=default_fetch):
    """Fetch + auto-detect: IRUG page, SOPRANO page, JCAMP-DX, CSV, or a
    spectrum embedded in an interactive JS chart (Plotly/Highcharts/Chart.js)."""
    url = resolve_source(id_or_url)
    raw, _ctype = fetch(url)
    text = raw.decode("utf-8", errors="replace")
    for parser in (parse_irug_jqplot, parse_soprano, parse_jcamp, parse_csv,
                   parse_embedded_arrays):
        spec = parser(text, source=url)
        if spec is not None and spec["x"].size >= 2:
            return spec
    raise ValueError(
        f"Could not extract a spectrum from {url}\n"
        "Supported: IRUG detail pages, SOPRANO pages, JCAMP-DX (AFFN), CSV/TSV, "
        "and spectra embedded in Plotly / Highcharts / Chart.js charts."
    )


# ================================================================ similarity
def _common_grid(xa, ya, xb, yb):
    lo, hi = max(xa[0], xb[0]), min(xa[-1], xb[-1])
    if hi <= lo:
        return None
    gx = xa[(xa >= lo) & (xa <= hi)]
    if gx.size < 3:
        return None
    return np.interp(gx, xa, ya), np.interp(gx, xb, yb)


def similarity_scores(xa, ya, xb, yb) -> dict:
    """correlation/cosine (higher = more similar), sam/euclid (lower = better).
    Same definitions as SpectraView's library search."""
    pair = _common_grid(np.asarray(xa, float), np.asarray(ya, float),
                        np.asarray(xb, float), np.asarray(yb, float))
    if pair is None:
        return {"correlation": 0.0, "cosine": 0.0,
                "sam": float("inf"), "euclid": float("inf")}
    va, vb = pair
    da, db = va - va.mean(), vb - vb.mean()
    denom = np.sqrt((da @ da) * (db @ db))
    corr = float(da @ db / denom) if denom > _EPS else 0.0
    na, nb = np.linalg.norm(va), np.linalg.norm(vb)
    cos = float(va @ vb / (na * nb)) if na > _EPS and nb > _EPS else 0.0
    sam = float(np.arccos(np.clip(cos, -1.0, 1.0)))
    ua = va / na if na > _EPS else va
    ub = vb / nb if nb > _EPS else vb
    return {"correlation": corr, "cosine": cos, "sam": sam,
            "euclid": float(np.linalg.norm(ua - ub))}


# ================================================================== library
def save_library(entries: list, path: str, name: str = "library") -> None:
    """Write a SpectraView-compatible .speclib (JSON) file."""
    def _entry(s):
        e = {"name": s["name"], "x_unit": _unit_from_label(s.get("x_label", "")),
             "y_unit": "intensity", "x": np.asarray(s["x"]).tolist(),
             "y": np.asarray(s["y"]).tolist()}
        for extra in ("color", "source"):
            if s.get(extra):
                e[extra] = s[extra]
        return e

    obj = {"name": name, "library": [_entry(s) for s in entries]}
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=1)


def load_library(path: str) -> list:
    """Read a SpectraView .speclib (JSON) file into spectrum dicts."""
    with open(path, "r", encoding="utf-8") as fh:
        obj = json.load(fh)
    out = []
    for e in obj.get("library", []):
        s = make_spectrum(
            e["x"], e["y"], name=e.get("name", "ref"),
            x_label=_label_from_unit(e.get("x_unit", "")),
            source=e.get("source", path))
        if e.get("color"):
            s["color"] = e["color"]
        out.append(s)
    return out


def builtin_libraries() -> dict:
    """Map display name -> path for the .speclib files shipped with the
    package (``orangespectra/libraries/``). The display name is the
    library's own "name" field, falling back to the file name."""
    lib_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "libraries")
    out = {}
    if not os.path.isdir(lib_dir):
        return out
    for fn in sorted(os.listdir(lib_dir)):
        if not fn.lower().endswith(".speclib"):
            continue
        path = os.path.join(lib_dir, fn)
        name = os.path.splitext(fn)[0]
        try:
            with open(path, "r", encoding="utf-8") as fh:
                name = json.load(fh).get("name") or name
        except Exception:  # noqa: BLE001 - unreadable file: keep file name
            pass
        out[name] = path
    return out


def _unit_from_label(label: str) -> str:
    low = label.lower()
    if "raman" in low:
        return "raman_cm-1"
    if "cm" in low:
        return "cm-1"
    if "nm" in low:
        return "nm"
    return "pixel"


def _label_from_unit(unit: str) -> str:
    return {"cm-1": "wavenumber (cm-1)", "raman_cm-1": "Raman shift (cm-1)",
            "nm": "wavelength (nm)"}.get(unit, unit or "x")


def search_library(query: dict, entries: list, rank_by: str = "correlation"):
    """Rank library entries by similarity to a query spectrum, best first."""
    hits = []
    for e in entries:
        scores = similarity_scores(query["x"], query["y"], e["x"], e["y"])
        hits.append({"entry": e, "name": e["name"], "scores": scores})
    reverse = rank_by in ("correlation", "cosine")
    hits.sort(key=lambda h: h["scores"].get(rank_by, 0.0), reverse=reverse)
    return hits



# =================================================== UCL pigment library
# The UCL Raman pigment library (Bell, Clark & Gibbs 1997) is downloaded from
# UCL's own website on first use and cached locally. The data itself is NOT
# redistributed with this package: UCL's site states its material may not be
# reproduced without permission, so each user downloads their own copy,
# exactly as they would in a browser.
UCL_BASE_URL = "https://www.chem.ucl.ac.uk/resources/raman"
# UCL's server blocks some regions/networks (403). The Internet Archive's
# copy of the same pages is world-readable, so it serves as a fallback.
UCL_WAYBACK_BASE = ("https://web.archive.org/web/2023/"
                    "https://www.chem.ucl.ac.uk/resources/raman")
UCL_LIBRARY_NAME = "UCL Raman Library of Pigments (Bell, Clark & Gibbs 1997)"
UCL_CITATION = ("I. M. Bell, R. J. H. Clark and P. J. Gibbs, Spectrochim. "
                "Acta A 53 (1997) 2159-2179, doi:10.1016/S1386-1425(97)00140-6")

# Canonical pigment names -> colour group, as organised on the UCL site.
UCL_PIGMENT_COLORS = {
    "Ivory Black": "black", "Lamp Black": "black",
    "Azurite": "blue", "Cerulean Blue": "blue", "Cobalt Blue": "blue",
    "Egyptian Blue": "blue", "Lazurite": "blue", "Prussian Blue": "blue",
    "Smalt": "blue",
    "Atacamite": "green", "Chromium(III) Oxide": "green",
    "Cobalt Green": "green", "Emerald Green": "green", "Malachite": "green",
    "Scheele's Green": "green", "Terre Verte": "green",
    "Verdigris (1)": "green", "Verdigris (2)": "green",
    "Verdigris (raw)": "green", "Viridian": "green",
    "Mars Orange": "orange",
    "Mars Red": "red", "Purpurin": "red", "Realgar": "red", "Red Lead": "red",
    "Red Ochre": "red", "Vermilion": "red",
    "Barium White": "white", "Bone White": "white", "Chalk": "white",
    "Gypsum": "white", "Lead White": "white", "Lithopone": "white",
    "Zinc White": "white",
    "Barium Yellow": "yellow", "Berberine": "yellow",
    "Cadmium Yellow": "yellow", "Chrome Yellow": "yellow",
    "Chrome Yellow (deep)": "yellow", "Chrome Yellow-Orange": "yellow",
    "Cobalt Yellow": "yellow", "Gamboge": "yellow", "Indian Yellow": "yellow",
    "Lead Tin Yellow (type I)": "yellow",
    "Lead Tin Yellow (type II)": "yellow", "Litharge": "yellow",
    "Mars Yellow": "yellow", "Massicot": "yellow", "Naples Yellow": "yellow",
    "Orpiment": "yellow", "Pararealgar": "yellow", "Saffron": "yellow",
    "Strontium Yellow": "yellow", "Yellow Ochre": "yellow",
    "Zinc Yellow": "yellow",
}

# Known lowercase file/page stems on the UCL site whose canonical name cannot
# be recovered from the file name alone. Anything not listed here is named
# from the pigment page's <title>.
UCL_STEM_NAMES = {
    "atacamit": "Atacamite", "barytes": "Barium White",
    "berberin": "Berberine", "bonewhit": "Bone White",
    "cerulblu": "Cerulean Blue", "chromoxi": "Chromium(III) Oxide",
    "cobaltgr": "Cobalt Green", "emeraldg": "Emerald Green",
    "ivoryblk": "Ivory Black", "leadwhit": "Lead White",
    "lithopon": "Lithopone", "malachit": "Malachite",
    "naplesyl": "Naples Yellow", "parareal": "Pararealgar",
    "prusblue": "Prussian Blue", "redlead": "Red Lead",
    "verdiraw": "Verdigris (raw)", "vermilio": "Vermilion",
    "zincwhit": "Zinc White",
}
_UCL_NAMES_LOWER = {n.lower(): n for n in UCL_PIGMENT_COLORS}


def read_spc(data: bytes, name: str = "spectrum", source: str = "") -> dict:
    """Parse a Thermo GRAMS .spc file (new format, LSB) into a spectrum dict.

    Supports the single-subfile files used by spectral libraries such as
    UCL's: evenly spaced x (or an explicit x array, TXVALS) with y stored as
    float32 or as int32 with a binary exponent. Multifile and the pre-1996
    "old format" are rejected with a clear error.
    """
    import struct

    if len(data) < 512 + 32:
        raise ValueError("Not an SPC file (too short).")
    ftflg, fversn, _fexper, fexp = struct.unpack_from("<BBBb", data, 0)
    if fversn == 0x4d:
        raise ValueError("Old-format (pre-1996) SPC files are not supported.")
    if fversn == 0x4c:
        raise ValueError("Big-endian SPC files are not supported.")
    if fversn != 0x4b:
        raise ValueError(f"Not an SPC file (version byte {fversn:#x}).")
    fnpts, ffirst, flast, fnsub = struct.unpack_from("<iddi", data, 4)
    if ftflg & 0x04 or fnsub > 1:                       # TMULTI
        raise ValueError("Multifile SPC files are not supported.")
    if ftflg & 0x40:                                    # TXYXYS
        raise ValueError("XYXY-format SPC files are not supported.")
    if not 2 <= fnpts <= 10_000_000:
        raise ValueError(f"Implausible SPC point count ({fnpts}).")
    off = 512
    if ftflg & 0x80:                                    # TXVALS: explicit x
        x = np.frombuffer(data, "<f4", fnpts, off).astype(float)
        off += 4 * fnpts
    else:
        x = np.linspace(ffirst, flast, fnpts)
    subexp = struct.unpack_from("<b", data, off + 1)[0]
    exp = fexp if fexp != 0 else subexp
    off += 32
    if len(data) < off + 4 * fnpts:
        raise ValueError("Truncated SPC file.")
    if exp == -128:                                     # IEEE float32 y
        y = np.frombuffer(data, "<f4", fnpts, off).astype(float)
    else:
        y = np.frombuffer(data, "<i4", fnpts, off) * 2.0 ** (exp - 32)
    return make_spectrum(x, y, name=name,
                         x_label="Raman shift (cm-1)", source=source)


def _find_links(html: str, base_url: str, suffix: str) -> list:
    """All distinct absolute href targets in ``html`` ending in ``suffix``."""
    hrefs = re.findall(r"""(?:href|src)\s*=\s*["']?([^"'\s>]+)""",
                       html, re.IGNORECASE)
    out, seen = [], set()
    for h in hrefs:
        h = _html.unescape(h).strip()
        if not h.lower().endswith(suffix.lower()):
            continue
        absu = urllib.parse.urljoin(base_url, h)
        if absu not in seen:
            seen.add(absu)
            out.append(absu)
    return out


def _page_title(html: str) -> str:
    m = re.search(r"<title[^>]*>(.*?)</title>", html, re.IGNORECASE | re.DOTALL)
    if not m:
        return ""
    return re.sub(r"\s+", " ", _html.unescape(m.group(1))).strip()


def _ucl_entry_name(stem: str, title: str, n_on_page: int) -> str:
    """Canonical pigment name for one .spc file found on a pigment page."""
    if stem in UCL_STEM_NAMES:
        return UCL_STEM_NAMES[stem]
    if title.lower() in _UCL_NAMES_LOWER:
        name = _UCL_NAMES_LOWER[title.lower()]
        # several files on one page (e.g. verdigris variants): disambiguate
        return name if n_on_page == 1 else f"{name} ({stem})"
    if stem.lower() in _UCL_NAMES_LOWER:
        return _UCL_NAMES_LOWER[stem.lower()]
    return title or stem


def _wayback_raw(url: str) -> str:
    """For a web.archive.org URL, ask for the original bytes (``id_`` flag) -
    needed for binary .spc files, harmless to skip for HTML."""
    return re.sub(r"(://web\.archive\.org/web/\d{1,14})/", r"\1id_/", url,
                  count=1)


def ucl_cache_path() -> str:
    """Where the downloaded UCL library is cached (override the directory
    with the ORANGE_SPECTRA_CACHE environment variable)."""
    d = os.environ.get("ORANGE_SPECTRA_CACHE") or os.path.join(
        os.path.expanduser("~"), ".orange-spectra")
    return os.path.join(d, "ucl_pigments_raman.speclib")


def _fetch_ucl_from(base_url: str, fetch, progress=None) -> list:
    """Crawl one source (UCL itself or a Wayback snapshot) for the library."""
    index_html = fetch(base_url)[0].decode("latin-1", "replace")
    linked = _find_links(index_html, base_url + "/", ".html")
    pages = [p for p in linked if "/pigfiles/" in p.lower()]
    if not pages:
        # frameset / palette layout: pigment pages sit one level deeper
        for nav in linked[:40]:
            try:
                nav_html = fetch(nav)[0].decode("latin-1", "replace")
            except Exception:  # noqa: BLE001 - skip unreachable nav pages
                continue
            pages.extend(p for p in _find_links(nav_html, nav, ".html")
                         if "/pigfiles/" in p.lower())
        pages = list(dict.fromkeys(pages))          # dedupe, keep order
    if not pages:
        raise ValueError("no pigment pages found")

    entries, seen_spc = [], set()
    for i, page in enumerate(pages):
        if progress:
            progress(i, len(pages), page.rsplit("/", 1)[-1])
        try:
            page_html = fetch(page)[0].decode("latin-1", "replace")
        except Exception:  # noqa: BLE001 - skip unreachable pages
            continue
        title = _page_title(page_html)
        spc_urls = [u for u in _find_links(page_html, page, ".spc")
                    if u not in seen_spc]
        for url in spc_urls:
            seen_spc.add(url)
            stem = url.rsplit("/", 1)[-1]
            stem = stem[:-4].lower() if stem.lower().endswith(".spc") else stem
            try:
                blob = fetch(_wayback_raw(url))[0]
                s = read_spc(blob,
                             name=_ucl_entry_name(stem, title, len(spc_urls)),
                             source=url)
            except Exception:  # noqa: BLE001 - skip one bad file, keep going
                continue
            color = UCL_PIGMENT_COLORS.get(s["name"], "")
            if color:
                s["color"] = color
            entries.append(s)
    if progress:
        progress(len(pages), len(pages), "done")
    if not entries:
        raise ValueError("pigment pages found but no readable .spc spectra")
    return entries


def fetch_ucl_library(fetch=None, cache_path: str = None,
                      progress=None, force: bool = False) -> list:
    """Download the UCL pigment Raman library from UCL's website (first use
    only; afterwards the local cache is read).

    The pigment pages are discovered from the site's own HTML, so no file
    names are hard-coded beyond the index URL. ``fetch(url) -> (bytes, str)``
    is injectable for tests and mirrors; ``progress(done, total, label)`` is
    called as files arrive. Returns a list of spectrum dicts and writes a
    .speclib cache at ``cache_path`` (default: :func:`ucl_cache_path`).
    """
    cache_path = cache_path or ucl_cache_path()
    if not force and os.path.isfile(cache_path):
        return load_library(cache_path)
    fetch = fetch or default_fetch

    # UCL first (or a mirror via ORANGE_SPECTRA_UCL_BASE), then the Internet
    # Archive's copy - UCL's server 403s whole regions (observed from Taiwan).
    bases = [os.environ.get("ORANGE_SPECTRA_UCL_BASE") or UCL_BASE_URL]
    if UCL_WAYBACK_BASE not in bases:
        bases.append(UCL_WAYBACK_BASE)
    entries, errors = [], []
    for base in bases:
        try:
            entries = _fetch_ucl_from(base, fetch, progress)
        except Exception as exc:  # noqa: BLE001 - try the next source
            errors.append(f"{base}: {exc}")
            entries = []
        if entries:
            break
    if not entries:
        raise ValueError(
            "Could not download the UCL library from any source:\n  "
            + "\n  ".join(errors) +
            "\nIf your network is blocked by UCL, download the .spc files "
            "another way and run build_ucl_library_from_folder(folder).")

    entries.sort(key=lambda s: s["name"].lower())
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    save_library(entries, cache_path, name=UCL_LIBRARY_NAME)
    return entries


def build_ucl_library_from_folder(folder: str, cache_path: str = None) -> list:
    """Offline alternative to :func:`fetch_ucl_library`: build the library
    from a folder of .spc files downloaded manually from UCL's site. Writes
    the same cache, so the widget's built-in entry works afterwards."""
    entries = []
    for fn in sorted(os.listdir(folder)):
        if not fn.lower().endswith(".spc"):
            continue
        stem = os.path.splitext(fn)[0].lower()
        path = os.path.join(folder, fn)
        with open(path, "rb") as fh:
            s = read_spc(fh.read(), name=_ucl_entry_name(stem, "", 1),
                         source=path)
        color = UCL_PIGMENT_COLORS.get(s["name"], "")
        if color:
            s["color"] = color
        entries.append(s)
    if not entries:
        raise ValueError(f"No .spc files found in {folder!r}.")
    entries.sort(key=lambda s: s["name"].lower())
    if cache_path is None:
        cache_path = ucl_cache_path()
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    save_library(entries, cache_path, name=UCL_LIBRARY_NAME)
    return entries


def available_libraries() -> dict:
    """Name -> zero-argument loader for every library the widget can offer:
    .speclib files shipped in ``orangespectra/libraries/`` plus the
    download-on-first-use UCL pigment library."""
    libs = {name: (lambda p=path: load_library(p))
            for name, path in builtin_libraries().items()}
    libs.setdefault(UCL_LIBRARY_NAME, fetch_ucl_library)
    return libs

# ================================================================== mixture
def mixture_nnls(mixture: dict, references: list, fit_offset: bool = True) -> dict:
    """Solve mixture ≈ Σ cᵢ·refᵢ (+ offset) with cᵢ ≥ 0 (scipy NNLS)."""
    from scipy.optimize import nnls

    if not references:
        raise ValueError("Need at least one reference spectrum.")
    xm, ym = mixture["x"], mixture["y"]
    lo = max([xm[0]] + [r["x"][0] for r in references])
    hi = min([xm[-1]] + [r["x"][-1] for r in references])
    if hi <= lo:
        raise ValueError("Mixture and references do not share an x-range.")
    gx = xm[(xm >= lo) & (xm <= hi)]
    if gx.size < len(references) + 2:
        raise ValueError("Too little overlap to fit that many components.")
    b = np.interp(gx, xm, ym)
    cols = [np.interp(gx, r["x"], r["y"]) for r in references]
    if fit_offset:
        cols.append(np.ones_like(gx))
    A = np.column_stack(cols)
    coeffs, _ = nnls(A, b)
    fit = A @ coeffs
    offset = float(coeffs[-1]) if fit_offset else 0.0
    comp = coeffs[:-1] if fit_offset else coeffs
    total = comp.sum()
    fractions = comp / total if total > _EPS else comp
    ss_res = float(np.sum((b - fit) ** 2))
    ss_tot = float(np.sum((b - b.mean()) ** 2)) or 1.0
    return {"names": [r["name"] for r in references],
            "coeffs": comp, "fractions": fractions, "offset": offset,
            "x": gx, "fit": fit, "residual": b - fit,
            "r_squared": 1.0 - ss_res / ss_tot}


# ===================================================== common-grid resampling
def merge_spectra_union(spectra: list, max_points: int = 16384):
    """Put spectra onto one shared x-grid covering the UNION of their ranges.

    Unlike :func:`merge_spectra` (which needs a common overlap and crops to
    it), this never fails for disjoint ranges - y is NaN wherever a spectrum
    was not measured. Suited to reference libraries such as UCL's, whose 55
    pigments have no shared x-interval at all. The grid step is the finest
    native step among the inputs (bounded by ``max_points``).
    """
    spectra = [s for s in spectra if s["x"].size]
    if not spectra:
        raise ValueError("No spectra to merge.")
    x0 = spectra[0]["x"]
    if all(s["x"].size == x0.size and np.allclose(s["x"], x0) for s in spectra):
        return x0, [np.asarray(s["y"], float) for s in spectra]
    lo = min(float(s["x"][0]) for s in spectra)
    hi = max(float(s["x"][-1]) for s in spectra)
    if hi <= lo:
        raise ValueError("Spectra have no x-extent to merge.")
    steps = [float(np.median(np.diff(s["x"])))
             for s in spectra if s["x"].size > 1]
    step = max(min(steps) if steps else (hi - lo) / max_points,
               (hi - lo) / float(max_points))
    n = int(np.floor((hi - lo) / step + 0.5)) + 1
    gx = lo + (hi - lo) * np.arange(n) / max(n - 1, 1)
    ys = [np.interp(gx, s["x"], s["y"], left=np.nan, right=np.nan)
          for s in spectra]
    return gx, ys


def merge_spectra(spectra: list):
    """Put spectra onto one shared x-grid (overlap region, max point count)."""
    spectra = [s for s in spectra if s["x"].size]
    if not spectra:
        raise ValueError("No spectra to merge.")
    x0 = spectra[0]["x"]
    if all(s["x"].size == x0.size and np.allclose(s["x"], x0) for s in spectra):
        return x0, [s["y"] for s in spectra]
    lo = max(s["x"][0] for s in spectra)
    hi = min(s["x"][-1] for s in spectra)
    if hi <= lo:
        raise ValueError("Spectra do not share a common x-range.")
    n = max(s["x"].size for s in spectra)
    gx = np.linspace(lo, hi, n)
    return gx, [np.interp(gx, s["x"], s["y"]) for s in spectra]


# ================================================================ aquagram
# The 12 water matrix coordinates (WAMACs, nm) — Tsenkova aquaphotomics.
WAMACS = [1342.0, 1364.0, 1372.0, 1382.0, 1398.0, 1410.0,
          1438.0, 1444.0, 1464.0, 1474.0, 1492.0, 1516.0]

AQUAGRAM_NORMS = ["raw", "snv", "aquagram"]


def _snv(y: np.ndarray) -> np.ndarray:
    """Standard Normal Variate: centre and scale each spectrum by its own stats."""
    mu = y.mean()
    sd = y.std()
    return (y - mu) / sd if sd > _EPS else y - mu


def aquagram_coordinates(spectra: list, wamacs=None, normalization: str = "aquagram"):
    """Compute Aquagram coordinates (absorbance at the 12 WAMACs) for spectra.

    normalization:
      * "raw"      — absorbance sampled at each WAMAC, as-is.
      * "snv"      — SNV each spectrum first, then sample at the WAMACs.
      * "aquagram" — SNV, then standardise **across the sample set** at each WAMAC
                     (value − column mean) / column std. This is the classic
                     "normalized aquagram" (nirpyresearch / Tsenkova): 0 = the
                     group average, ± = above/below average water absorbance.

    Returns {wamacs, names, values (n×12), normalization}.
    """
    if normalization not in AQUAGRAM_NORMS:
        raise ValueError(f"normalization must be one of {AQUAGRAM_NORMS}")
    bands = np.asarray(wamacs if wamacs is not None else WAMACS, float)
    if bands.size < 3:
        raise ValueError("Need at least 3 WAMACs bands.")
    specs = [s for s in spectra if s["x"].size >= 2]
    if not specs:
        raise ValueError("No usable spectra.")

    # Warn (softly) if the WAMACs fall outside a spectrum's measured range: np.interp
    # clamps at the ends, so out-of-range bands would all read the edge value.
    lo = max(s["x"][0] for s in specs)
    hi = min(s["x"][-1] for s in specs)
    covered = bool(bands.min() >= lo - 1e-9 and bands.max() <= hi + 1e-9)

    rows = []
    for s in specs:
        y = _snv(s["y"]) if normalization in ("snv", "aquagram") else s["y"]
        rows.append(np.interp(bands, s["x"], y))
    values = np.vstack(rows)

    if normalization == "aquagram":
        mu = values.mean(axis=0)
        sd = values.std(axis=0)
        sd = np.where(sd > _EPS, sd, 1.0)
        values = (values - mu) / sd

    return {"wamacs": bands, "names": [s["name"] for s in specs],
            "values": values, "normalization": normalization, "covered": covered}


# ============================================================ peak finding
def find_spectrum_peaks(x, y, min_height_frac: float = 0.05,
                        min_prominence_frac: float = 0.03,
                        min_distance: float = 0.0,
                        smooth_window: int = 0) -> list:
    """Detect peaks and measure FWHM (same method as SpectraView's Analyze menu).

    Fractions are relative to the signal's full range; ``min_distance`` is in
    x-units; ``smooth_window`` (odd, >=3) applies a Savitzky-Golay pre-smooth
    for detection only (heights are read from the raw signal). Returns dicts
    (center, height, fwhm, prominence, area, index) sorted by center.
    """
    from scipy.signal import find_peaks as _sp_find_peaks
    from scipy.signal import peak_widths, savgol_filter

    x = np.asarray(x, float)
    y0 = np.asarray(y, float)
    order = np.argsort(x)
    x, y0 = x[order], y0[order]
    y = y0
    if smooth_window and 3 <= smooth_window < y0.size:
        w = smooth_window + (smooth_window + 1) % 2      # force odd
        y = savgol_filter(y0, w, min(3, w - 1), mode="interp")
    rng = float(y.max() - y.min()) or 1.0
    height = y.min() + min_height_frac * rng
    prominence = max(min_prominence_frac * rng, _EPS)
    dx = float(np.median(np.abs(np.diff(x)))) or 1.0
    distance = max(1, int(round(min_distance / dx))) if min_distance > 0 else 1

    idx, props = _sp_find_peaks(y, height=height, prominence=prominence,
                                distance=distance)
    if idx.size == 0:
        return []
    widths, _, lips, rips = peak_widths(y, idx, rel_height=0.5)
    axis = np.arange(x.size)
    fwhm = np.abs(np.interp(rips, axis, x) - np.interp(lips, axis, x))

    gauss_area = np.sqrt(np.pi / (4.0 * np.log(2.0)))    # ~1.0645 * h * FWHM
    peaks = []
    for k, i in enumerate(idx):
        h, f = float(y0[i]), float(fwhm[k])
        peaks.append({"center": float(x[i]), "height": h, "fwhm": f,
                      "prominence": float(props["prominences"][k]),
                      "area": abs(h) * f * gauss_area, "index": int(i)})
    peaks.sort(key=lambda p: p["center"])
    return peaks


# ============================================================ PLS-DA
def plsda_fit(X, labels, n_components: int = 2) -> dict:
    """PLS-DA: PLS2 (NIPALS) regression of one-hot classes on spectra.

    Pure-numpy NIPALS keeps this dependency-free and deterministic. Returns
    scores T (n x A), x-loadings P (p x A), weights W, y-loadings Q (c x A),
    VIP scores (p), per-sample predicted class + soft y-hat, training accuracy
    and confusion matrix (rows = true class).
    """
    X = np.asarray(X, float)
    labels = [str(v) for v in labels]
    if X.ndim != 2 or X.shape[0] != len(labels):
        raise ValueError("X must be (n_samples, n_features) matching labels.")
    classes = sorted(set(labels))
    if len(classes) < 2:
        raise ValueError("PLS-DA needs at least 2 classes.")
    n, p = X.shape
    A = int(max(1, min(n_components, n - 1, p)))
    Y = np.zeros((n, len(classes)))
    for i, lab in enumerate(labels):
        Y[i, classes.index(lab)] = 1.0

    x_mean, y_mean = X.mean(axis=0), Y.mean(axis=0)
    Xc, Yc = X - x_mean, Y - y_mean
    T = np.zeros((n, A))
    W = np.zeros((p, A))
    P = np.zeros((p, A))
    Q = np.zeros((len(classes), A))
    Xa, Ya = Xc.copy(), Yc.copy()
    for a in range(A):
        u = Ya[:, int(np.argmax(Ya.var(axis=0)))]
        for _ in range(500):
            w = Xa.T @ u
            w /= (np.linalg.norm(w) or 1.0)
            t = Xa @ w
            q = Ya.T @ t / max(t @ t, _EPS)
            u_new = Ya @ q / max(q @ q, _EPS)
            if np.linalg.norm(u_new - u) <= 1e-10 * max(np.linalg.norm(u), 1.0):
                u = u_new
                break
            u = u_new
        t = Xa @ w
        pa = Xa.T @ t / max(t @ t, _EPS)
        q = Ya.T @ t / max(t @ t, _EPS)
        Xa = Xa - np.outer(t, pa)
        Ya = Ya - np.outer(t, q)
        T[:, a], W[:, a], P[:, a], Q[:, a] = t, w, pa, q

    # Regression coefficients B = W (P'W)^-1 Q'
    B = W @ np.linalg.solve(P.T @ W, Q.T)
    y_hat = Xc @ B + y_mean
    pred_idx = np.argmax(y_hat, axis=1)
    predicted = [classes[i] for i in pred_idx]
    truth_idx = np.array([classes.index(lab) for lab in labels])
    accuracy = float(np.mean(pred_idx == truth_idx))
    confusion = np.zeros((len(classes), len(classes)), int)
    for ti, pi in zip(truth_idx, pred_idx):
        confusion[ti, pi] += 1

    # VIP_j = sqrt( p * sum_a(ssy_a * (w_ja/||w_a||)^2) / sum_a ssy_a )
    ssy = np.array([(T[:, a] @ T[:, a]) * (Q[:, a] @ Q[:, a]) for a in range(A)])
    wnorm2 = np.maximum((W ** 2).sum(axis=0), _EPS)
    vip = np.sqrt(p * ((W ** 2) / wnorm2 @ ssy) / max(ssy.sum(), _EPS))

    # Explained X variance per component
    total = float((Xc ** 2).sum()) or 1.0
    xvar = np.array([(np.outer(T[:, a], P[:, a]) ** 2).sum() / total
                     for a in range(A)])
    return {"classes": classes, "n_components": A, "scores": T,
            "loadings": P, "weights": W, "y_loadings": Q, "vip": vip,
            "coef": B, "x_mean": x_mean, "y_mean": y_mean,
            "y_hat": y_hat, "predicted": predicted, "accuracy": accuracy,
            "confusion": confusion, "explained_x_variance": xvar}
