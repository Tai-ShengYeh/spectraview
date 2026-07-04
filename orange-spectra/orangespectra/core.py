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
    lm = re.search(r"labels\s*:\s*(\[[^\]]+\])", text, re.S)
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
def parse_csv(text: str, source: str = ""):
    rows = []
    header = None
    for line in text.splitlines():
        parts = re.split(r"[,\t;]+", line.strip())
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


def load_spectrum_url(id_or_url, fetch=default_fetch):
    """Fetch + auto-detect: IRUG page, SOPRANO page, JCAMP-DX, or CSV."""
    url = resolve_source(id_or_url)
    raw, _ctype = fetch(url)
    text = raw.decode("utf-8", errors="replace")
    for parser in (parse_irug_jqplot, parse_soprano, parse_jcamp, parse_csv):
        spec = parser(text, source=url)
        if spec is not None and spec["x"].size >= 2:
            return spec
    raise ValueError(
        f"Could not extract a spectrum from {url}\n"
        "Supported: IRUG detail pages, SOPRANO pages, JCAMP-DX (AFFN), CSV/TSV."
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
    obj = {"name": name, "library": [
        {"name": s["name"], "x_unit": _unit_from_label(s.get("x_label", "")),
         "y_unit": "intensity", "x": np.asarray(s["x"]).tolist(),
         "y": np.asarray(s["y"]).tolist()} for s in entries]}
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=1)


def load_library(path: str) -> list:
    """Read a SpectraView .speclib (JSON) file into spectrum dicts."""
    with open(path, "r", encoding="utf-8") as fh:
        obj = json.load(fh)
    out = []
    for e in obj.get("library", []):
        out.append(make_spectrum(
            e["x"], e["y"], name=e.get("name", "ref"),
            x_label=_label_from_unit(e.get("x_unit", "")), source=path))
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
    T = np.zeros((n, A)); W = np.zeros((p, A))
    P = np.zeros((p, A)); Q = np.zeros((len(classes), A))
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
