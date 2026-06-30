"""Reader for SOPRANO spectral-library pages.

SOPRANO embeds its plotted spectrum directly in the page as a Dygraph data
array.  We parse that array instead of automating the download modal.
"""
from __future__ import annotations

import ast
import html
import re
import urllib.parse
import urllib.request

import numpy as np

from ..spectrum import Spectrum


_DYGRAPH_RE = re.compile(r"new\s+Dygraph\s*\(", re.I)


def _strip_tags(text: str) -> str:
    text = re.sub(r"<[^>]+>", "", text)
    return html.unescape(text).strip()


def _find_balanced(text: str, start: int, opener: str = "[", closer: str = "]") -> str:
    depth = 0
    in_str: str | None = None
    escape = False
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
    raise ValueError("Could not find complete SOPRANO Dygraph data array.")


def _soprano_units(text: str) -> tuple[str, str]:
    lower = text.lower()
    if "raman" in lower:
        x_unit = "raman_cm-1"
    elif "cm" in lower and "-1" in lower:
        x_unit = "cm-1"
    elif "nm" in lower:
        x_unit = "nm"
    else:
        x_unit = "pixel"
    y_unit = "intensity" if "intensity" in lower else "a.u."
    return x_unit, y_unit


def parse_soprano_html(text: str, source: str = "") -> Spectrum:
    """Parse one SOPRANO page into a :class:`Spectrum`."""
    match = _DYGRAPH_RE.search(text)
    if not match:
        raise ValueError("No SOPRANO Dygraph spectrum found.")
    data_start = text.find("[[", match.end())
    if data_start < 0:
        raise ValueError("No SOPRANO spectrum data array found.")
    data_literal = _find_balanced(text, data_start)
    rows = ast.literal_eval(data_literal)
    arr = np.asarray(rows, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 2 or arr.shape[0] < 2:
        raise ValueError("SOPRANO spectrum data has an unexpected shape.")

    labels = []
    labels_match = re.search(r"labels\s*:\s*(\[[^\]]+\])", text, re.S)
    if labels_match:
        try:
            labels = [_strip_tags(str(v)) for v in ast.literal_eval(labels_match.group(1))]
        except Exception:
            labels = []
    axis_options = [_strip_tags(v) for _, v in re.findall(
        r"\b(?:x|y)label\s*:\s*([\"'])(.*?)\1", text, re.S | re.I)]

    title_bits = [_strip_tags(v) for v in re.findall(
        r"<h[12][^>]*>(.*?)</h[12]>", text, re.S | re.I)]
    name = labels[1] if len(labels) > 1 and labels[1] else ""
    subtitle = next((v for v in title_bits if "raman" in v.lower()), "")
    if not name and title_bits:
        name = title_bits[0]
    if subtitle and subtitle not in name:
        name = f"{name} {subtitle}".strip()
    if not name:
        parsed = urllib.parse.urlparse(source)
        qs = urllib.parse.parse_qs(parsed.query)
        name = qs.get("id", ["SOPRANO spectrum"])[0]

    axis_text = " ".join(labels + axis_options + title_bits)
    x_unit, y_unit = _soprano_units(axis_text)
    parsed_url = urllib.parse.urlparse(source)
    qs = urllib.parse.parse_qs(parsed_url.query)
    meta = {
        "source": source,
        "format": "SOPRANO page",
        "soprano_id": qs.get("id", [""])[0],
        "dataset": qs.get("ds", ["baseline corrected"])[0],
        "library": qs.get("lib", [""])[0],
    }
    return Spectrum(x=arr[:, 0], y=arr[:, 1], name=name,
                    x_unit=x_unit, y_unit=y_unit, meta=meta)


def load_soprano_url(url: str, timeout: float = 20.0) -> list[Spectrum]:
    """Fetch and read a SOPRANO spectral-library URL."""
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "SpectraView/1.0 (+https://github.com/Tai-ShengYeh/spectraview)"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        charset = response.headers.get_content_charset() or "utf-8"
        text = response.read().decode(charset, errors="replace")
    return [parse_soprano_html(text, source=url)]
