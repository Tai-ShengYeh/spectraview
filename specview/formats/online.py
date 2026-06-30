"""Import spectra straight from the web (IRUG database or any direct URL).

The single entry point used by the UI is :func:`load_online`. Give it:

  * an IRUG spectrum id            -> ``"3537"`` or ``3537``
  * an IRUG detail-page URL        -> ``"http://www.irug.org/jcamp-details?id=3537"``
  * a direct JCAMP-DX / CSV URL    -> ``"https://example.org/spectrum.jdx"``

It downloads the data, writes it to a temporary file with a sensible
extension, and hands it to :func:`specview.formats.load_any` so the existing
JCAMP-DX / ASCII readers do the parsing. Network access uses only the Python
standard library (``urllib``) — no extra dependency.

The actual HTTP fetch is injectable (the ``fetch`` argument) so the importer
can be unit-tested offline.
"""
from __future__ import annotations

import os
import re
import tempfile
import urllib.request

from ..spectrum import Spectrum

USER_AGENT = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/124.0 Safari/537.36"
)

#: IRUG detail-page URL for a numeric spectrum id.
IRUG_DETAIL_URL = "http://www.irug.org/jcamp-details?id={id}"


# --------------------------------------------------------------------------- #
# Source resolution
# --------------------------------------------------------------------------- #
def resolve_source(id_or_url) -> str:
    """Turn a user-supplied id / URL into a fetchable URL.

    A bare number (or ``irug:123``) becomes the IRUG detail-page URL; anything
    that already looks like a URL is returned unchanged.
    """
    s = str(id_or_url).strip()
    if not s:
        raise ValueError("Empty source. Enter an IRUG id or a URL.")
    m = re.fullmatch(r"(?:irug:)?(\d+)", s, re.IGNORECASE)
    if m:
        return IRUG_DETAIL_URL.format(id=m.group(1))
    if not re.match(r"https?://", s, re.IGNORECASE):
        raise ValueError(
            f"Not an IRUG id or URL: {s!r}. "
            "Enter a number (IRUG id) or a full http(s) URL."
        )
    return s


# --------------------------------------------------------------------------- #
# HTTP
# --------------------------------------------------------------------------- #
def _default_fetch(url: str) -> tuple[bytes, str]:
    """GET ``url`` with a browser-like User-Agent.

    Returns ``(body_bytes, content_type)``. Uses ``requests`` if it happens to
    be installed, otherwise falls back to :mod:`urllib`.
    """
    try:
        import requests  # type: ignore

        r = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=60)
        r.raise_for_status()
        return r.content, r.headers.get("Content-Type", "")
    except ImportError:
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=60) as resp:  # noqa: S310
            return resp.read(), resp.headers.get("Content-Type", "") or ""


# --------------------------------------------------------------------------- #
# Content sniffing
# --------------------------------------------------------------------------- #
def _looks_like_jcamp(text: str) -> bool:
    head = text.lstrip()[:4000].upper()
    return "##TITLE=" in head or "##XYDATA=" in head or "##XYPOINTS=" in head


def _looks_like_html(text: str, content_type: str) -> bool:
    if "html" in content_type.lower():
        return True
    head = text.lstrip()[:200].lower()
    return head.startswith("<!doctype") or head.startswith("<html") or "<head" in head


def _extract_inline_jcamp(html: str) -> str | None:
    """Return a JCAMP-DX block embedded directly in a page, if present."""
    m = re.search(r"(##TITLE\s*=.*?##END\s*=\s*)", html, re.DOTALL | re.IGNORECASE)
    return m.group(1) if m else None


def _find_data_link(html: str, base_url: str) -> str | None:
    """Find a link to a downloadable spectrum file in page HTML."""
    candidates = re.findall(
        r"""(?:href|src|data-[\w-]+)\s*=\s*['"]([^'"]+?\."""
        r"""(?:jdx|dx|jcamp|jcm|csv|txt|spc))(?:\?[^'"]*)?['"]""",
        html,
        re.IGNORECASE,
    )
    if not candidates:
        return None
    return urllib.request.urljoin(base_url, candidates[0])


def _ext_for(text: str, url: str) -> str:
    """Pick a temp-file extension so ``load_any`` dispatches to the right reader."""
    if _looks_like_jcamp(text):
        return ".jdx"
    path = urllib.request.urlparse(url).path.lower()
    for ext in (".jdx", ".dx", ".jcamp", ".jcm", ".csv", ".txt", ".spc", ".asc"):
        if path.endswith(ext):
            return ext
    return ".csv"  # default: let the ASCII reader try


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #
def load_online(id_or_url, fetch=_default_fetch) -> list[Spectrum]:
    """Download a spectrum from IRUG (or any URL) and parse it into Spectra.

    Parameters
    ----------
    id_or_url:
        IRUG spectrum id, IRUG detail-page URL, or a direct file URL.
    fetch:
        Callable ``url -> (bytes, content_type)``. Injectable for testing.

    Returns the parsed list of :class:`Spectrum`, with each spectrum's
    ``meta['source']`` set to the URL it came from.
    """
    # Local import avoids a circular import (formats/__init__ imports nothing
    # from here, but load_any lives there and pulls in every reader).
    from . import load_any

    url = resolve_source(id_or_url)
    raw, ctype = fetch(url)
    text = raw.decode("utf-8", errors="replace")
    source_url = url

    # If we landed on an HTML detail page, dig out the spectrum data.
    if _looks_like_html(text, ctype) and not _looks_like_jcamp(text):
        inline = _extract_inline_jcamp(text)
        if inline:
            text, raw = inline, inline.encode("utf-8")
        else:
            link = _find_data_link(text, url)
            if not link:
                raise ValueError(
                    "Could not find a downloadable spectrum on the page:\n"
                    f"  {url}\n"
                    "Open it in a browser, use its Download button to save the "
                    "JCAMP-DX/CSV file, then load that file via File ▸ Open."
                )
            raw, ctype = fetch(link)
            text = raw.decode("utf-8", errors="replace")
            source_url = link

    ext = _ext_for(text, source_url)
    tmp_dir = tempfile.mkdtemp(prefix="specview_online_")
    base = re.sub(r"[^\w.-]+", "_", os.path.basename(
        urllib.request.urlparse(source_url).path) or "download") or "download"
    if not base.lower().endswith(ext):
        base += ext
    tmp_path = os.path.join(tmp_dir, base)
    with open(tmp_path, "wb") as fh:
        fh.write(raw)

    spectra = load_any(tmp_path)
    for s in spectra:
        s.meta = dict(s.meta or {})
        s.meta["source"] = source_url
    if not spectra:
        raise ValueError(f"No spectral data parsed from {source_url}")
    return spectra


__all__ = ["load_online", "resolve_source", "IRUG_DETAIL_URL"]
