"""Headless tests for the download-on-first-use UCL pigment library
(no Orange / Qt, no network - the UCL site is simulated with a fake fetch).

Run:  python tests/test_builtin_library.py
"""
import os
import struct
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from orangespectra import core  # noqa: E402

PASS = FAIL = 0


def check(name, cond):
    global PASS, FAIL
    PASS, FAIL = PASS + cond, FAIL + (not cond)
    print(f"  {'ok  ' if cond else 'FAIL'} {name}")


# ------------------------------------------------ tiny SPC writer for tests
def make_spc(x, y, float_y=True, descending=False, exp=13):
    """Minimal new-format (0x4b) single-sub SPC with an even x grid."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if descending:
        x, y = x[::-1], y[::-1]
    head = bytearray(512)
    struct.pack_into("<BBBb", head, 0, 0, 0x4b, 0,
                     -128 if float_y else exp)
    struct.pack_into("<iddi", head, 4, x.size, x[0], x[-1], 1)
    sub = bytes(32)
    if float_y:
        body = np.asarray(y, "<f4").tobytes()
    else:
        body = np.asarray(np.round(y / 2.0 ** (exp - 32)), "<i4").tobytes()
    return bytes(head) + sub + body


print("== read_spc round-trip ==")
gx = np.linspace(100, 1200, 551)
gy = 40 + 900 * np.exp(-0.5 * ((gx - 1008) / 4.0) ** 2)
s = core.read_spc(make_spc(gx, gy), name="t")
check("float32 y round-trips", np.allclose(s["y"], gy, rtol=1e-6)
      and np.allclose(s["x"], gx))
s2 = core.read_spc(make_spc(gx, gy, float_y=False))
check("int+exponent y round-trips", np.allclose(s2["y"], gy, atol=1e-4))
s3 = core.read_spc(make_spc(gx, gy, descending=True))
check("descending x comes back ascending",
      np.all(np.diff(s3["x"]) > 0) and np.allclose(s3["y"], gy, rtol=1e-6))
try:
    core.read_spc(b"\x00\x4d" + bytes(600))
    check("old-format SPC rejected", False)
except ValueError as exc:
    check("old-format SPC rejected", "Old-format" in str(exc))

print("== pigment name/colour tables ==")
check("55 canonical pigments", len(core.UCL_PIGMENT_COLORS) == 55)
check("7 colour groups",
      set(core.UCL_PIGMENT_COLORS.values())
      == {"black", "blue", "green", "orange", "red", "white", "yellow"})
check("stem map targets are canonical",
      all(v in core.UCL_PIGMENT_COLORS for v in core.UCL_STEM_NAMES.values()))

# ------------------------------------------------------- fake UCL web site
gy2 = 30 + 500 * np.exp(-0.5 * ((gx - 1086) / 5.0) ** 2)
SITE = {
    core.UCL_BASE_URL: (b"""<html><body>
        <a href="pigfiles/gypsum.html">gypsum</a>
        <a href="pigfiles/verdigri.html">verdigris</a>
        <a href="refs.html">references</a></body></html>""", "text/html"),
    core.UCL_BASE_URL + "/pigfiles/gypsum.html": (
        b"<html><head><title>Gypsum</title></head>"
        b"<body><a href='../spectra/gypsum.spc'>spectrum</a></body></html>",
        "text/html"),
    core.UCL_BASE_URL + "/pigfiles/verdigri.html": (
        b"<html><head><title>Verdigris</title></head><body>"
        b"<a href='../spectra/verdiraw.spc'>raw</a></body></html>",
        "text/html"),
    core.UCL_BASE_URL + "/spectra/gypsum.spc": (make_spc(gx, gy), "app"),
    core.UCL_BASE_URL + "/spectra/verdiraw.spc": (make_spc(gx, gy2), "app"),
}
CALLS = []


def fake_fetch(url):
    CALLS.append(url)
    if url not in SITE:
        raise OSError(f"404 {url}")
    return SITE[url]


print("== fetch_ucl_library against the fake site ==")
cache = os.path.join(tempfile.mkdtemp(), "ucl.speclib")
seen_progress = []
entries = core.fetch_ucl_library(fetch=fake_fetch, cache_path=cache,
                                 progress=lambda d, t, s: seen_progress.append(d))
check("two spectra found", len(entries) == 2)
check("names resolved canonically (title + stem map)",
      [e["name"] for e in entries] == ["Gypsum", "Verdigris (raw)"])
check("colours attached", [e.get("color") for e in entries]
      == ["white", "green"])
check("source is the UCL url",
      entries[0]["source"].endswith("/spectra/gypsum.spc"))
check("data survived download", np.allclose(entries[0]["y"], gy, rtol=1e-6))
check("progress was reported", len(seen_progress) >= 2)
check("cache file written", os.path.isfile(cache))

print("== cache short-circuits the network ==")
def no_net(url):
    raise AssertionError("network touched despite cache")
cached = core.fetch_ucl_library(fetch=no_net, cache_path=cache)
check("loaded from cache, no fetch", len(cached) == 2)
check("colour survives the .speclib cache",
      [e.get("color") for e in cached] == ["white", "green"])
check("y identical through cache",
      np.allclose(cached[0]["y"], entries[0]["y"], atol=1e-6))

print("== offline folder fallback ==")
folder = tempfile.mkdtemp()
open(os.path.join(folder, "Gypsum.spc"), "wb").write(make_spc(gx, gy))
open(os.path.join(folder, "Leadwhit.spc"), "wb").write(make_spc(gx, gy2))
cache2 = os.path.join(tempfile.mkdtemp(), "ucl2.speclib")
offline = core.build_ucl_library_from_folder(folder, cache_path=cache2)
check("folder build finds both", [e["name"] for e in offline]
      == ["Gypsum", "Lead White"])
check("folder build writes cache", os.path.isfile(cache2))

print("== frameset / palette indirection ==")
SITE2 = {
    core.UCL_BASE_URL: (
        b"<html><frameset><frame src='palette.html'></frameset></html>",
        "text/html"),
    core.UCL_BASE_URL + "/palette.html": (
        b"<a href='pigfiles/gypsum.html'>gypsum</a>", "text/html"),
    core.UCL_BASE_URL + "/pigfiles/gypsum.html": (
        SITE[core.UCL_BASE_URL + "/pigfiles/gypsum.html"][0], "text/html"),
    core.UCL_BASE_URL + "/spectra/gypsum.spc": (make_spc(gx, gy), "app"),
}


def fake_fetch2(url):
    if url not in SITE2:
        raise OSError(f"404 {url}")
    return SITE2[url]


cache3 = os.path.join(tempfile.mkdtemp(), "ucl3.speclib")
deep = core.fetch_ucl_library(fetch=fake_fetch2, cache_path=cache3)
check("pigment pages found one level below a frameset",
      [e["name"] for e in deep] == ["Gypsum"])

print("== union merge: disjoint ranges must not crash (UCL crash fix) ==")
a = core.make_spectrum(np.linspace(100, 400, 301), np.ones(301), name="low")
b = core.make_spectrum(np.linspace(1100, 1900, 801), 2 * np.ones(801),
                       name="high")
try:
    core.merge_spectra([a, b])
    check("overlap merge rejects disjoint ranges (unchanged)", False)
except ValueError:
    check("overlap merge rejects disjoint ranges (unchanged)", True)
ugx, uys = core.merge_spectra_union([a, b])
check("union grid spans both ranges", ugx[0] <= 100 and ugx[-1] >= 1900)
ua, ub = uys
in_a = (ugx >= 100) & (ugx <= 400)
check("values kept inside measured range",
      np.allclose(ua[in_a], 1.0) and np.isnan(ua[~in_a]).all())
check("second spectrum likewise",
      np.allclose(ub[(ugx >= 1100) & (ugx <= 1900)], 2.0))
same_gx, same_ys = core.merge_spectra_union([a, a])
check("identical grids pass through exactly",
      same_gx.size == 301 and np.allclose(same_gx, a["x"]))

print("== wayback fallback when UCL blocks the region ==")
WB = core.UCL_WAYBACK_BASE
TS = "https://web.archive.org/web/20230601120000"
ORIG = "https://www.chem.ucl.ac.uk/resources/raman"
SITE3 = {
    WB: (("<a href='/web/20230601120000/" + ORIG +
          "/pigfiles/gypsum.html'>gypsum</a>").encode(), "text/html"),
    f"{TS}/{ORIG}/pigfiles/gypsum.html": (
        (b"<html><head><title>Gypsum</title></head><body>"
         b"<a href='/web/20230601120000/" + ORIG.encode() +
         b"/spectra/gypsum.spc'>spc</a></body></html>"), "text/html"),
    # binary is only served with the id_ (raw bytes) flag
    f"{TS}id_/{ORIG}/spectra/gypsum.spc": (make_spc(gx, gy), "app"),
}


def fake_fetch3(url):
    if url.startswith(ORIG):
        raise OSError("403 Forbidden (nginx)")       # UCL blocks this region
    if url not in SITE3:
        raise OSError(f"404 {url}")
    return SITE3[url]


cache4 = os.path.join(tempfile.mkdtemp(), "ucl4.speclib")
wb = core.fetch_ucl_library(fetch=fake_fetch3, cache_path=cache4)
check("library built from the Wayback snapshot",
      [e["name"] for e in wb] == ["Gypsum"])
check("binary fetched with id_ raw-bytes flag",
      wb[0]["source"].endswith("gypsum.spc"))
check("data intact via wayback", np.allclose(wb[0]["y"], gy, rtol=1e-6))

print("== both sources down -> helpful error ==")
def all_fail(url):
    raise OSError("403")
try:
    core.fetch_ucl_library(fetch=all_fail,
                           cache_path=os.path.join(tempfile.mkdtemp(), "x.speclib"))
    check("raises with guidance", False)
except ValueError as exc:
    check("raises with guidance",
          "build_ucl_library_from_folder" in str(exc))

print("== available_libraries ==")
libs = core.available_libraries()
check("UCL library is offered", core.UCL_LIBRARY_NAME in libs)
check("loaders are callable", all(callable(v) for v in libs.values()))

print("== search sanity on the fetched entries ==")
hits = core.search_library(entries[0], entries, rank_by="correlation")
check("self-match first", hits[0]["name"] == "Gypsum"
      and hits[0]["scores"]["correlation"] > 0.999999)

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)
