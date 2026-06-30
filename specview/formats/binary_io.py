"""Readers for binary instrument formats (optional dependencies).

These degrade gracefully: if the needed package is not installed, a clear
ImportError with a ``pip install`` hint is raised instead of crashing.
"""
from __future__ import annotations

import os
import struct

import numpy as np

from ..spectrum import Spectrum


class MissingDependency(ImportError):
    """Raised when an optional reader's backend package is not installed."""


def _is_shimadzu_spc(head: bytes) -> bool:
    """Shimadzu UVProbe ``.SPC`` (UV-Vis) starts with 0x00 0x16; GRAMS/Galactic
    ``.spc`` instead carries a version byte (0x4B/0x4C/0x4D) at offset 1."""
    return len(head) >= 2 and head[0] == 0x00 and head[1] == 0x16


def load_shimadzu_spc(path: str) -> list[Spectrum]:
    """Shimadzu UV-Vis ``.SPC`` (e.g. UV-1900). Pure-Python, no dependency.

    Little-endian layout: a 120-byte header (firstX float32 @10, lastX float32
    @14) then the float32 values to end-of-file. The point count is taken from
    the file size — there is no explicit NPOINTS field. The .SPC carries no
    unit tag, so X is read as nm and Y as absorbance (or %T if values exceed
    ~10, which absorbance never does but transmittance-% routinely does).
    """
    with open(path, "rb") as fh:
        raw = fh.read()
    if len(raw) < 124 or not _is_shimadzu_spc(raw[:2]):
        raise ValueError(f"{os.path.basename(path)} is not a Shimadzu .SPC file.")
    first_x = float(struct.unpack_from("<f", raw, 10)[0])
    last_x = float(struct.unpack_from("<f", raw, 14)[0])
    body = len(raw) - 120
    if body <= 0 or body % 4:
        raise ValueError(
            f"Unexpected Shimadzu .SPC body length in {os.path.basename(path)}.")
    n = body // 4
    y = np.frombuffer(raw, "<f4", n, 120).astype(float)
    x = np.linspace(first_x, last_x, n)
    y_unit = "%T" if float(np.nanmax(y)) > 10.0 else "absorbance"
    base = os.path.splitext(os.path.basename(path))[0]
    return [Spectrum(x=x, y=y, name=base, x_unit="nm", y_unit=y_unit,
                     meta={"source": path, "format": "Shimadzu .SPC"})]


def load_spc(path: str) -> list[Spectrum]:
    """GRAMS/Galactic .spc (``pip install spc-spectra``) or Shimadzu UV-Vis .SPC.

    The ``.spc`` extension is shared by two unrelated formats; Shimadzu UVProbe
    files are detected by signature and parsed natively (no dependency needed).
    """
    with open(path, "rb") as fh:
        if _is_shimadzu_spc(fh.read(8)):
            return load_shimadzu_spc(path)
    try:
        import spc_spectra as spc  # type: ignore
    except ImportError:
        try:
            import spc  # type: ignore  # older package name
        except ImportError as exc:
            raise MissingDependency(
                "Reading .spc files needs the 'spc-spectra' package.\n"
                "Install it with:  pip install spc-spectra"
            ) from exc

    f = spc.File(path)
    base = os.path.splitext(os.path.basename(path))[0]
    # Map SPC unit codes to our axis units (best-effort).
    xlabel = (getattr(f, "xlabel", "") or "").lower()
    ylabel = (getattr(f, "ylabel", "") or "").lower()
    x_unit = ("cm-1" if "cm-1" in xlabel or "wavenumber" in xlabel else
              "nm" if "nm" in xlabel or "nanomet" in xlabel else
              "raman_cm-1" if "raman" in xlabel else "pixel")
    y_unit = ("absorbance" if "absorb" in ylabel else
              "transmittance" if "trans" in ylabel else
              "counts" if "count" in ylabel else "intensity")

    spectra = []
    for i, sub in enumerate(getattr(f, "sub", [])):
        x = np.asarray(getattr(f, "x", getattr(sub, "x", [])), dtype=float)
        y = np.asarray(sub.y, dtype=float)
        if x.size != y.size:
            x = np.arange(y.size, dtype=float)
        name = base if len(f.sub) == 1 else f"{base} [{i + 1}]"
        spectra.append(Spectrum(x=x, y=y, name=name, x_unit=x_unit, y_unit=y_unit,
                                meta={"source": path}))
    if not spectra:
        raise ValueError(f"No subfiles found in {os.path.basename(path)}.")
    return spectra


_OPUS_HEADER_LEN = 504
_OPUS_FIRST_CURSOR = 24
_OPUS_META_BLOCK_SIZE = 12
_OPUS_PARAM_TYPES = {0: "int", 1: "float", 2: "str", 3: "str", 4: "str"}


def _opus_entries(raw: bytes) -> list[dict]:
    """Read the OPUS block directory from the fixed 504-byte header."""
    if len(raw) < _OPUS_HEADER_LEN or raw[:4] != b"\x0a\x0a\xfe\xfe":
        raise ValueError("Not a Bruker OPUS file.")
    entries = []
    for cursor in range(_OPUS_FIRST_CURSOR, _OPUS_HEADER_LEN, _OPUS_META_BLOCK_SIZE):
        data_type = raw[cursor]
        channel_type = raw[cursor + 1]
        text_type = raw[cursor + 2]
        chunk_size = struct.unpack_from("<I", raw, cursor + 4)[0]
        offset = struct.unpack_from("<I", raw, cursor + 8)[0]
        if offset <= 0:
            break
        stop = offset + 4 * chunk_size
        if stop > len(raw):
            raise ValueError("Invalid OPUS block directory.")
        entries.append({
            "data_type": data_type, "channel_type": channel_type,
            "text_type": text_type, "chunk_size": chunk_size,
            "offset": offset, "stop": stop,
        })
    return entries


def _opus_parse_param(raw: bytes, entry: dict) -> dict:
    """Parse an OPUS parameter block into a small dict."""
    chunk = raw[entry["offset"]:entry["stop"]]
    params = {}
    cursor = 0
    while cursor + 8 <= len(chunk):
        name = chunk[cursor:cursor + 3].decode("latin1", errors="replace")
        if name == "END":
            break
        type_index = struct.unpack_from("<H", chunk, cursor + 4)[0]
        param_type = _OPUS_PARAM_TYPES.get(type_index)
        param_size = struct.unpack_from("<H", chunk, cursor + 6)[0]
        param_bytes = chunk[cursor + 8:cursor + 8 + 2 * param_size]
        if param_type == "int" and len(param_bytes) >= 4:
            value = struct.unpack_from("<i", param_bytes, 0)[0]
        elif param_type == "float" and len(param_bytes) >= 8:
            value = struct.unpack_from("<d", param_bytes, 0)[0]
        elif param_type == "str":
            value = param_bytes.split(b"\x00", 1)[0].decode("latin1", errors="replace")
        else:
            value = param_bytes
        params[name] = value
        cursor += 8 + 2 * param_size
    return params


def _opus_data_param_key(data_type: int, channel_type: int, text_type: int) -> tuple:
    """Pair data blocks with their OPUS parameter block."""
    if data_type == 15:
        return (31, channel_type, text_type)      # AB Data Parameter
    if data_type == 7:
        return (23, channel_type, text_type)      # sample spectrum
    if data_type == 11:
        return (27, channel_type, text_type)      # reference spectrum
    return (-1, -1, -1)


def _opus_block_name(data_type: int, channel_type: int, text_type: int) -> str:
    if data_type == 15:
        return "AB"
    prefix = {7: "ScSm", 11: "ScRf"}.get(data_type, f"Block{data_type}")
    # Bruker channel type 132 appears in Raman OPUS files; keep it readable.
    return prefix if text_type == 0 else f"{prefix} {text_type + 1}"


def _load_opus_native(path: str) -> list[Spectrum]:
    """Small native OPUS reader for processed 1-D spectra, including Raman."""
    with open(path, "rb") as fh:
        raw = fh.read()
    entries = _opus_entries(raw)
    param_blocks = {
        (e["data_type"], e["channel_type"], e["text_type"]): _opus_parse_param(raw, e)
        for e in entries if e["data_type"] in (23, 27, 31)
    }

    spectra = []
    base = os.path.basename(path)
    for entry in entries:
        data_type = entry["data_type"]
        if data_type not in (7, 11, 15):
            continue
        params = param_blocks.get(_opus_data_param_key(
            data_type, entry["channel_type"], entry["text_type"]), {})
        n_meta = int(params.get("NPT", entry["chunk_size"]))
        n = min(max(n_meta, 0), entry["chunk_size"])
        if n <= 0:
            continue
        y = np.frombuffer(raw, "<f4", n, entry["offset"]).astype(float)
        fx, lx = params.get("FXV"), params.get("LXV")
        x = (np.linspace(float(fx), float(lx), n)
             if fx is not None and lx is not None else np.arange(n, dtype=float))
        block = _opus_block_name(data_type, entry["channel_type"], entry["text_type"])
        name = f"{base} ({block})"
        lower_name = base.lower()
        is_raman = "raman" in lower_name
        x_unit = "raman_cm-1" if is_raman else "cm-1"
        y_unit = "intensity" if is_raman else (
            "absorbance" if block == "AB" else "intensity")
        spectra.append(Spectrum(x=x, y=y, name=name, x_unit=x_unit, y_unit=y_unit,
                                meta={"source": path, "format": "Bruker OPUS",
                                      "opus_block": block, **params}))
    if not spectra:
        raise ValueError(f"No recognised spectral block in {base}.")
    # Prefer the processed AB block when present; raw sample/reference blocks
    # remain available in unusual files without AB.
    ab = [s for s in spectra if s.meta.get("opus_block") == "AB"]
    return ab or spectra


def load_opus(path: str) -> list[Spectrum]:
    """Bruker OPUS (.0, .1, ...), with a native fallback for Raman files."""
    try:
        return _load_opus_native(path)
    except Exception as native_exc:
        native_error = native_exc

    try:
        from .ascii_io import load_ascii
        spectra = load_ascii(path)
        if "raman" in os.path.basename(path).lower():
            for spec in spectra:
                spec.x_unit = "raman_cm-1"
                spec.y_unit = "intensity"
                spec.meta["format"] = "ASCII spectrum"
        return spectra
    except Exception:
        pass

    try:
        from brukeropusreader import read_file  # type: ignore
    except ImportError as exc:
        raise MissingDependency(
            "Reading this Bruker OPUS file needs the 'brukeropusreader' package.\n"
            "Install it with:  pip install brukeropusreader"
        ) from exc

    try:
        data = read_file(path)
    except Exception as exc:
        raise ValueError(f"Could not read {os.path.basename(path)} as Bruker OPUS: "
                         f"{native_error}") from exc
    base = os.path.basename(path)
    spectra = []
    # OPUS files contain several blocks; pick the spectral ones we recognise.
    for block in ("AB", "Absorbance", "Transmittance", "ScSm", "ScRf", "Raman"):
        if block not in data:
            continue
        y = np.asarray(data[block], dtype=float)
        try:
            x = np.asarray(data.get_range(block), dtype=float)
        except Exception:
            params = data.get(block + " Data Parameter", {})
            fx, lx = params.get("FXV"), params.get("LXV")
            x = (np.linspace(fx, lx, y.size) if fx is not None and lx is not None
                 else np.arange(y.size, dtype=float))
        n = min(x.size, y.size)
        y_unit = ("absorbance" if block in ("AB", "Absorbance") else
                  "transmittance" if block == "Transmittance" else "intensity")
        spectra.append(Spectrum(x=x[:n], y=y[:n], name=f"{base} ({block})",
                                x_unit="cm-1", y_unit=y_unit, meta={"source": path}))
    if not spectra:
        raise ValueError(f"No recognised spectral block in {base}.")
    return spectra


# ---------------------------------------------------------------------------
# Shimadzu LabSolutions / IRTracer spectrum project data (.ispd)
# ---------------------------------------------------------------------------
# ISPD files are B-tree containers. The FTIR files produced by IRTracer store the
# processed absorbance spectrum as two little-endian float64 pages.
_ISPD_NPOINTS_OFFSET = 15305
_ISPD_FIRST_X_OFFSET = 15325
_ISPD_LAST_X_OFFSET = 15341
_ISPD_STEP_OFFSET = 15357
_ISPD_Y_PAGE1_OFFSET = 49178
_ISPD_Y_PAGE1_POINTS = 1020
_ISPD_Y_PAGE2_OFFSET = 40980


def load_ispd(path: str) -> list[Spectrum]:
    """Shimadzu IRTracer/LabSolutions ``.ispd`` FTIR spectrum project data."""
    with open(path, "rb") as fh:
        raw = fh.read()
    base = os.path.splitext(os.path.basename(path))[0]
    if b"_BTREE_DATA" not in raw or b"IRTracer" not in raw:
        raise ValueError(f"{os.path.basename(path)} is not a recognised ISPD file.")
    if len(raw) < _ISPD_Y_PAGE1_OFFSET + 8:
        raise ValueError(f"{os.path.basename(path)} is too small for ISPD spectrum data.")

    npoints = int(struct.unpack_from("<I", raw, _ISPD_NPOINTS_OFFSET)[0])
    first_x = float(struct.unpack_from("<d", raw, _ISPD_FIRST_X_OFFSET)[0])
    last_x = float(struct.unpack_from("<d", raw, _ISPD_LAST_X_OFFSET)[0])
    step = float(struct.unpack_from("<d", raw, _ISPD_STEP_OFFSET)[0])
    if not (10 <= npoints <= 200000 and np.isfinite(first_x) and np.isfinite(last_x)):
        raise ValueError(f"Invalid ISPD spectral metadata in {os.path.basename(path)}.")

    first_count = min(npoints, _ISPD_Y_PAGE1_POINTS)
    second_count = npoints - first_count
    if _ISPD_Y_PAGE1_OFFSET + 8 * first_count > len(raw):
        raise ValueError(f"Truncated ISPD first data page in {os.path.basename(path)}.")
    y_parts = [np.frombuffer(raw, "<f8", first_count, _ISPD_Y_PAGE1_OFFSET)]
    if second_count:
        if _ISPD_Y_PAGE2_OFFSET + 8 * second_count > len(raw):
            raise ValueError(f"Truncated ISPD second data page in {os.path.basename(path)}.")
        y_parts.append(np.frombuffer(raw, "<f8", second_count, _ISPD_Y_PAGE2_OFFSET))
    y = np.concatenate(y_parts).astype(float)

    x = first_x + step * np.arange(npoints, dtype=float)
    if abs(x[-1] - last_x) > max(abs(step), 1e-9):
        x = np.linspace(first_x, last_x, npoints)

    return [Spectrum(x=x, y=y, name=base, x_unit="cm-1", y_unit="absorbance",
                     meta={"source": path, "format": "Shimadzu ISPD",
                           "npoints": npoints, "first_x": first_x,
                           "last_x": last_x, "step": step})]


# ---------------------------------------------------------------------------
# Bruker handheld XRF (.pdz)  --  pure-Python, no external dependency.
# ---------------------------------------------------------------------------
# Bruker PDZ is a proprietary binary format. The common PDZ25 variant is a
# sequence of blocks, each starting with <block-type:int16, block-size:int32>.
# XRF spectra live in type-3 blocks; their last N int32 values are detector
# counts, and the fixed metadata header contains eV/channel, eV start and N.
_PDZ25_HEADER_FMT = "<hi3i9f7hfhfhfhf8hfhi"
_PDZ25_HEADER_SIZE = struct.calcsize(_PDZ25_HEADER_FMT)


def _pdz_blocks(raw: bytes) -> list[tuple[int, int, int, int]]:
    """Return ``(start, stop, block_type, block_size)`` for PDZ25 blocks."""
    blocks: list[tuple[int, int, int, int]] = []
    start = 0
    n = len(raw)
    while start + 6 <= n:
        block_type, block_size = struct.unpack_from("<hi", raw, start)
        stop = start + block_size + 6
        if block_size < 0 or stop <= start or stop > n:
            raise ValueError("Invalid PDZ block table.")
        blocks.append((start, stop, block_type, block_size))
        start = stop
    if start != n:
        raise ValueError("Trailing bytes after PDZ block table.")
    return blocks


def _pdz_spectrum_from_pdz25_block(raw: bytes, path: str, block: tuple[int, int, int, int],
                                   index: int, total: int) -> Spectrum:
    start, stop, block_type, _block_size = block
    if block_type != 3 or stop - start < _PDZ25_HEADER_SIZE:
        raise ValueError("PDZ block is not a spectral data block.")

    vals = struct.unpack_from(_PDZ25_HEADER_FMT, raw, start)
    ev_per_channel = float(vals[25])
    ev_start = float(vals[27])
    n_channels = int(vals[37])
    xray_voltage_kv = float(vals[12])
    tube_current_ua = float(vals[13])
    live_time_s = float(vals[11])
    active_time_s = float(vals[7])
    raw_counts = int(vals[3])
    valid_counts = int(vals[4])

    available_counts = (stop - start) // 4
    if n_channels <= 0 or n_channels > available_counts:
        for fallback in (2048, 1024):
            if stop - fallback * 4 >= start:
                n_channels = fallback
                break
        else:
            raise ValueError("No PDZ spectral counts found.")

    y_offset = stop - n_channels * 4
    y = np.frombuffer(raw, "<i4", n_channels, y_offset).astype(float)
    if ev_per_channel > 0 and abs(ev_start) < 1e6:
        x = (ev_start + ev_per_channel * np.arange(n_channels, dtype=float)) / 1000.0
        x_unit = "keV"
    else:
        x = np.arange(n_channels, dtype=float)
        x_unit = "pixel"

    base = os.path.splitext(os.path.basename(path))[0]
    name = base if total == 1 else f"{base} [{index + 1}]"
    return Spectrum(
        x=x, y=y, name=name, x_unit=x_unit, y_unit="counts",
        meta={
            "source": path,
            "format": "Bruker PDZ25 XRF",
            "n_channels": n_channels,
            "ev_per_channel": ev_per_channel,
            "ev_start": ev_start,
            "xray_voltage_kv": xray_voltage_kv,
            "tube_current_ua": tube_current_ua,
            "active_time_s": active_time_s,
            "live_time_s": live_time_s,
            "raw_counts": raw_counts,
            "valid_counts": valid_counts,
        },
    )


def _load_pdz25(raw: bytes, path: str) -> list[Spectrum]:
    blocks = _pdz_blocks(raw)
    spectral_blocks = [b for b in blocks if b[2] == 3]
    if not spectral_blocks:
        raise ValueError(f"No spectral data block found in {os.path.basename(path)}.")
    return [_pdz_spectrum_from_pdz25_block(raw, path, b, i, len(spectral_blocks))
            for i, b in enumerate(spectral_blocks)]


def _load_pdz11(raw: bytes, path: str) -> list[Spectrum]:
    """Read older PDZ11 files with 1024 or 2048 int32 channels."""
    if len(raw) < 16:
        raise ValueError(f"{os.path.basename(path)} is too small for a PDZ file.")
    n_channels = int(struct.unpack_from("<h", raw, 6)[0])
    if n_channels not in (1024, 2048):
        raise ValueError(f"Unsupported PDZ11 channel count: {n_channels}.")

    # Legacy files store the counts as a fixed-width int32 array after a fixed
    # metadata header. The energy offset is not reliable, so start at 0 keV.
    header_fmt = ("<2xih34x2d86x2i10x2f188x"
                  if n_channels == 2048
                  else "<2xih34x2d86x2i10x2f24x")
    header_size = struct.calcsize(header_fmt)
    if len(raw) < header_size + n_channels * 4:
        raise ValueError(f"Truncated PDZ11 spectral data in {os.path.basename(path)}.")
    ev_per_channel = float(struct.unpack_from("<d", raw, 42)[0])
    y = np.frombuffer(raw, "<i4", n_channels, header_size).astype(float)
    x = (ev_per_channel * np.arange(n_channels, dtype=float) / 1000.0
         if ev_per_channel > 0 else np.arange(n_channels, dtype=float))
    base = os.path.splitext(os.path.basename(path))[0]
    return [Spectrum(x=x, y=y, name=base, x_unit="keV" if ev_per_channel > 0 else "pixel",
                     y_unit="counts",
                     meta={"source": path, "format": "Bruker PDZ11 XRF",
                           "n_channels": n_channels,
                           "ev_per_channel": ev_per_channel, "ev_start": 0.0})]


def load_pdz(path: str) -> list[Spectrum]:
    """Bruker handheld XRF ``.pdz`` files (PDZ25 and common legacy PDZ11)."""
    with open(path, "rb") as fh:
        raw = fh.read()
    if len(raw) < 8:
        raise ValueError(f"{os.path.basename(path)} is too small for a PDZ file.")

    signature = struct.unpack_from("<h", raw, 0)[0]
    if signature == 25:
        return _load_pdz25(raw, path)
    if signature == 257:
        return _load_pdz11(raw, path)
    raise ValueError(
        f"{os.path.basename(path)} is not a recognised Bruker PDZ XRF file."
    )


# ---------------------------------------------------------------------------
# PerkinElmer Spectrum (.sp)  --  pure-Python, no external dependency.
# ---------------------------------------------------------------------------
# The ".sp" file is a tagged-block binary ("PEPE2D constant interval DataSet
# file"). After a NUL-terminated text signature comes a sequence of blocks,
# each:  <block-id : uint16 LE>  <length : int32 LE>  <content[length]>.
# Every real block-id has high byte 0x8b. The blocks we need:
#   0x8b72  X range   -> member tag (2 B) + 2 float64  (firstX, lastX)
#   0x8b74  delta X   -> member tag (2 B) + 1 float64
#   0x8b75  N points  -> member tag (2 B) + 1 int32
#   0x8b77  X units   -> member tag (2 B) + uint16 len + ASCII
#   0x8b78  Y units   -> member tag (2 B) + uint16 len + ASCII
#   0x8b7b  data type -> member tag (2 B) + uint16 len + ASCII ("Spectrum")
#   0x8b7c  Y data    -> member tag (2 B) + int32 byte-count + float64[N]
_PE_SIG = b"PEPE"

# PerkinElmer unit strings -> SpectraView canonical units.
_PE_XUNITS = {
    "cm-1": "cm-1", "1/cm": "cm-1", "cm^-1": "cm-1", "wavenumber": "cm-1",
    "nm": "nm", "nanometers": "nm", "wavelength": "nm",
    "um": "um", "µm": "um", "micrometers": "um", "micron": "um",
    "ev": "eV", "raman": "raman_cm-1",
}
_PE_YUNITS = {
    "%t": "%T", "%transmittance": "%T", "transmittance": "transmittance",
    "a": "absorbance", "abs": "absorbance", "absorbance": "absorbance",
    "%r": "%R", "reflectance": "reflectance", "log(1/r)": "log1R",
    "kubelka-munk": "KM", "km": "KM", "counts": "counts",
}


def _pe_blocks(raw: bytes) -> dict:
    """Return {block_id: (content_offset, length)} for a PerkinElmer .sp file."""
    start = raw.find(b"\x00")
    i = (start + 1) if start != -1 else 4
    blocks: dict[int, tuple[int, int]] = {}
    n = len(raw)
    while i < n - 6:
        if raw[i + 1] != 0x8B:          # second byte of every block-id is 0x8b
            i += 1
            continue
        bid = struct.unpack_from("<H", raw, i)[0]
        length = struct.unpack_from("<i", raw, i + 2)[0]
        if length < 0 or i + 6 + length > n:
            i += 1                       # spurious 0x8b inside data: resync
            continue
        blocks.setdefault(bid, (i + 6, length))
        i += 6 + length
    return blocks


def load_sp(path: str) -> list[Spectrum]:
    """PerkinElmer Spectrum ``.sp`` (FTIR/UV-Vis). No extra package required."""
    with open(path, "rb") as fh:
        raw = fh.read()
    if not raw.startswith(_PE_SIG):
        raise ValueError(
            f"{os.path.basename(path)} is not a PerkinElmer .sp file "
            "(missing 'PEPE' signature)."
        )
    blocks = _pe_blocks(raw)

    def _doubles(bid, count):
        if bid not in blocks:
            return None
        off, ln = blocks[bid]
        off += 2                         # skip the 2-byte member tag
        if off + 8 * count > len(raw):
            return None
        return struct.unpack_from("<%dd" % count, raw, off)

    def _string(bid):
        if bid not in blocks:
            return ""
        off, ln = blocks[bid]
        slen = struct.unpack_from("<H", raw, off + 2)[0]   # tag(2) + len(2)
        return raw[off + 4:off + 4 + slen].decode("latin1").strip()

    # --- Y data (block 0x8b7c): member tag (2 B) + int32 byte-count + doubles
    if 0x8B7C not in blocks:
        raise ValueError(f"No data block found in {os.path.basename(path)}.")
    doff, _ = blocks[0x8B7C]
    nbytes = struct.unpack_from("<i", raw, doff + 2)[0]
    npts = nbytes // 8
    y = np.frombuffer(raw, "<f8", npts, doff + 6).astype(float)

    # --- X axis from (firstX, lastX, npoints); fall back to delta if needed
    xr = _doubles(0x8B72, 2)
    np_block = _doubles(0x8B75, 1)
    n_meta = int(struct.unpack_from("<i", raw, blocks[0x8B75][0] + 2)[0]) \
        if 0x8B75 in blocks else npts
    if xr is not None:
        first_x, last_x = xr
        x = np.linspace(first_x, last_x, npts)
    else:                                # reconstruct from firstX + deltaX
        dx = _doubles(0x8B74, 1)
        first_x = (xr[0] if xr else 0.0)
        step = dx[0] if dx else 1.0
        x = first_x + step * np.arange(npts)

    x_unit = _PE_XUNITS.get(_string(0x8B77).lower(), "cm-1")
    y_unit = _PE_YUNITS.get(_string(0x8B78).lower(), "transmittance")

    base = os.path.splitext(os.path.basename(path))[0]
    return [Spectrum(x=x, y=y, name=base, x_unit=x_unit, y_unit=y_unit,
                     meta={"source": path, "format": "PerkinElmer .sp",
                           "npoints_meta": n_meta})]
