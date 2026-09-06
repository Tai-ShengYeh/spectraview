# Changelog

All notable changes to orange-spectra are documented here.

Format: [Keep a Changelog](https://keepachangelog.com/), versioning: [SemVer](https://semver.org/).

## [0.7.3] - 2026-09-06

### Added

- Spectrometer: **Auto-centre strip** button — moves the ROI strip onto the
  brightest row of the (rotated) photo, so a thin spectrum band in a mostly
  black hand-held grating photo is found without dragging the slider.
- Spectrometer: the Channel combo explains when to prefer `blue` /
  `sum(RGB)` (violet/blue lines such as Hg 405 / 436 nm are nearly invisible
  in `luminance`, which weights blue at only 11 %).

### Fixed

- Spectrometer: writing every peak into the calibration table with `λ =`
  still on the previous value produced a flat fit and the bare
  "calibration error" placeholder. The widget now refuses to add a peak whose
  wavelength (or pixel) is already in the table, `fit_calibration` rejects
  identical wavelengths / identical pixels with a specific message, and the
  placeholder plot echoes the actual reason.

## [0.7.2] - 2026-08-25

### Added

- Spectral Library: the 9-spectrum **sugars & food additives NIR library**
  (`sugars_nir.speclib`) now ships inside the package and appears in the
  "Add built-in" dropdown, so the tutorial's demos work offline without
  downloading anything from the repo. (The UCL pigment library remains
  download-on-first-use for licensing reasons.)

## [0.7.1] - 2026-08-25

### Added

- Spectra Similarity: new **Similarity Matrix** output — a wide
  queries × references matrix of the "Sort scores by" metric (Data-only
  mode gives the full n×n symmetric matrix, diagonal = 1). Feed it to
  Orange's native Heat Map widget to reproduce the tutorial's 9×9
  correlation heatmap entirely inside Orange.

## [0.7.0] - 2026-08-25

### Added

- Spectrometer: image rotation (90/180/270° clockwise) — photos taken with a
  vertical dispersion axis can now use the horizontal ROI. Rotation is
  lossless (an axis swap, no interpolation).
- Spectrometer: interactive peak cursor — automatic peak detection with
  prominence / min-distance / smoothing controls, prev/next-peak jumping,
  sub-pixel snap to the nearest peak, a pixel/intensity/wavelength readout,
  and one-click "write into the calibration table" for the line under the
  cursor.

### Fixed

- Spectrometer: the cursor slider's value label was sized before any image
  was loaded, so 3-4 digit pixel positions rendered clipped (e.g. "335"
  showed as "35").
- Spectrometer: calibration input is validated with specific error messages
  (non-numeric pairs, inf/nan, wavelengths outside 10-100000 nm, pixels far
  outside the image, fewer than 2 distinct pixel positions).
- Spectrometer: a calibration that cannot be drawn (non-finite or degenerate
  axis range) reports an error instead of crashing inside Qt's paint loop.

## [0.6.0] - 2026-08-16

### Added

- Spectral Library: "Add built-in" — the UCL Raman Library of Pigments
  (Bell, Clark & Gibbs 1997, 55 pigments) is downloaded from UCL's own
  website on first use and cached locally (`~/.orange-spectra/`). The data
  itself is deliberately NOT redistributed with this package.
- `core.read_spc()` — Thermo GRAMS .spc parser (new-format LSB, float32 and
  integer-exponent y, TXVALS, descending x).
- `core.fetch_ucl_library()` / `core.build_ucl_library_from_folder()` /
  `core.available_libraries()`; `.speclib` files now round-trip the
  optional `color` / `source` fields.
- `core.merge_spectra_union()` and `table_from_spectra(union=True)` —
  shared grid over the union of x-ranges, NaN where not measured.
- CJK-capable preview fonts (`orangespectra.mplfonts`): spectrum names in
  Chinese/Japanese/Korean no longer render as boxes in matplotlib previews.
- Library download falls back to the Internet Archive's copy when UCL's
  server is unreachable or blocks the user's region (403), and honours the
  `ORANGE_SPECTRA_UCL_BASE` mirror override.

### Fixed

- Spectral Library: `commit()` crashed when library entries had no common
  x-range (true for the whole UCL library), leaving the status at
  "Library: 0 entries" and never sending Hits/Best/Library outputs.
