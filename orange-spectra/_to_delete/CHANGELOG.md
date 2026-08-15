# Changelog

All notable changes to `orange-spectra` are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Because this package produces numbers that end up in publications, any change that alters
the output of an existing analysis is listed under a **Changed numerics** heading and is
never treated as a patch release.

## [Unreleased]

### Added
- `split_fields()` in `core.py`: data lines are split on comma / tab / semicolon /
  pipe when present, otherwise on runs of whitespace, so space-delimited instrument
  exports (Raman `.txt`, `.dat`, `.asc`, `.xy`) parse like CSV. Blank and comment
  lines (`#`, `%`, `!`, `'`) are skipped. Used by `parse_csv()` and `files.py`
- `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md` and this changelog
- `.gitattributes` normalising line endings (a Windows clone previously showed every
  file as modified)
- `[tool.ruff]` section pinning the lint rule set, and a `dev` extra for
  `pip install -e ".[dev]"`
- `LICENSE` is now included in the sdist and wheel (`license-files`)
- Continuous integration on Linux / macOS / Windows across Python 3.9–3.12
- Cross-implementation validation of `plsda_fit()` against R (`pls::plsr`,
  `method = "oscorespls"`) on a public 10-spectrum gypsum/calcite Raman dataset;
  agreement in VIP, loadings and scores to < 1e-6 (see `tests/test_core.py`)

### Documentation
- Worked example: gypsum vs. calcite Raman discrimination with PCA and PLS-DA,
  reproduced identically in R, Python and the Orange PLS-DA widget

## [0.5.0] — 2026-07-05

### Added
- **PLS-DA** widget — partial least squares discriminant analysis with a class-coloured
  score plot. Outputs Scores, Loadings, VIP and Predictions. Setting: number of components
  (1–20). Implements PLS2 (NIPALS) on one-hot encoded classes with mean-centering
- **XRF Element ID** widget — element identification from X-ray fluorescence spectra
- **Aquagram** widget — aquaphotomics radar charts from water absorption bands

## [0.4.1] — 2026-07-05

### Fixed
- Widget discovery when installed through the Orange Add-ons dialog

## [0.4.0] — 2026-07-05

### Added
- **Peak Finder** widget — peak detection with position, height, FWHM and area
- **Spectral Library** widget — build and rank against `.speclib` reference libraries

## [0.3.0] — 2026-07-04

### Added
- **Mixture Analysis** widget — non-negative least squares (NNLS) unmixing of a measured
  spectrum against a set of reference components
- **Spectra Similarity** widget — correlation, cosine, SAM and Euclidean metrics

## [0.2.0] — 2026-07-04

### Added
- **Merge Spectra** widget — combine spectra from several sources into one table
- **Spectrometer** widget — convert diffraction-grating photographs into calibrated spectra

## [0.1.4] — 2026-07-04

### Added
- **Load Spectra Files** widget — bulk import of JCAMP-DX, CSV and NetCDF

## [0.1.3] — 2026-07-03

### Fixed
- Packaging metadata so the Orange Add-ons dialog lists the package

## [0.1.2] — 2026-07-03

### Added
- First public release
- **Import Spectrum URL** widget — fetch spectra from IRUG, SOPRANO and direct links

[Unreleased]: https://github.com/Tai-ShengYeh/spectraview/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/Tai-ShengYeh/spectraview/releases/tag/v0.5.0
[0.4.1]: https://github.com/Tai-ShengYeh/spectraview/releases/tag/v0.4.1
[0.4.0]: https://github.com/Tai-ShengYeh/spectraview/releases/tag/v0.4.0
[0.3.0]: https://github.com/Tai-ShengYeh/spectraview/releases/tag/v0.3.0
[0.2.0]: https://github.com/Tai-ShengYeh/spectraview/releases/tag/v0.2.0
[0.1.4]: https://github.com/Tai-ShengYeh/spectraview/releases/tag/v0.1.4
[0.1.3]: https://github.com/Tai-ShengYeh/spectraview/releases/tag/v0.1.3
[0.1.2]: https://github.com/Tai-ShengYeh/spectraview/releases/tag/v0.1.2
