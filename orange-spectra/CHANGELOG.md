\# Changelog



All notable changes to orange-spectra are documented here.

Format: \[Keep a Changelog](https://keepachangelog.com/), versioning: \[SemVer](https://semver.org/).



\## \[0.6.0] - 2026-08-16



\### Added

\- Spectral Library: "Add built-in" — the UCL Raman Library of Pigments

&#x20; (Bell, Clark \& Gibbs 1997, 55 pigments) is downloaded from UCL's own

&#x20; website on first use and cached locally (`\~/.orange-spectra/`). The data

&#x20; itself is deliberately NOT redistributed with this package.

\- `core.read\_spc()` — Thermo GRAMS .spc parser (new-format LSB, float32 and

&#x20; integer-exponent y, TXVALS, descending x).

\- `core.fetch\_ucl\_library()` / `core.build\_ucl\_library\_from\_folder()` /

&#x20; `core.available\_libraries()`; `.speclib` files now round-trip the

&#x20; optional `color` / `source` fields.

\- `core.merge\_spectra\_union()` and `table\_from\_spectra(union=True)` —

&#x20; shared grid over the union of x-ranges, NaN where not measured.

\- CJK-capable preview fonts (`orangespectra.mplfonts`): spectrum names in

&#x20; Chinese/Japanese/Korean no longer render as boxes in matplotlib previews.

\- Library download falls back to the Internet Archive's copy when UCL's

&#x20; server is unreachable or blocks the user's region (403), and honours the

&#x20; `ORANGE\_SPECTRA\_UCL\_BASE` mirror override.



\### Fixed

\- Spectral Library: `commit()` crashed when library entries had no common

&#x20; x-range (true for the whole UCL library), leaving the status at

&#x20; "Library: 0 entries" and never sending Hits/Best/Library outputs.

