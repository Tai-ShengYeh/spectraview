# Aquagram validation on real NIR spectra

This folder validates the **Aquagram** widget core (`orangespectra.core.aquagram_coordinates`)
on real near-infrared spectra, as requested in `HANDOFF.md` §6 item 5
(previously the widget was only checked against synthetic spectra).

## Data

Public **Eigenvector corn dataset** — 80 corn samples measured on the `m5`
spectrometer, 1100–2498 nm at 2 nm (700 channels), with reference moisture /
oil / protein / starch. Source:
<https://eigenvector.com/resources/data-sets/>. The range fully covers the 12
WAMACs water bands (1342–1516 nm).

## What is checked

Run:

```bash
pip install numpy scipy matplotlib
python orange-spectra/validation/validate_aquagram_real_nir.py
```

1. **Coverage** — all 12 WAMACs fall inside the corn 1100–2498 nm range (`covered=True`).
2. **Math invariants** — for `normalization="aquagram"` the n×12 coordinate
   columns are zero-mean and unit-std across the sample set (classic normalized
   aquagram, 0 = group average).
3. **Physical meaning** — samples are split into low- vs high-moisture groups
   (means 9.74 % vs 10.72 %). The aquagram separates them at the water bands:
   the high-moisture group sits **above** average at the short WAMACs
   (1342–1382 nm) and **below** average at 1492–1516 nm — consistent with
   moisture-driven O–H overtone/hydrogen-bonding changes.

## Figures

- `aquagram_real_corn_nir.png` — the classic normalized aquagram (left) vs
  SNV-only sampled at the WAMACs (right), low vs high moisture.
- `aquagram_wamacs_coverage.png` — group-mean corn spectra with the 12 WAMACs
  marked, showing every band is inside the measured range.

## Result

`VALIDATION PASSED`: the widget produces correct math and physically sensible
aquagrams on real NIR data.
