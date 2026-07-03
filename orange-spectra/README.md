# orange-spectra — spectroscopy widgets for Orange Data Mining

Five spectroscopy widgets for [Orange Data Mining](https://orangedatamining.com/).
They share their algorithms and the `.speclib` library format with
[SpectraView](https://github.com/Tai-ShengYeh/spectraview), a desktop
spectroscopy viewer. Fetch spectra from public databases by URL, compare and
search them, build reusable reference libraries, decompose mixtures, and draw
aquaphotomics aquagrams — all inside Orange's visual workflow canvas.

| Widget | What it does |
|---|---|
| **Import Spectrum URL** | Paste an **IRUG id/URL** or a **SOPRANO URL** (direct JCAMP-DX/CSV links work too), download and plot the spectrum, and output an Orange `Table` (one spectrum per row, wavenumbers as column names). |
| **Spectra Similarity** | Score similarity between two sets of spectra: correlation, cosine, spectral angle (SAM), and Euclidean distance. |
| **Spectral Library** | Build a reference library, save it as **`.speclib`** (interoperable with SpectraView), and rank an unknown spectrum against the library. |
| **Mixture Analysis** | Decompose a mixed spectrum with non‑negative least squares (NNLS): solve `mixture ≈ Σ cᵢ·refᵢ` and report coefficients, proportions, and R². |
| **Aquagram** | Aquaphotomics: read normalized absorbance at water's 12 characteristic bands (WAMACs) and draw a 12‑axis radar chart (raw / SNV / aquagram normalization). |

The output `Table` follows the
[Orange-Spectroscopy](https://orange-spectroscopy.readthedocs.io/) convention
(column names = wavelength/wavenumber, one spectrum per row), so it plugs
straight into the Spectra viewer widget or into PCA / PLS chemometrics
pipelines.

## Install

> ⚠️ Know **which Orange you run** first — the desktop App and a pip‑installed
> Orange are separate Python environments. Installing into the wrong one means
> the widgets won't appear.

**A. Desktop App** (the standalone program from orangedatamining.com):
`Options ▸ Add-ons… ▸ Add more…`, type **`orange-spectra`**, tick it, **OK**,
and restart.

**B. pip Orange** (started with `python -m Orange.canvas`):

```bash
pip install orange-spectra
python -m Orange.canvas
```

If Orange fails to start with `ImportError: PyQt5 … not available`, it's
missing a Qt binding — install one:

```bash
pip install PyQt5 PyQtWebEngine
```

After (re)starting Orange, a **Spectra** category with the five widgets appears
in the toolbox.

## Quick start

1. Drop in **Import Spectrum URL**, enter `4119` (IRUG's PB15 phthalocyanine
   blue Raman spectrum) → **Fetch**.
2. Fetch a few reference spectra → feed them to the **Spectral Library**
   *Spectra* input → *Add input spectra to library* → *Save…* as a `.speclib`.
3. Feed an unknown spectrum to the library's *Query* input → the *Hits* output
   is the ranked match table.
4. Feed a mixed spectrum to **Mixture Analysis** *Mixture* and the references
   (or the library's *Library* output) to *References* → get component
   proportions and R².

Every widget has an **ℹ How to use** box and a **📖 Open tutorial** button in
its top‑left corner.

## Supported URL formats

- **IRUG** detail pages (jqPlot‑embedded data)
- **SOPRANO** pages (Dygraph‑embedded data)
- **JCAMP-DX** (AFFN plain‑number format)
- Two‑column **CSV/TSV**

Compressed JCAMP (SQZ/DIF) is not parsed here — open it in SpectraView and
export first.

## Documentation

Full tutorial with real‑data demos:
<https://tai-shengyeh.github.io/spectraview/orange.html>
([English](https://tai-shengyeh.github.io/spectraview/orange_en.html)).

## Notes

- You only need to install **once**; updating requires a reinstall
  (`pip install --upgrade orange-spectra`) and an Orange restart.
- Source, issues, and the desktop SpectraView app:
  <https://github.com/Tai-ShengYeh/spectraview>.

## License

MIT.
