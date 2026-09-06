---
title: 'orange-spectra: visual-programming widgets for heritage and materials spectroscopy in Orange'
tags:
  - Python
  - Orange
  - spectroscopy
  - Raman
  - FTIR
  - XRF
  - chemometrics
  - cultural heritage
  - conservation science
  - near-infrared
  - food analysis
authors:
  - name: Tai-Sheng Yeh
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: PLEASE FILL IN — institution, city, country
    index: 1
date: DD Month YYYY
bibliography: paper.bib
---

# Summary

`orange-spectra` is an add-on for the Orange Data Mining platform [@Demsar2013] that brings
spectroscopic workflows to conservation scientists, archaeometrists, materials
researchers and food-analysis teaching laboratories, without programming. It provides
eleven widgets covering the steps that sit
between a spectrometer and a conclusion: importing spectra from public heritage databases
and instrument files, building and searching reference libraries, detecting and
characterising peaks, unmixing composite spectra, identifying elements in X-ray
fluorescence data, and performing partial least squares discriminant analysis (PLS-DA) with
variable importance in projection (VIP) scores.

Every widget is a node on Orange's visual canvas, so an analysis is assembled by connecting
boxes rather than by writing code, and the resulting workflow file is itself a
human-readable, shareable record of what was done.

The numerical layer (`orangespectra.core`) is a dependency-light NumPy/SciPy module with no
Orange or Qt imports. It can therefore be imported, tested and cited independently of the
graphical interface, and the same functions that the widgets call can be used inside a
script or notebook.

# Statement of need

Heritage and materials laboratories routinely combine Raman, FTIR and XRF measurements to
identify pigments, binders, stones and corrosion products. The analytical steps involved —
baseline correction, library matching, spectral unmixing, discriminant analysis — are
individually standard, but the software situation is not. Instrument vendors ship closed
packages tied to one file format and one machine. General-purpose chemometrics libraries
require programming. Domain databases such as IRUG publish reference spectra that must be
downloaded and reformatted by hand before they can be compared with a measurement.

Orange, with its Spectroscopy add-on [@Toplak2021], has become an effective answer to the
programming barrier in biospectroscopy: it offers preprocessing, machine learning and
visualisation on a drag-and-drop canvas, and its workflows double as documentation.
`orange-spectra` extends the same approach to the heritage and materials domain by adding
capabilities the existing ecosystem lacks:

- **Direct import from heritage spectral databases.** The *Import Spectrum URL* widget
  retrieves spectra from IRUG and SOPRANO by URL, removing a manual reformatting step that
  is both tedious and a common source of transcription error.
- **Reference libraries with ranked matching.** *Spectral Library* builds portable
  `.speclib` files and ranks an unknown against them using several similarity metrics,
  which is the everyday task in pigment and mineral identification.
- **Mixture analysis.** Heritage samples are almost always mixtures. *Mixture Analysis*
  performs non-negative least squares unmixing of a measurement against chosen reference
  components.
- **PLS-DA with VIP.** Discriminant analysis is the standard tool for comparing groups of
  spectra, but interpreting *which wavenumbers* drive a separation requires VIP scores,
  which are not available in a point-and-click form elsewhere in the Orange ecosystem.
- **Element identification for XRF**, and **aquaphotomics** radar charts, extending coverage
  beyond vibrational spectroscopy.

The same four steps — build a reference library, rank an unknown against it, unmix a
mixture, and ask a discriminant model *which bands* separate two groups — recur outside
heritage science. Food authentication with hand-held near-infrared (NIR) spectrometers is
one such case, and one where the users are often undergraduates in a laboratory course
rather than chemometricians. The second example below is drawn from that setting.

# Implementation and validation

Widgets follow Orange's input/output signal conventions. The *PLS-DA* widget accepts a data
table with a categorical target, one-hot encodes the classes, fits PLS2 by the NIPALS
algorithm with mean-centering, and emits four tables: component scores, X-loadings, VIP
scores sorted in descending order, and per-sample predictions. Because loading and VIP
tables carry the original variable names in a metadata column, results can be traced
directly back to wavenumbers.

VIP is computed as

$$\mathrm{VIP}_j = \sqrt{\frac{p \sum_a \mathrm{SSY}_a \left( w_{ja} / \lVert \mathbf{w}_a \rVert \right)^2}{\sum_a \mathrm{SSY}_a}}$$

where $\mathrm{SSY}_a = (\mathbf{t}_a^\top \mathbf{t}_a)(\mathbf{q}_a^\top \mathbf{q}_a)$ is
the response variance explained by component $a$ [@Chong2005].

Software that produces numbers for publication should be checked against an independent
implementation. We validated `plsda_fit()` against the R package `pls`
[@Mevik2007] using `method = "oscorespls"` on a ten-spectrum Raman dataset of gypsum
(CaSO$_4\cdot$2H$_2$O) and calcite marble (CaCO$_3$). Scores, X-loadings and VIP agree to
better than $10^{-6}$ in absolute value across all 664 wavenumbers, the residual being
attributable to decimal truncation on export. The comparison is included in the test suite
and, together with the dataset, in the accompanying tutorial.

# Example 1: interpreting a discriminant model

TODO — one short paragraph plus one figure. Suggested content, drawn from the accompanying
tutorial:

The gypsum/calcite dataset separates perfectly (leave-one-out accuracy 1.000; exact
permutation test over all $\binom{10}{5} = 252$ label assignments, $p = 2/252 = 0.0079$).
The VIP output ranks the calcite $\nu_1$ band at 1086 cm$^{-1}$ first and the gypsum
$\nu_1$ band at 1008 cm$^{-1}$ eleventh, both chemically expected. It also ranks a feature
at 419 cm$^{-1}$ second — a band present in every spectrum, including the calcite samples,
which have no vibrational mode there. Inspecting the loadings confirms that the four
largest contributions to the second latent variable all lie in this band. The feature is
instrumental, not chemical.

This is the argument for exposing loadings and VIP in a point-and-click tool rather than
only reporting accuracy: a model can be statistically significant and still be leaning on
an artefact, and in vibrational spectroscopy the artefact is detectable because every band
can be checked against a known assignment.

*(Figure to add: VIP plot with the 419 cm$^{-1}$ band marked, or the two-panel loading
plot.)*

# Example 2: sugars and sweeteners by hand-held NIR in a food-analysis course

The second example uses data that ship with the package. Nine food ingredients — five
sugars (fructose, glucose, lactose, maltose, sucrose) and four additives (aspartame,
benzoic acid, caffeine, sucralose) — were each scanned five times as powders with a
hand-held DLP Hadamard-scan NIR spectrometer (InnoSpectra NIR-S-R14; 1600–2400 nm, 200
points, 12 nm pattern width). The 45 absorbance spectra are provided as a single table
(`examples/sugars_nir_replicates.csv`, column names = wavelengths, one row per scan), and
the per-substance means as the built-in *Sugars & food additives (NIR)* library that the
*Spectral Library* widget loads with one click. The whole example runs on Orange's canvas;
the script `examples/sugars_nir_plsda.py` reproduces every number below and
\autoref{fig:sugars} from the same `core` functions the widgets call.

**Identification.** With the library rebuilt from the other four scans of every substance,
each held-out scan was ranked by Pearson correlation. All 45 were identified correctly
(top-1 correlation $\geq$ 0.9999). The margin to the runner-up was smallest for sucrose,
whose nearest neighbour is fructose (0.986), and the full 9 $\times$ 9 correlation matrix
tells a chemically sensible story: the sugars correlate at 0.83–0.99 with one another,
while caffeine correlates with the sugars at only 0.71–0.74 and benzoic acid at
0.76–0.85. *Mixture Analysis* (non-negative least squares [@Lawson1974]) on the same
held-out scans assigned 97.6–100 % (median 99.8 %) of each to the correct reference with
$R^2 \geq$ 0.9999; the only spill-over is between the sugars themselves, never into an
additive. That is the sanity check a student should run before trusting the widget on a
real blend.

**Discrimination and interpretation.** A *PLS-DA* model of sugar versus additive reaches
leave-one-out accuracy 40/45 with two components and 45/45 with three; the nine-class model
needs eight components for 45/45. The first component carries 97 % of the spectral variance
but is essentially the overall absorbance level — a scatter and particle-size effect visible
in \autoref{fig:sugars}a — so the class information lives in components that explain about
1 % of $X$. This is a useful classroom demonstration that explained $X$ variance is not
discriminative power, and it motivates scatter correction (SNV, available from the
Orange-Spectroscopy *Preprocess* widget) as the next step.

The VIP output of the three-component model has three bands above 1
(\autoref{fig:sugars}b): 1962–2126 nm, 2221–2271 nm and 1600–1640 nm. The first two are
the O–H stretch/deformation and C–H stretch/deformation combination bands of carbohydrates
[@Workman2012], which is what a sugar-versus-non-sugar model *should* use. The third sits
at the short-wavelength edge of the scan, exactly where an instrumental artefact would
appear, and it is the case Example 1 warns about. Here the widget's tables settle the
question in the student's favour: the sugars absorb 0.72 on average in that band against
0.37 for the additives, while the replicate standard deviation there is 0.0007. The band is
the long-wavelength tail of the hydroxyl first overtone, real and reproducible, not noise.
Having the loadings, VIP and the original spectra side by side on the canvas is what makes
that check a two-minute exercise instead of a script.

![(a) Mean NIR absorbance spectra of the nine ingredients, five replicate scans each.
(b) VIP scores of the three-component PLS-DA sugar-versus-additive model; shaded bands are
VIP $>$ 1 regions wider than 20 nm.\label{fig:sugars}](fig_sugars_nir.png)

# Availability

`orange-spectra` is available on PyPI (`pip install orange-spectra`) and through Orange's
Add-ons dialog. Source, issue tracker and documentation are at
<https://github.com/Tai-ShengYeh/spectraview>; releases are archived on Zenodo
[@Yeh2026]. The package is released under the MIT licence and requires Python $\geq$ 3.9
and Orange $\geq$ 3.34. The NIR dataset of Example 2 is distributed with the repository
under the same licence; the UCL pigment library is fetched from its publisher on first use
and is not redistributed.

# Acknowledgements

TODO — funding, institutions, colleagues who tested the widgets, sources of reference data.

# References
