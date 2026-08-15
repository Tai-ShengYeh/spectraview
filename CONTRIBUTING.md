# Contributing to orange-spectra

Thanks for taking the time to contribute. This project provides Orange Data Mining
widgets for spectroscopy, with a focus on cultural-heritage and materials analysis.

Contributions of any size are welcome — bug reports, documentation fixes, new
reference data, and new widgets.

## Ways to contribute

| I want to… | Do this |
|---|---|
| Report a bug | [Open an issue](https://github.com/Tai-ShengYeh/spectraview/issues/new) with the bug template |
| Request a feature | Open an issue describing the analysis you cannot currently do |
| Fix a typo or improve docs | Send a pull request directly — no issue needed |
| Add a widget | **Open an issue first** so we can discuss scope and the input/output signals |
| Report a wrong number | Open an issue and include the input file plus the expected value and its source |

## Reporting a bug

Please include:

1. Orange version (`Help → About`) and `orange-spectra` version (`pip show orange-spectra`)
2. Operating system and Python version
3. The widget involved and its settings
4. A **minimal** input file that reproduces the problem, if you can share one
5. What you expected versus what happened; paste the full traceback if there was an error

Numerical bugs are the highest priority. If a widget produces a value you believe is
wrong, say what the correct value should be and where that reference comes from
(a textbook formula, another package, a published dataset).

## Development setup

```bash
git clone https://github.com/Tai-ShengYeh/spectraview.git
cd spectraview/orange-spectra
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
```

Run the test suite. The two files are written differently, so run them separately:

```bash
# Numerical core -- a self-contained script. Prints "N passed, M failed"
# and exits 1 on failure. No Orange or Qt needed.
python tests/test_core.py

# Widgets -- unittest based, needs a Qt platform.
QT_QPA_PLATFORM=offscreen python -m unittest tests.test_widgets -v
```

On a headless Linux machine the widget tests need a virtual display:

```bash
QT_QPA_PLATFORM=offscreen xvfb-run -a python -m unittest tests.test_widgets -v
```

> **Do not run `pytest tests/`.** `tests/test_core.py` calls `sys.exit()` at module
> scope, which makes pytest abort during collection with an `INTERNALERROR`.
> Converting these to pytest would be a welcome contribution -- see the issue tracker.

## Pull requests

1. Branch from `main`: `git checkout -b fix/short-description`
2. Keep each pull request focused on one thing
3. **Add a test.** Every bug fix should come with a test that fails before the fix
4. Run `pytest tests/` before pushing
5. Update `CHANGELOG.md` under `## [Unreleased]`
6. Describe *why* the change is needed, not only what it does

## Code style

- Follow PEP 8; keep lines under 100 characters
- Every public function needs a docstring stating what it returns
- **Numerical code must cite its source.** If you implement a published algorithm, put the
  formula in a comment with the reference. See the VIP formula in `core.py` for the pattern
- Keep `core.py` free of Orange and Qt imports — it is the layer other people can reuse
- New widgets go in `orangespectra/widgets/ow<name>.py` and need an icon in `widgets/icons/`

## Widget conventions

New widgets should follow the existing pattern:

```python
class OWMyWidget(OWWidget):
    name = "My Widget"
    description = "One sentence, shown in the widget tooltip."
    icon = "icons/mywidget.svg"
    priority = 90                      # controls position in the toolbar

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        result = Output("Result", Table, default=True)

    class Error(OWWidget.Error):
        bad_input = Msg("Explain what the user should do about it.")
```

Guidelines:

- Error messages should tell the user **what to do**, not just what went wrong.
  Compare `"Invalid data"` with
  `"Data needs a categorical target (class) variable — use Select Columns to set one."`
- Put the maths in `core.py` as a plain function taking and returning numpy arrays.
  The widget should only handle input/output conversion and drawing
- Output tables should carry the original variable names in a `variable` meta column so
  results can be traced back to wavenumbers

## Adding reference data

Reference spectra and element line tables must state their source and licence.
Do not add data you are not permitted to redistribute. If a database allows programmatic
access but not redistribution, add an importer rather than a copy of the data.

## Scientific correctness

This package is used to produce numbers that end up in publications. Two rules follow:

1. **Validate against an independent implementation.** When you add or change an algorithm,
   check the output against a reference implementation in another package (R, MATLAB,
   scikit-learn) on a dataset you can share, and put that comparison in the tests
2. **Do not silently change results.** If a fix changes the numbers a previous version
   produced, say so prominently in `CHANGELOG.md` under a `### Changed numerics` heading

## Releasing

Maintainers only:

1. Update `CHANGELOG.md` — move `[Unreleased]` items under a new version heading with a date
2. Bump `version` in `pyproject.toml`
3. `git tag -a v0.x.0 -m "v0.x.0"` and push the tag
4. The release workflow builds and uploads to PyPI, and Zenodo mints a DOI for the tag

## Code of conduct

This project follows the [Contributor Covenant](CODE_OF_CONDUCT.md). By participating you
agree to abide by it.

## Licence

Contributions are accepted under the project's [MIT licence](LICENSE).
