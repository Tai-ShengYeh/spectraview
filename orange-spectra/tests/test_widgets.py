"""GUI smoke tests for the four Orange widgets (offscreen).

Run:  QT_QPA_PLATFORM=offscreen python -m unittest tests.test_widgets -v
"""
import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Orange.widgets.tests.base import WidgetTest  # noqa: E402

from orangespectra import core  # noqa: E402
from orangespectra.table_io import table_from_spectra  # noqa: E402
from orangespectra.widgets import owimporturl, owlibrary  # noqa: E402
from orangespectra.widgets.owimporturl import OWImportSpectrumURL  # noqa: E402
from orangespectra.widgets.owlibrary import OWSpectralLibrary  # noqa: E402
from orangespectra.widgets.owmixture import OWMixtureAnalysis  # noqa: E402
from orangespectra.widgets.owsimilarity import OWSpectraSimilarity  # noqa: E402
from orangespectra.widgets.owaquagram import OWAquagram  # noqa: E402

X = np.linspace(400, 1800, 200)


def _g(c, w):
    return np.exp(-((X - c) ** 2) / (2.0 * w * w))


def _table(*name_y):
    return table_from_spectra(
        [core.make_spectrum(X, y, name=n, x_label="wavenumber (cm-1)")
         for n, y in name_y])


_IRUG_PAGE = ("<html><script>var jqPlotData={};jqPlotData.series=[{data:["
              + ",".join(f'"{x:.1f}:{y:.4f}"' for x, y in
                         zip(np.linspace(1900, 100, 50), np.cos(np.linspace(0, 6, 50))))
              + "]}];</script></html>").encode()


class TestImportURL(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWImportSpectrumURL)
        self._orig = owimporturl.load_spectrum_url
        owimporturl.load_spectrum_url = lambda u: core.load_spectrum_url(
            u, fetch=lambda _url: (_IRUG_PAGE, "text/html"))

    def tearDown(self):
        owimporturl.load_spectrum_url = self._orig
        super().tearDown()

    def test_fetch_outputs_table(self):
        self.widget.url_edit.setText("4119")
        self.widget.fetch()
        out = self.get_output(self.widget.Outputs.data)
        self.assertIsNotNone(out)
        self.assertEqual(len(out), 1)
        self.assertEqual(len(out.domain.attributes), 50)
        self.assertEqual(str(out.metas[0][0]), "IRUG 4119")
        # second fetch accumulates
        self.widget.fetch()
        out2 = self.get_output(self.widget.Outputs.data)
        self.assertEqual(len(out2), 2)
        self.widget.clear_all()
        self.assertIsNone(self.get_output(self.widget.Outputs.data))

    def test_fetch_error_reported(self):
        owimporturl.load_spectrum_url = lambda u: (_ for _ in ()).throw(
            ValueError("boom"))
        self.widget.url_edit.setText("https://x.org/none")
        self.widget.fetch()
        self.assertTrue(self.widget.Error.fetch_failed.is_shown())


class TestSimilarity(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWSpectraSimilarity)

    def test_data_vs_references(self):
        data = _table(("q", _g(600, 30)))
        refs = _table(("match", _g(600, 30) * 2), ("other", _g(1200, 40)))
        self.send_signal(self.widget.Inputs.data, data)
        self.send_signal(self.widget.Inputs.references, refs)
        out = self.get_output(self.widget.Outputs.scores)
        self.assertEqual(len(out), 2)
        self.assertEqual(str(out.metas[0][1]), "match")   # best first
        self.assertAlmostEqual(float(out.X[0][0]), 1.0, places=6)

    def test_within_data_pairs(self):
        data = _table(("a", _g(600, 30)), ("b", _g(600, 30)), ("c", _g(1200, 40)))
        self.send_signal(self.widget.Inputs.data, data)
        out = self.get_output(self.widget.Outputs.scores)
        self.assertEqual(len(out), 3)                     # C(3,2) pairs

    def test_matrix_within_data(self):
        data = _table(("a", _g(600, 30)), ("b", _g(600, 30)), ("c", _g(1200, 40)))
        self.send_signal(self.widget.Inputs.data, data)
        mat = self.get_output(self.widget.Outputs.matrix)
        self.assertEqual(len(mat), 3)                          # 3x3
        self.assertEqual(len(mat.domain.attributes), 3)
        self.assertEqual([v.name for v in mat.domain.attributes],
                         ["a", "b", "c"])
        for i in range(3):                                     # self-corr = 1
            self.assertAlmostEqual(float(mat.X[i][i]), 1.0, places=6)
        self.assertAlmostEqual(float(mat.X[0][1]), float(mat.X[1][0]),
                               places=10)                      # symmetric

    def test_matrix_data_vs_references(self):
        data = _table(("q", _g(600, 30)))
        refs = _table(("match", _g(600, 30) * 2), ("other", _g(1200, 40)))
        self.send_signal(self.widget.Inputs.data, data)
        self.send_signal(self.widget.Inputs.references, refs)
        mat = self.get_output(self.widget.Outputs.matrix)
        self.assertEqual(len(mat), 1)                          # 1x2
        self.assertEqual(len(mat.domain.attributes), 2)
        self.assertEqual(str(mat.metas[0][0]), "q")
        self.assertAlmostEqual(float(mat.X[0][0]), 1.0, places=6)

    def test_non_spectral_table_errors(self):
        from Orange.data import Table
        self.send_signal(self.widget.Inputs.data, Table("iris"))
        self.assertTrue(self.widget.Error.bad_table.is_shown())


class TestLibrary(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWSpectralLibrary)

    def test_builtin_sugars_library(self):
        # The 9-spectrum sugars NIR library ships inside the package, so the
        # tutorial's demos work without downloading anything from the repo.
        names = [n for n, _ in self.widget._builtin]
        sugars = "Sugars & food additives (NIR, Hadamard)"
        self.assertIn(sugars, names)
        self.widget.builtin_idx = names.index(sugars)
        self.widget.add_builtin()
        lib_out = self.get_output(self.widget.Outputs.library)
        self.assertEqual(len(lib_out), 9)
        self.assertFalse(self.widget.Error.file_error.is_shown())

    def test_build_search_save_load(self):
        refs = _table(("R1", _g(600, 30)), ("R2", _g(1200, 40)))
        self.send_signal(self.widget.Inputs.spectra, refs)
        self.widget.add_pending()
        lib_out = self.get_output(self.widget.Outputs.library)
        self.assertEqual(len(lib_out), 2)

        query = _table(("unknown", _g(600, 30) + 0.01))
        self.send_signal(self.widget.Inputs.query, query)
        hits = self.get_output(self.widget.Outputs.hits)
        self.assertEqual(len(hits), 2)
        self.assertEqual(str(hits.metas[0][1]), "R1")     # best match first
        best = self.get_output(self.widget.Outputs.best)
        self.assertEqual(str(best.metas[0][0]), "R1")

        # save via core (file dialog is interactive) and reload through widget API
        path = os.path.join(tempfile.mkdtemp(), "t.speclib")
        core.save_library(self.widget._library, path)
        self.widget.clear()
        self.assertIsNone(self.get_output(self.widget.Outputs.library))
        loaded = core.load_library(path)
        self.assertEqual([e["name"] for e in loaded], ["R1", "R2"])

    def test_add_builtin_library(self):
        # the UCL entry is offered (download-on-first-use, so not exercised
        # over the network here); widget mechanics are tested with a stub
        names = [n for n, _ in self.widget._builtin]
        self.assertTrue(any(n.startswith("UCL ") for n in names))

        spectra = [core.make_spectrum([1, 2, 3], [1.0, 5.0, 1.0], name="R1"),
                   core.make_spectrum([1, 2, 3], [4.0, 1.0, 4.0], name="R2")]
        calls = []

        def stub_loader(progress=None):
            calls.append(progress)
            if progress:
                progress(1, 1, "done")
            return spectra

        self.widget._builtin = [("Stub library", stub_loader)]
        self.widget.builtin_idx = 0
        self.widget.add_builtin()
        self.assertEqual(len(calls), 1)                 # progress passed in
        self.assertEqual(len(self.widget._library), 2)
        self.assertEqual(self.widget.listing.count(), 2)

        query = table_from_spectra([dict(spectra[0], name="unknown")])
        self.send_signal(self.widget.Inputs.query, query)
        hits = self.get_output(self.widget.Outputs.hits)
        self.assertEqual(len(hits), 2)
        self.assertEqual(str(hits.metas[0][1]), "R1")

    def test_disjoint_ranges_do_not_break_commit(self):
        # regression: UCL's 55 pigments share no common x-range; commit()
        # used to crash in table_from_spectra and never send Hits
        lo = _table(("low", _g(600, 30)))
        xhi = np.linspace(2400, 3000, 200)
        hi = table_from_spectra([core.make_spectrum(
            xhi, np.exp(-((xhi - 2700) ** 2) / (2.0 * 40 * 40)),
            name="high", x_label="wavenumber (cm-1)")])
        self.send_signal(self.widget.Inputs.spectra, lo)
        self.widget.add_pending()
        self.send_signal(self.widget.Inputs.spectra, hi)
        self.widget.add_pending()

        lib_out = self.get_output(self.widget.Outputs.library)
        self.assertIsNotNone(lib_out)                   # union table, no crash
        self.assertEqual(len(lib_out), 2)

        query = _table(("unknown", _g(600, 30) + 0.02))
        self.send_signal(self.widget.Inputs.query, query)
        hits = self.get_output(self.widget.Outputs.hits)
        self.assertIsNotNone(hits)                      # search now runs
        self.assertEqual(len(hits), 2)
        self.assertEqual(str(hits.metas[0][1]), "low")

    def test_add_builtin_error_is_reported(self):
        def boom(progress=None):
            raise OSError("no network")
        self.widget._builtin = [("Broken", boom)]
        self.widget.builtin_idx = 0
        self.widget.add_builtin()
        self.assertTrue(self.widget.Error.file_error.is_shown())
        self.assertEqual(len(self.widget._library), 0)


class TestMixture(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWMixtureAnalysis)

    def test_nnls_composition(self):
        mix = _table(("mix", 3.0 * _g(600, 30) + 1.0 * _g(1200, 40) + 0.05))
        refs = _table(("R1", _g(600, 30)), ("R2", _g(1200, 40)),
                      ("R3", _g(900, 25)))
        self.send_signal(self.widget.Inputs.mixture, mix)
        self.send_signal(self.widget.Inputs.references, refs)
        comp = self.get_output(self.widget.Outputs.composition)
        self.assertEqual(len(comp), 3)
        # sorted by fraction: R1 (3x) first, then R2 (1x)
        self.assertEqual(str(comp.metas[0][0]), "R1")
        self.assertEqual(str(comp.metas[1][0]), "R2")
        self.assertAlmostEqual(float(comp.X[0][0]), 3.0, places=1)
        self.assertGreater(comp.attributes["r_squared"], 0.999)
        fit = self.get_output(self.widget.Outputs.fit)
        self.assertEqual(len(fit), 3)                     # mixture / fit / residual

    def test_no_overlap_error(self):
        mix = _table(("mix", _g(600, 30)))
        x2 = np.linspace(5000, 6000, 50)
        refs = table_from_spectra(
            [core.make_spectrum(x2, np.ones(50), name="far")])
        self.send_signal(self.widget.Inputs.mixture, mix)
        self.send_signal(self.widget.Inputs.references, refs)
        self.assertTrue(self.widget.Error.analysis_failed.is_shown())
        self.assertIsNone(self.get_output(self.widget.Outputs.composition))


class TestAquagram(WidgetTest):
    def setUp(self):
        self.widget = self.create_widget(OWAquagram)

    def _nir_table(self):
        wl = np.linspace(1300, 1600, 200)
        base = 0.5 + 0.1 * np.exp(-((wl - 1440) ** 2) / (2 * 30 ** 2))
        specs = []
        for i in range(3):
            specs.append(core.make_spectrum(wl, base + 0.01 * i, name=f"A{i}",
                                            x_label="wavelength (nm)"))
        for i in range(3):
            y = base + 0.05 * np.exp(-((wl - 1492) ** 2) / (2 * 8 ** 2)) + 0.01 * i
            specs.append(core.make_spectrum(wl, y, name=f"B{i}",
                                            x_label="wavelength (nm)"))
        return table_from_spectra(specs)

    def test_aquagram_output(self):
        self.send_signal(self.widget.Inputs.data, self._nir_table())
        out = self.get_output(self.widget.Outputs.coordinates)
        self.assertIsNotNone(out)
        self.assertEqual(len(out), 6)
        self.assertEqual(len(out.domain.attributes), 12)   # 12 WAMACs
        # default normalization = aquagram -> each band column mean ~0
        self.assertTrue(np.allclose(out.X.mean(axis=0), 0, atol=1e-6))
        self.assertEqual(out.attributes["normalization"], "aquagram")

    def test_normalization_switch(self):
        self.send_signal(self.widget.Inputs.data, self._nir_table())
        self.widget.normalization = 0          # raw
        self.widget._recompute()
        out = self.get_output(self.widget.Outputs.coordinates)
        self.assertFalse(np.allclose(out.X.mean(axis=0), 0, atol=1e-6))

    def test_custom_bands(self):
        self.send_signal(self.widget.Inputs.data, self._nir_table())
        self.widget.bands_edit.setText("1400, 1450, 1500")
        self.widget._bands_changed()
        out = self.get_output(self.widget.Outputs.coordinates)
        self.assertEqual(len(out.domain.attributes), 3)

    def test_out_of_range_warns(self):
        wl = np.linspace(1400, 1450, 60)         # far narrower than the WAMACs span
        t = table_from_spectra([core.make_spectrum(wl, np.ones(60), name="n",
                                                   x_label="wavelength (nm)")])
        self.send_signal(self.widget.Inputs.data, t)
        self.assertTrue(self.widget.Warning.out_of_range.is_shown())

    def test_non_spectral_errors(self):
        from Orange.data import Table
        self.send_signal(self.widget.Inputs.data, Table("iris"))
        self.assertTrue(self.widget.Error.bad_table.is_shown())



class TestOWPeakFinder(WidgetTest):
    def setUp(self):
        from orangespectra.widgets.owpeaks import OWPeakFinder
        self.widget = self.create_widget(OWPeakFinder)

    def _table(self):
        x = np.linspace(400, 1800, 1401)
        y = (np.exp(-((x - 1000) / 15) ** 2)
             + 0.5 * np.exp(-((x - 1450) / 10) ** 2) + 0.01)
        return table_from_spectra([core.make_spectrum(x, y, "s")])

    def test_finds_and_outputs_peaks(self):
        self.send_signal(self.widget.Inputs.data, self._table())
        out = self.get_output(self.widget.Outputs.peaks)
        self.assertEqual(len(out), 2)
        pos = sorted(float(r["position"]) for r in out)
        self.assertAlmostEqual(pos[0], 1000, delta=2)
        self.assertAlmostEqual(pos[1], 1450, delta=2)
        names = [str(v.name) for v in out.domain.attributes]
        self.assertIn("fwhm", names)

    def test_thresholds(self):
        self.send_signal(self.widget.Inputs.data, self._table())
        self.widget.min_height = 80.0
        self.widget._recompute()
        out = self.get_output(self.widget.Outputs.peaks)
        self.assertEqual(len(out), 1)

    def test_clear(self):
        self.send_signal(self.widget.Inputs.data, self._table())
        self.send_signal(self.widget.Inputs.data, None)
        self.assertIsNone(self.get_output(self.widget.Outputs.peaks))

    def test_non_spectral_errors(self):
        from Orange.data import Table
        self.send_signal(self.widget.Inputs.data, Table("iris"))
        self.assertTrue(self.widget.Error.bad_table.is_shown())


class TestOWXRFElementID(WidgetTest):
    def setUp(self):
        from orangespectra.widgets.owxrfid import OWXRFElementID
        self.widget = self.create_widget(OWXRFElementID)

    def _table(self):
        x = np.linspace(1, 20, 1901)
        y = (np.exp(-((x - 6.404) / 0.08) ** 2) * 100
             + np.exp(-((x - 8.048) / 0.08) ** 2) * 60
             + np.exp(-((x - 10.551) / 0.09) ** 2) * 40 + 1.0)
        return table_from_spectra(
            [core.make_spectrum(x, y, "xrf", x_label="energy (keV)")])

    def test_identifies_elements(self):
        self.send_signal(self.widget.Inputs.data, self._table())
        out = self.get_output(self.widget.Outputs.elements)
        syms = {str(r["symbol"]) for r in out}
        self.assertTrue({"Fe", "Cu", "Pb"} <= syms)

    def test_line_filter(self):
        self.send_signal(self.widget.Inputs.data, self._table())
        self.widget.line_filter = 1                      # K lines only
        self.widget._recompute()
        out = self.get_output(self.widget.Outputs.elements)
        self.assertNotIn("Pb", {str(r["symbol"]) for r in out})

    def test_warns_when_not_kev(self):
        x = np.linspace(400, 1800, 500)
        y = np.exp(-((x - 1000) / 15) ** 2) + 0.01
        t = table_from_spectra([core.make_spectrum(x, y, "ir")])
        self.send_signal(self.widget.Inputs.data, t)
        self.assertTrue(self.widget.Warning.not_kev.is_shown())


class TestOWPLSDA(WidgetTest):
    def setUp(self):
        from orangespectra.widgets.owplsda import OWPLSDA
        self.widget = self.create_widget(OWPLSDA)

    def _table(self):
        from Orange.data import ContinuousVariable, DiscreteVariable, Domain, Table
        rs = np.random.RandomState(0)
        grid = np.linspace(0, 1, 50)
        X, y = [], []
        for ci, c in enumerate((0.3, 0.7)):
            for _ in range(8):
                X.append(np.exp(-((grid - c) / 0.08) ** 2) + 0.02 * rs.randn(50))
                y.append(ci)
        dom = Domain([ContinuousVariable(f"{v:.3f}") for v in grid],
                     DiscreteVariable("class", values=("A", "B")))
        return Table.from_numpy(dom, np.array(X), np.array(y, float))

    def test_outputs(self):
        self.send_signal(self.widget.Inputs.data, self._table())
        scores = self.get_output(self.widget.Outputs.scores)
        self.assertEqual(len(scores), 16)
        self.assertEqual(len(scores.domain.attributes), 2)
        vip = self.get_output(self.widget.Outputs.vip)
        self.assertEqual(len(vip), 50)
        vals = [float(r["VIP"]) for r in vip]
        self.assertEqual(vals, sorted(vals, reverse=True))
        pred = self.get_output(self.widget.Outputs.predictions)
        self.assertEqual(str(pred.domain.metas[-1].name), "PLS-DA prediction")
        self.assertIn("100.0%", self.widget.info_label.text())

    def test_no_class_errors(self):
        from Orange.data import ContinuousVariable, Domain, Table
        dom = Domain([ContinuousVariable(f"x{i}") for i in range(5)])
        t = Table.from_numpy(dom, np.random.RandomState(0).rand(6, 5))
        self.send_signal(self.widget.Inputs.data, t)
        self.assertTrue(self.widget.Error.no_class.is_shown())

    def test_components_setting(self):
        self.send_signal(self.widget.Inputs.data, self._table())
        self.widget.n_components = 5
        self.widget._recompute()
        scores = self.get_output(self.widget.Outputs.scores)
        self.assertEqual(len(scores.domain.attributes), 5)



class TestOWMergeSpectra(WidgetTest):
    def setUp(self):
        from orangespectra.widgets.owmerge import OWMergeSpectra
        self.widget = self.create_widget(OWMergeSpectra)

    def _tables(self):
        x1 = np.linspace(400, 1800, 400)
        x2 = np.linspace(600, 2000, 400)
        t1 = table_from_spectra([core.make_spectrum(x1, np.exp(-((x1 - 1000) / 40) ** 2), "A"),
                                 core.make_spectrum(x1, np.exp(-((x1 - 900) / 40) ** 2), "A2")])
        t2 = table_from_spectra([core.make_spectrum(x2, np.exp(-((x2 - 1200) / 50) ** 2), "B")])
        return t1, t2

    def test_merge_two_inputs(self):
        t1, t2 = self._tables()
        self.widget.set_data(t1, 1)
        self.widget.set_data(t2, 2)
        self.widget.handleNewSignals()
        out = self.get_output(self.widget.Outputs.spectra)
        self.assertEqual(len(out), 3)
        attrs = [float(v.name) for v in out.domain.attributes]
        self.assertGreaterEqual(min(attrs), 600)      # overlap region
        self.assertLessEqual(max(attrs), 1800)

    def test_remove_input(self):
        t1, t2 = self._tables()
        self.widget.set_data(t1, 1)
        self.widget.set_data(t2, 2)
        self.widget.handleNewSignals()
        self.widget.set_data(None, 2)
        self.widget.handleNewSignals()
        self.assertEqual(len(self.get_output(self.widget.Outputs.spectra)), 2)

    def test_normalization(self):
        t1, _ = self._tables()
        self.widget.set_data(t1, 1)
        self.widget.normalization = 1                 # max = 1
        self.widget.handleNewSignals()
        out = self.get_output(self.widget.Outputs.spectra)
        self.assertLessEqual(float(np.nanmax(np.abs(out.X))), 1.0 + 1e-6)

    def test_empty(self):
        self.widget.handleNewSignals()
        self.assertIsNone(self.get_output(self.widget.Outputs.spectra))



class TestOWLoadSpectraFiles(WidgetTest):
    def setUp(self):
        from orangespectra.widgets.owloadfiles import OWLoadSpectraFiles
        self.widget = self.create_widget(OWLoadSpectraFiles)

    def _folder(self):
        import tempfile
        d = tempfile.mkdtemp()
        with open(os.path.join(d, "a.csv"), "w") as fh:
            fh.write("wl,i\n400,1\n401,2\n402,4\n")
        with open(os.path.join(d, "wide.csv"), "w") as fh:
            fh.write("name,400,401,402\nS1,1,2,3\nS2,3,2,1\n")
        return d

    def test_folder_source(self):
        self.widget.sources = [self._folder()]
        self.widget._sync_list()
        self.widget._reload()
        out = self.get_output(self.widget.Outputs.spectra)
        self.assertEqual(len(out), 3)

    def test_single_file(self):
        d = self._folder()
        self.widget.sources = [os.path.join(d, "wide.csv")]
        self.widget._sync_list()
        self.widget._reload()
        self.assertEqual(len(self.get_output(self.widget.Outputs.spectra)), 2)

    def test_clear(self):
        self.widget.sources = [self._folder()]
        self.widget._sync_list()
        self.widget._reload()
        self.widget._clear()
        self.assertIsNone(self.get_output(self.widget.Outputs.spectra))
        self.assertTrue(self.widget.Information.nothing.is_shown())



class TestOWSpectrometer(WidgetTest):
    def setUp(self):
        from orangespectra.widgets.owspectrometer import OWSpectrometer
        self.widget = self.create_widget(OWSpectrometer)
        H, W = 40, 640
        img = np.zeros((H, W, 3))
        for col, rgb in [(100, (0.2, 0.4, 1.0)), (500, (1.0, 0.3, 0.1))]:
            band = np.exp(-((np.arange(W) - col) / 4.0) ** 2)
            for c, ch in enumerate(rgb):
                img[:, :, c] += band * ch * 255
        self.widget._rgb = np.clip(img, 0, 255)

    def test_pixel_axis(self):
        self.widget.cal_text = ""
        self.widget._recompute()
        out = self.get_output(self.widget.Outputs.spectrum)
        self.assertEqual(len(out), 1)
        attrs = [float(v.name) for v in out.domain.attributes]
        self.assertEqual(attrs[0], 0)
        self.assertEqual(attrs[-1], 639)

    def test_calibration(self):
        self.widget.cal_text = "100=435.8, 500=611.6"
        self.widget._recompute()
        out = self.get_output(self.widget.Outputs.spectrum)
        attrs = [float(v.name) for v in out.domain.attributes]
        self.assertLess(min(attrs), 440)
        self.assertGreater(max(attrs), 600)
        self.assertIn("R", self.widget.info_label.text())

    def test_bad_calibration_errors(self):
        self.widget.cal_text = "not-a-pair"
        self.widget._recompute()
        self.assertTrue(self.widget.Error.bad_calibration.is_shown())

    def test_near_duplicate_wavelengths_still_make_a_table(self):
        # A quadratic fit near its turning point maps neighbouring pixels to
        # wavelengths equal to 6 significant digits; the output table must
        # still have unique column names instead of crashing in Orange.
        from orangespectra.table_io import unique_axis_names
        x = np.array([500.0, 500.0000001, 500.0000002, 600.0])
        names = unique_axis_names(x)
        self.assertEqual(len(set(names)), 4)
        t = table_from_spectra([core.make_spectrum(x, np.ones(4), name="s")])
        self.assertEqual(len(t.domain.attributes), 4)

    def test_stale_calibration_warns(self):
        w = self.widget
        w.cal_text = ""
        w._recompute()                     # peaks at ~100 and ~500
        w.cal_text = "300=435.8, 620=611.6"  # nowhere near a peak
        w._recompute()
        self.assertTrue(w.Warning.stale_calibration.is_shown())
        w.cal_text = "100=435.8, 500=611.6"
        w._recompute()
        self.assertFalse(w.Warning.stale_calibration.is_shown())

    def test_render_failure_updates_status(self):
        w = self.widget
        w.cal_text = "100=435.8, 500=611.6"
        w._recompute()
        w.table_from_spectra = None
        import orangespectra.widgets.owspectrometer as mod
        orig = mod.table_from_spectra
        mod.table_from_spectra = lambda *a, **k: (_ for _ in ()).throw(
            ValueError("boom"))
        try:
            w._recompute()
        finally:
            mod.table_from_spectra = orig
        self.assertTrue(w.Error.render_failed.is_shown())
        self.assertIn("error", w.info_label.text().lower())
        self.assertIsNone(self.get_output(w.Outputs.spectrum))

    def test_same_wavelength_twice_is_a_clear_error(self):
        # The classic slip: both peaks written with lambda still on 546.1.
        self.widget.cal_text = "100=546.1, 500=546.1"
        self.widget._recompute()
        self.assertTrue(self.widget.Error.bad_calibration.is_shown())
        self.assertIn("546.1", str(self.widget.Error.bad_calibration))
        self.assertFalse(self.widget.Error.render_failed.is_shown())
        self.assertIsNone(self.get_output(self.widget.Outputs.spectrum))

    def test_assign_cursor_refuses_duplicates(self):
        w = self.widget
        w.cal_text = ""
        w._recompute()
        w.cursor_px = 100
        w.assign_nm = "435.8"
        w._assign_cursor()
        self.assertIn("=435.8", w.cal_text)
        # same lambda on another peak -> refused, table unchanged
        before = w.cal_text
        w.cursor_px = 500
        w._assign_cursor()
        self.assertEqual(w.cal_text, before)
        self.assertTrue(w.Error.bad_calibration.is_shown())
        # same peak again with a new lambda -> refused too
        w.cursor_px = 100
        w.assign_nm = "611.6"
        w._assign_cursor()
        self.assertEqual(w.cal_text, before)
        # a genuinely new line is accepted and the fit succeeds
        w.cursor_px = 500
        w._assign_cursor()
        self.assertIn("=611.6", w.cal_text)
        self.assertFalse(w.Error.bad_calibration.is_shown())
        self.assertIsNotNone(self.get_output(w.Outputs.spectrum))

    def test_auto_centre_finds_the_band(self):
        w = self.widget
        img = np.zeros((100, 640, 3))
        img[70:74, :, :] = self.widget._rgb[:4, :, :]     # band at ~72 %
        w._rgb = img
        w.row_center_pct = 10
        w._auto_centre()
        self.assertTrue(68 <= w.row_center_pct <= 76, w.row_center_pct)

    def test_no_image(self):
        w = self.create_widget(OWSpectrometer) if False else self.widget
        w._rgb = None
        w._recompute()
        self.assertIsNone(self.get_output(w.Outputs.spectrum))


if __name__ == "__main__":
    unittest.main()
