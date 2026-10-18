"""mu2ecombine tests on synthetic inputs: python3 -m unittest discover -s tests (from combine/)."""

import itertools
import os
import shutil
import sys
import tempfile
import unittest
import warnings

import numpy as np
import ROOT

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import mu2ecombine as mc  # noqa: E402
from mu2ecombine.hist import arrays, integral  # noqa: E402
from mu2ecombine.model import matches  # noqa: E402

_ids = itertools.count()


def toy(x, expected, nbins=300, lo=90., hi=120.):
    """A TH1D of the values x, normalized to `expected` events."""
    h = ROOT.TH1D(f"toy{next(_ids)}", "", nbins, lo, hi)
    h.SetDirectory(0)
    h.Sumw2()
    x = np.ascontiguousarray(x, dtype=float)
    h.FillN(len(x), x, np.full(len(x), expected / len(x)))
    return h


def scaled(h, k):
    out = h.Clone(f"{h.GetName()}_x{next(_ids)}")
    out.SetDirectory(0)
    out.Scale(k)
    return out


class HistTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(1)
        self.h = toy(self.rng.normal(105, 1, 10000), 10.)

    def test_rebin(self):
        r = mc.rebin_to(self.h, np.linspace(100, 110, 11))
        self.assertAlmostEqual(integral(r), self.h.Integral(self.h.FindBin(100.01), self.h.FindBin(109.99)), places=9)
        self.assertAlmostEqual(r.GetBinError(5) ** 2, sum(self.h.GetBinError(i) ** 2 for i in range(141, 151)),
                               places=12)
        with self.assertRaises(ValueError):
            mc.rebin_to(self.h, [100.05, 101])

    def test_transformed_keeps_contents(self):
        t = mc.transformed(self.h, shift=0.3)
        self.assertAlmostEqual(integral(t), integral(self.h), places=6)
        self.assertGreater(t.GetMean(), self.h.GetMean())
        self.assertFalse(t.GetDirectory())

    def test_replace_tail_and_convolution(self):
        flat = toy(np.linspace(90.05, 119.95, 300), 300.)
        c = mc.convolved(flat, toy(self.rng.normal(0, 0.5, 5000), 1., 60, -3., 3.))
        self.assertAlmostEqual(integral(c), integral(flat), places=6)
        out = mc.replace_tail(self.h, 107., flat, rel_error=0.2)
        e, v, _ = arrays(out)
        _, v0, _ = arrays(self.h)
        above = 0.5 * (e[1:] + e[:-1]) > 107.
        self.assertAlmostEqual(v[above].sum(), v0[above].sum(), places=9)
        self.assertAlmostEqual(v[~above].sum(), v0[~above].sum(), places=9)

    def test_matches(self):
        self.assertTrue(matches("rmc_ext", "rmc_*"))
        self.assertFalse(matches("cosmic", "!cosmic"))
        self.assertTrue(matches("dio", ["*", "!cosmic"]))
        self.assertTrue(matches("x", None))


class ShapeTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(2)
        self.sr = toy(rng.uniform(90, 120, 200), 4.)
        self.cr = toy(rng.uniform(90, 120, 5000), 50.)
        self.edges = np.linspace(100, 110, 21)

    def test_shape_from_control_region(self):
        t = mc.Template(self.sr, shape=[self.sr, self.cr])
        nominal = arrays(t.nominal(self.edges))[1]
        self.assertAlmostEqual(nominal.sum(), integral(mc.rebin_to(self.sr, self.edges)), places=9)
        expected = arrays(mc.rebin_to(mc.total([self.sr, self.cr]), self.edges))[1]
        np.testing.assert_allclose(nominal / nominal.sum(), expected / expected.sum())

    def test_variation_modes(self):
        up, down = scaled(self.sr, 1.2), scaled(mc.transformed(self.sr, shift=0.5), 0.9)
        for mode, k_up in (("shape", 1.2), ("shape_only", 1.0), ("rate", 1.2)):
            t = mc.Template(self.sr, variations={"v": (up, down, mode)})
            u, d = t.varied("v", self.edges)
            self.assertAlmostEqual(integral(u) / t.rate(self.edges), k_up, places=6, msg=mode)
            if mode == "rate":
                dv, nv = arrays(d)[1], arrays(t.nominal(self.edges))[1]
                np.testing.assert_allclose(dv / dv.sum(), nv / nv.sum())

    def test_bin_fractions_are_integrals(self):
        from mu2ecombine.pdfs import bin_fractions
        x = mc.observable("xf", 0, 1, bins=4)
        slope = ROOT.RooRealVar("slope_f", "", -2.)
        pdf = ROOT.RooExponential("expo_f", "", x, slope)
        cdf = lambda a: (1 - np.exp(-2 * a)) / (1 - np.exp(-2))
        expected = np.diff([cdf(a) for a in (0, 0.25, 0.5, 0.75, 1)])
        np.testing.assert_allclose(bin_fractions(pdf, x, np.linspace(0, 1, 5)), expected, rtol=1e-6)


class ModelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(3)
        cls.sig = toy(rng.normal(105, 0.4, 20000), 2.)
        cls.bkg = toy(90 + rng.exponential(8, 20000), 50.)
        cls.tmp = tempfile.mkdtemp(prefix="m2c_test_")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def build(self, name, fit=False):
        model = mc.Model(name, poi_range=(0, 10))
        ch = model.channel("sr", mc.observable(f"p_{name}", 100, 110, bin_width=0.5))
        ch.add("signal", mc.Fit("gauss", self.sig) if fit else self.sig, signal=True)
        ch.add("bkg", mc.Template(self.bkg, variations={"bscale": (scaled(self.bkg, 1.1), scaled(self.bkg, 0.95))}))
        model.lnN("lumi", (1.1, 0.95), processes="!bkg")
        model.shape_shift("pscale", shift=0.1)
        model.rate_param("k_bkg", "bkg")
        out = os.path.join(self.tmp, name)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model.write(out)
        with open(os.path.join(out, "combined.txt")) as f:
            return model, out, f.read()

    def test_binned_card(self):
        _, out, card = self.build("binned")
        self.assertIn("shapes * sr shapes/sr.root $PROCESS $PROCESS_$SYSTEMATIC", card)
        self.assertIn("lumi lnN", card)
        self.assertIn("0.9500/1.1000", card)
        self.assertIn("bscale shape", card)
        self.assertIn("pscale shape", card)
        self.assertIn("k_bkg rateParam sr bkg 1 [0,10]", card)
        self.assertNotIn("autoMCStats", card)
        self.assertIn("mcstat_sr_bkg lnN", card)
        f = ROOT.TFile(os.path.join(out, "shapes", "sr.root"))
        self.assertAlmostEqual(f.Get("bkg").Integral(), integral(mc.rebin_to(self.bkg, np.linspace(100, 110, 21))),
                               places=6)
        self.assertTrue(f.Get("signal_pscaleUp"))
        f.Close()

    def test_workspace_channel_links_nuisances(self):
        _, out, card = self.build("workspace", fit=True)
        self.assertIn("shapes signal sr shapes/sr.root w:sr_signal", card)
        self.assertIn("pscale shape", card)       # the template's shape line declares the shared nuisance
        self.assertNotIn("pscale param", card)
        self.assertIn("mcstat_sr_bkg lnN", card)
        f = ROOT.TFile(os.path.join(out, "shapes", "sr.root"))
        w = f.Get("w")
        mean = w.function("sr_signal_mean_eff")
        theta = w.var("pscale")
        nominal = mean.getVal()
        theta.setVal(1.)
        self.assertAlmostEqual(mean.getVal() - nominal, 0.1, delta=0.02)   # refit to the shifted histogram
        f.Close()

    def test_param_nuisance_and_user_pdf(self):
        model = mc.Model("userpdf")
        x = mc.observable("p_user", 100, 110, bin_width=0.5)
        ch = model.channel("sr", x)
        ch.add("signal", mc.Fit("gauss", self.sig, nuisances={"pres": {"sigma": "10%"}}), signal=True)
        # a user RooFit PDF with a nuisance inside: 'escale' is held during the fit and floats in Combine
        escale = ROOT.RooRealVar("escale", "", 0., -5., 5.)
        slope = ROOT.RooRealVar("bkg_slope", "", -0.1, -1., 0.)
        rate = ROOT.RooFormulaVar("bkg_slope_eff", "@0*(1 + 0.05*@1)", ROOT.RooArgList(slope, escale))
        bkg = ROOT.RooExponential("bkg_pdf", "", x, rate)
        ch.add("bkg", mc.Fit(bkg, self.bkg))
        model.param("escale", 0., 1.)
        out = os.path.join(self.tmp, "userpdf")
        model.write(out)
        with open(os.path.join(out, "combined.txt")) as f:
            card = f.read()
        self.assertIn("pres param 0 1", card)
        self.assertIn("escale param 0 1", card)
        self.assertEqual(card.count("escale param"), 1)
        self.assertEqual(escale.getVal(), 0.)
        self.assertFalse(escale.isConstant())
        self.assertTrue(slope.isConstant())
        self.assertLess(slope.getVal(), -0.05)   # fitted to the falling background

    def test_conflicting_nuisance_types(self):
        model = mc.Model("conflict")
        ch = model.channel("sr", mc.observable("p_c", 100, 110, bin_width=0.5))
        ch.add("signal", self.sig, signal=True)
        model.lnN("pscale", 1.1)
        model.shape_shift("pscale", shift=0.1)
        with self.assertRaises(mc.ModelError):
            model.write(os.path.join(self.tmp, "conflict"))

    def test_bad_binning_and_observables(self):
        model = mc.Model("badbins")
        ch = model.channel("sr", mc.observable("p_b", 100, 110, bin_width=0.25))
        ch.add("signal", self.sig.Rebin(2, f"sig_rebin{next(_ids)}"), signal=True)
        with self.assertRaises(mc.ModelError):
            model.write(os.path.join(self.tmp, "badbins"))
        model = mc.Model("badobs")
        model.channel("a", mc.observable("p_o", 100, 110, bins=10)).add("signal", self.sig, signal=True)
        model.channel("b", mc.observable("p_o", 100, 120, bins=10)).add("signal", self.sig, signal=True)
        with self.assertRaises(mc.ModelError):
            model.write(os.path.join(self.tmp, "badobs"))

    def test_commands_select_pois(self):
        model, out, _ = self.build("pois")
        model.add_poi("r2")
        text = model.commands({})
        self.assertNotIn("--expectSignal", text)
        self.assertIn("--redefineSignalPOIs r2 --setParameters r2=0,r=0 --freezeParameters r", text)

    @unittest.skipUnless(shutil.which("text2workspace.py"), "Combine is not set up")
    def test_text2workspace(self):
        model, out, _ = self.build("t2w", fit=True)
        model.text2workspace(out)
        self.assertTrue(os.path.exists(os.path.join(out, "model.root")))


if __name__ == "__main__":
    unittest.main()
