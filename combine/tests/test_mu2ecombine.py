"""mu2ecombine tests on synthetic inputs: python3 -m unittest discover -s tests (from combine/)."""

import itertools
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
import warnings

import numpy as np
import ROOT

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import mu2ecombine as mc  # noqa: E402
from mu2ecombine.hist import arrays, arrays2, integral  # noqa: E402
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

    def test_counting_card(self):
        model = mc.Model("counting", poi_range=(0, 10), workspace=True)
        ch = model.channel("sr", mc.observable("p_counting", 104, 105, bins=1))
        ch.add("signal", mc.hist.make([104, 105], [1.5], [0.01], "signal"), signal=True)
        ch.add("bkg", mc.hist.make([104, 105], [0.2], [0.0004], "bkg"))
        model.lnN("bscale", (1.2, 0.9), processes="bkg")
        out = os.path.join(self.tmp, "counting")
        model.write(out, shapes=False)
        with open(os.path.join(out, "combined.txt")) as f:
            card = f.read()
        self.assertNotIn("shapes", card)
        self.assertFalse(os.path.exists(os.path.join(out, "shapes")))
        self.assertIn("0.9000/1.2000", card)
        self.assertIn("mcstat_sr_bkg lnN", card)
        ch.add("other", mc.Template(mc.hist.make([104, 105], [0.1], None, "other"),
                                    variations={"oshape": (mc.hist.make([104, 105], [0.12], None, "u"),
                                                           mc.hist.make([104, 105], [0.09], None, "d"))}))
        with self.assertRaises(mc.ModelError):
            model.write(os.path.join(self.tmp, "counting_shape"), shapes=False)

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

    def test_template_channel_as_workspace(self):
        model = mc.Model("tmplws", workspace=True)
        ch = model.channel("sr", mc.observable("p_tw", 100, 110, bin_width=0.5))
        ch.add("signal", self.sig, signal=True)
        ch.add("bkg", mc.Template(self.bkg, variations={"bscale": (scaled(self.bkg, 1.1), scaled(self.bkg, 0.95))}))
        out = os.path.join(self.tmp, "tmplws")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            summary = model.write(out)
        with open(os.path.join(out, "combined.txt")) as f:
            card = f.read()
        self.assertEqual(summary["channels"]["sr"]["mode"], "workspace")
        self.assertIn("shapes bkg sr shapes/sr.root w:sr_bkg w:sr_bkg_$SYSTEMATIC", card)
        f = ROOT.TFile(os.path.join(out, "shapes", "sr.root"))
        w = f.Get("w")
        self.assertTrue(w.var("p_tw"))
        self.assertAlmostEqual(w.data("sr_bkg").sumEntries(), integral(mc.rebin_to(self.bkg, np.linspace(100, 110, 21))),
                               places=6)
        self.assertTrue(w.data("sr_bkg_bscaleUp"))
        f.Close()

    @unittest.skipUnless(shutil.which("combine"), "Combine is not set up")
    def test_template_model_generates_for_pdf_model(self):
        """Toys of a template model written as a workspace live on the observable of a PDF model of the channel."""
        def write(name, fit):
            model = mc.Model(name, workspace=True)
            ch = model.channel("sr", mc.observable("p_gen", 100, 110, bin_width=0.5))
            ch.add("signal", mc.Fit("gauss", self.sig) if fit else self.sig, signal=True)
            ch.add("bkg", mc.Fit("exponential", self.bkg) if fit else self.bkg)
            out = os.path.join(self.tmp, name)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model.write(out, text2workspace=True)
            return out
        gen, fit = write("gen_tmpl", False), write("gen_pdf", True)
        proc = subprocess.run("combine -M GenerateOnly model.root -t -1 --saveToys -n .gen", shell=True, cwd=gen,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        self.assertEqual(proc.returncode, 0, proc.stdout[-2000:])
        toys = ROOT.TFile(os.path.join(gen, "higgsCombine.gen.GenerateOnly.mH120.123456.root"))
        names = {v.GetName() for v in toys.Get("toys/toy_asimov").get()}
        f = ROOT.TFile(os.path.join(fit, "model.root"))
        observables = {v.GetName() for v in f.Get("w").set("observables")}
        self.assertIn("p_gen", names)
        self.assertTrue(observables <= names, (observables, names))
        toys.Close()
        f.Close()

    def test_unconstrained_nuisances(self):
        model = mc.Model("free", workspace=True)
        ch = model.channel("sr", mc.observable("p_free", 100, 110, bin_width=0.5))
        ch.add("signal", mc.Fit("gauss", self.sig, nuisances={"pres": {"sigma": "10%"}}), signal=True)
        ch.add("bkg", mc.Template(self.bkg, variations={"bshape": (scaled(self.bkg, 1.1), scaled(self.bkg, 0.95))}))
        model.lnN("lumi", 1.1)
        model.unconstrained("bshape", "pres")
        out = os.path.join(self.tmp, "free")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model.write(out)
        card = open(os.path.join(out, "combined.txt")).read()
        self.assertRegex(card, r"(?m)^bshape shapeU ")
        self.assertIn("pres flatParam", card)
        self.assertNotIn("pres param", card)
        self.assertRegex(card, r"(?m)^lumi lnN ")
        model.unconstrained("lumi")
        with self.assertRaises(mc.ModelError):
            model.write(os.path.join(self.tmp, "free_lnn"))
        if shutil.which("text2workspace.py"):
            proc = subprocess.run("text2workspace.py combined.txt -o model.root", shell=True, cwd=out,
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            self.assertEqual(proc.returncode, 0, proc.stdout[-2000:])

    def test_fit_process_alone(self):
        model = mc.Model("fitonly", workspace=True)
        ch = model.channel("sr", mc.observable("p_fo", 100, 110, bin_width=0.5))
        ch.add("signal", mc.Fit("gauss", self.sig), signal=True)
        ch.add("bkg", mc.Fit("exponential", self.bkg))
        figdir = os.path.join(self.tmp, "fitonly")
        r = model.fit_process("sr", "signal", figdir=figdir)
        self.assertEqual(list(r.fits), ["signal"])                  # the other process is not built
        self.assertAlmostEqual(r.fits["signal"].values()["sr_signal_mean"], 105., delta=0.05)
        self.assertTrue(os.path.exists(os.path.join(figdir, "sr_signal_fit.png")))
        self.assertFalse(os.path.exists(os.path.join(figdir, "combined.txt")))
        with self.assertRaises(mc.ModelError):
            model.fit_process("sr", "nothing")

    def _narrow_pair(self, name, fit):
        model = mc.Model(name, workspace=True)
        ch = model.channel("sr", mc.observable(f"p_{name}", 100, 110, bin_width=0.5))
        ch.add("signal", mc.Fit("gauss", self.sig) if fit else self.sig, signal=True)
        ch.add("bkg", mc.Fit("exponential", self.bkg) if fit else self.bkg)
        model.lnN("lumi", 1.1)
        full = os.path.join(self.tmp, name)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model.write(full)
        narrow = os.path.join(self.tmp, f"{name}_narrow")
        used = mc.narrow_model(full, narrow, {"sr": (102.1, 107.)})
        return full, narrow, used

    def test_narrow_templates(self):
        full, narrow, used = self._narrow_pair("narrow_t", fit=False)
        self.assertEqual(used, {"sr": (102., 107.)})                  # moved to the nearest bin edges
        info = json.load(open(os.path.join(narrow, "model.json")))["channels"]["sr"]
        self.assertEqual(info["observable"]["bins"], 10)
        expected = integral(mc.rebin_to(self.bkg, np.linspace(100, 110, 21)).Clone())
        h = mc.rebin_to(self.bkg, np.linspace(100, 110, 21))
        inside = sum(h.GetBinContent(i) for i in range(5, 15))
        self.assertAlmostEqual(info["processes"]["bkg"]["rate"], inside, places=6)
        self.assertLess(inside, expected)
        card = open(os.path.join(narrow, "combined.txt")).read()
        self.assertRegex(card, rf"(?m)^rate\s+\S+\s+{float(f'{inside:.6g}')}")
        self.assertIn("lumi lnN", card)
        f = ROOT.TFile(os.path.join(narrow, "shapes", "sr.root"))
        w = f.Get("w")
        self.assertAlmostEqual(w.var("p_narrow_t").getMin(), 102.)
        self.assertAlmostEqual(w.data("sr_bkg").sumEntries(), inside, places=6)
        self.assertAlmostEqual(w.data("sr_data_obs").sumEntries(), inside, places=6)
        f.Close()

    def test_narrow_pdfs(self):
        full, narrow, _ = self._narrow_pair("narrow_p", fit=True)
        before = json.load(open(os.path.join(full, "model.json")))["channels"]["sr"]["processes"]["signal"]["rate"]
        after = json.load(open(os.path.join(narrow, "model.json")))["channels"]["sr"]["processes"]["signal"]["rate"]
        self.assertGreater(after, 0.8 * before)                       # the peak at 105 is inside 102-107
        self.assertLess(after, before)
        f = ROOT.TFile(os.path.join(narrow, "shapes", "sr.root"))
        w = f.Get("w")
        card = open(os.path.join(narrow, "combined.txt")).read()
        name = re.search(r"(?m)^shapes signal sr \S+ w:(\S+)", card).group(1)
        pdf = w.pdf(name)
        self.assertTrue(pdf)
        self.assertTrue(pdf.dependsOn(w.var("p_narrow_p")))           # attached to the narrowed observable
        self.assertAlmostEqual(w.var("p_narrow_p").getMax(), 107.)
        f.Close()
        if shutil.which("text2workspace.py"):
            proc = subprocess.run("text2workspace.py combined.txt -o model.root", shell=True, cwd=narrow,
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            self.assertEqual(proc.returncode, 0, proc.stdout[-2000:])

    @unittest.skipUnless(shutil.which("text2workspace.py"), "Combine is not set up")
    def test_text2workspace(self):
        model, out, _ = self.build("t2w", fit=True)
        model.text2workspace(out)
        self.assertTrue(os.path.exists(os.path.join(out, "model.root")))


def toy2(x, y, expected, nx=20, ny=8, xlim=(100., 110.), ylim=(0., 4.)):
    """A TH2D of the points (x, y), normalized to `expected` events."""
    h = ROOT.TH2D(f"toy2_{next(_ids)}", "", nx, *xlim, ny, *ylim)
    h.SetDirectory(0)
    h.Sumw2()
    x, y = np.ascontiguousarray(x, dtype=float), np.ascontiguousarray(y, dtype=float)
    h.FillN(len(x), x, y, np.full(len(x), expected / len(x)))
    return h


def edges_of(obs):
    from mu2ecombine.model import binning
    return binning(obs)


def unrolled_values(path, name, nx, ny):
    """The template `name` of a binned channel shapes file as an array [nx, ny]."""
    f = ROOT.TFile(path)
    h = f.Get(name)
    values = np.array([h.GetBinContent(k + 1) for k in range(nx * ny)]).reshape(ny, nx).T
    f.Close()
    return values


class TwoDimensionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(7)
        n = 40000
        cls.sig = toy2(rng.normal(105, 0.5, n), rng.exponential(1., n), 3.)
        cls.bkg = toy2(100 + rng.exponential(4., n), rng.uniform(0, 4, n), 60.)
        cls.gauss2 = toy2(rng.normal(105, 1., n), rng.normal(2., 0.3, n), 10.)
        cls.hx = toy(rng.normal(105, 1, n), 5., 100, 100., 110.)
        cls.hy = toy(rng.exponential(1.2, n), 5., 40, 0., 4.)
        cls.hyg = toy(rng.normal(2., 0.3, n), 5., 40, 0., 4.)
        cls.tmp = tempfile.mkdtemp(prefix="m2c_test2d_")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def axes(self, tag):
        return mc.observable(f"p2_{tag}", 100, 110, bins=10), mc.observable(f"t2_{tag}", 0, 4, bins=4)

    def write(self, model, name):
        out = os.path.join(self.tmp, name)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            summary = model.write(out)
        return out, summary

    def test_rebin2_unroll_roll(self):
        ex, ey = np.linspace(100, 110, 11), np.linspace(0, 4, 5)
        r = mc.rebin2_to(self.bkg, ex, ey)
        self.assertEqual((r.GetNbinsX(), r.GetNbinsY()), (10, 4))
        self.assertFalse(r.GetDirectory())
        self.assertAlmostEqual(r.Integral(), self.bkg.Integral(), places=9)
        _, _, v, var = arrays2(r)
        _, _, fv, fvar = arrays2(self.bkg)
        np.testing.assert_allclose(v[1, 2], fv[2:4, 4:6].sum())
        np.testing.assert_allclose(var[1, 2], fvar[2:4, 4:6].sum())
        u = mc.unroll(r)
        self.assertEqual(u.GetNbinsX(), 40)
        self.assertEqual(u.GetXaxis().GetXmin(), 0.)
        self.assertEqual(u.GetBinContent(2 * 10 + 3 + 1), r.GetBinContent(4, 3))   # k = iy * nx + ix
        self.assertEqual(u.GetBinError(2 * 10 + 3 + 1), r.GetBinError(4, 3))
        back = mc.roll(u, ex, ey)
        np.testing.assert_array_equal(arrays2(back)[2], arrays2(r)[2])
        np.testing.assert_array_equal(arrays2(back)[3], arrays2(r)[3])
        with self.assertRaises(ValueError):
            mc.rebin2_to(self.bkg, ex, [0, 0.3, 4])
        with self.assertRaises(ValueError):
            mc.rebin2_to(self.bkg, [100.05, 110], ey)
        with self.assertRaises(ValueError):
            mc.roll(u, ex, ey[:-1])
        with self.assertRaises(ValueError):
            mc.unroll(self.hx)

    def test_transformed2_along_axes(self):
        v0 = arrays2(self.gauss2)[2]
        for axis in (0, 1):
            t = mc.transformed2(self.gauss2, axis, shift=0.3)
            v = arrays2(t)[2]
            np.testing.assert_allclose(v.sum(axis=axis), v0.sum(axis=axis), rtol=1e-2)    # other projection kept
            self.assertGreater(abs(t.GetMean(axis + 1) - self.gauss2.GetMean(axis + 1)), 0.1)
        with self.assertRaises(ValueError):
            mc.transformed2(self.gauss2, 2)

    def test_2d_template_channel(self):
        px, tx = self.axes("tmpl")
        ex, ey = edges_of(px), edges_of(tx)
        model = mc.Model("tmpl2d", poi_range=(0, 10))
        data = mc.rebin2_to(self.bkg, ex, ey)
        ch = model.channel("sr", (px, tx), data=data)
        self.assertEqual((ch.dim, len(ch.edges)), (2, 41))
        self.assertEqual(ch.obs.GetName(), "sr_bin")
        ch.add("signal", self.sig, signal=True)
        ch.add("bkg", mc.Template(self.bkg, variations={"bscale": (self.bkg, self.bkg)}))
        model.lnN("lumi", 1.1, processes="!bkg")
        out, summary = self.write(model, "tmpl2d")
        info = summary["channels"]["sr"]
        self.assertEqual(info["mode"], "binned")
        self.assertEqual(info["observable"]["unroll"], "k = iy * nx + ix")
        self.assertEqual([a["name"] for a in info["observable"]["axes"]], [px.GetName(), tx.GetName()])
        self.assertEqual(info["observable"]["bins"], 40)
        self.assertAlmostEqual(info["processes"]["bkg"]["rate"], self.bkg.Integral(), places=6)
        self.assertEqual(info["data"]["source"], "data")
        with open(os.path.join(out, "combined.txt")) as f:
            card = f.read()
        self.assertRegex(card, rf"(?m)^observation\s+{float(f'{self.bkg.Integral():.6g}')}")
        sig = unrolled_values(os.path.join(out, "shapes", "sr.root"), "signal", 10, 4)
        np.testing.assert_allclose(sig, arrays2(mc.rebin2_to(self.sig, ex, ey))[2])

    def test_product_of_template_and_fit(self):
        px, tx = self.axes("prod")
        model = mc.Model("prod2d", poi_range=(0, 10))
        ch = model.channel("sr", (px, tx))
        shape = mc.Product(self.hx, mc.Fit("exponential", self.hy), rate=4.)
        ch.add("signal", shape, signal=True)
        ch.add("bkg", self.bkg)
        out, summary = self.write(model, "prod2d")
        info = summary["channels"]["sr"]["processes"]
        self.assertEqual(info["signal"]["kind"], "product")
        self.assertAlmostEqual(info["signal"]["rate"], 4., places=9)
        self.assertEqual(list(info["signal"]["components"]), [tx.GetName()])
        v = unrolled_values(os.path.join(out, "shapes", "sr.root"), "signal", 10, 4)
        self.assertAlmostEqual(v.sum(), 4., places=9)
        fx, fy = v.sum(axis=1), v.sum(axis=0)
        np.testing.assert_allclose(v, np.outer(fx, fy) / v.sum(), rtol=1e-9)     # factorized: the outer product
        hx = mc.rebin_to(self.hx, edges_of(px))
        np.testing.assert_allclose(fx / fx.sum(), arrays(hx)[1] / integral(hx), rtol=1e-9)
        r = model.build()["sr"]
        self.assertIs(r.components["signal"][tx.GetName()], shape.y)
        self.assertEqual(r.mc_rel["signal"], 0.)         # a numeric rate has no MC uncertainty
        ch.add("free", mc.Product(self.hx, self.hy))
        self.assertGreater(model.build()["sr"].mc_rel["free"], 0.)     # default: the x component's histogram

    def test_shape_shift_along_axis(self):
        px, tx = self.axes("shift")
        model = mc.Model("shift2d", poi_range=(0, 10))
        ch = model.channel("sr", (px, tx))
        ch.add("prod", mc.Product(self.hx, self.hyg), signal=True)
        ch.add("th2", self.gauss2)
        model.shape_shift("pscale", shift=0.4, axis=px.GetName())
        model.shape_shift("tscale", shift=0.3, axis=tx.GetName())
        model.shape_shift("other", shift=0.3, axis="nothing")     # matches no axis: skipped
        out, _ = self.write(model, "shift2d")
        path = os.path.join(out, "shapes", "sr.root")
        for proc in ("prod", "th2"):
            nom = unrolled_values(path, proc, 10, 4)
            for nuis, axis in (("pscale", 0), ("tscale", 1)):
                for sign in ("Up", "Down"):
                    v = unrolled_values(path, f"{proc}_{nuis}{sign}", 10, 4)
                    np.testing.assert_allclose(v.sum(axis=axis), nom.sum(axis=axis), rtol=1e-2, atol=2e-3 * nom.sum())
                    self.assertGreater(np.abs(v.sum(axis=1 - axis) - nom.sum(axis=1 - axis)).max(), 1e-3 * nom.sum())
        f = ROOT.TFile(path)
        self.assertFalse(f.Get("prod_otherUp"))
        f.Close()
        up, nom = unrolled_values(path, "prod_pscaleUp", 10, 4), unrolled_values(path, "prod", 10, 4)
        np.testing.assert_allclose(up, np.outer(up.sum(axis=1), nom.sum(axis=0)) / nom.sum(), rtol=1e-6)
        # a 1D channel on the same observable name is shifted by the same nuisance, one on another observable is not
        m1 = mc.Model("shift1d")
        m1.channel("a", mc.observable(px.GetName(), 100, 110, bins=10)).add("bkg", self.hx)
        m1.channel("b", mc.observable("q_other", 100, 110, bins=10)).add("bkg", self.hx)
        m1.shape_shift("pscale", shift=0.4, axis=px.GetName())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = m1.build()
        self.assertIn("bkg_pscaleUp", r["a"].templates)
        self.assertNotIn("bkg_pscaleUp", r["b"].templates)

    def test_shared_nuisance_on_both_components(self):
        px, tx = self.axes("shared")
        hx, hy = mc.rebin_to(self.hx, edges_of(px)), mc.rebin_to(self.hy, edges_of(tx))
        vx = {"nu": (scaled(self.hx, 1.2), scaled(self.hx, 0.9)), "only_x": (scaled(self.hx, 1.1), scaled(self.hx, 0.9))}
        vy = {"nu": (scaled(self.hy, 1.5), scaled(self.hy, 0.7))}
        model = mc.Model("shared2d", poi_range=(0, 10))
        ch = model.channel("sr", (px, tx))
        ch.add("signal", mc.Product(mc.Template(self.hx, variations=vx), mc.Template(self.hy, variations=vy)),
               signal=True)
        self.assertTrue({"nu", "only_x"} <= model.nuisance_names())
        out, summary = self.write(model, "shared2d")
        self.assertEqual(summary["channels"]["sr"]["processes"]["signal"]["nuisances"],
                         ["mcstat_sr_signal", "nu", "only_x"])
        path = os.path.join(out, "shapes", "sr.root")
        nom, up = unrolled_values(path, "signal", 10, 4), unrolled_values(path, "signal_nuUp", 10, 4)
        nx, ny = integral(hx), integral(hy)
        np.testing.assert_allclose(nom, nx * np.outer(arrays(hx)[1] / nx, arrays(hy)[1] / ny), rtol=1e-9)
        np.testing.assert_allclose(up.sum(), nx * 1.2 * 1.5, rtol=1e-9)    # both components vary with the nuisance
        ox = unrolled_values(path, "signal_only_xDown", 10, 4)
        np.testing.assert_allclose(ox.sum(axis=0), 0.9 * nom.sum(axis=0), rtol=1e-9)     # the y shape stays nominal

    def test_errors(self):
        px, tx = self.axes("err")
        model = mc.Model("err2d")
        ch = model.channel("sr", (px, tx))
        with self.assertRaisesRegex(mc.ModelError, "TH1"):
            ch.add("a", self.hx)
        with self.assertRaisesRegex(mc.ModelError, "TH1"):
            ch.add("a", mc.Template(self.hx))
        with self.assertRaisesRegex(mc.ModelError, "Product"):
            ch.add("a", mc.Fit("exponential", self.hy))
        gauss = ROOT.RooGaussian("g_err", "", px, ROOT.RooRealVar("m_err", "", 1.), ROOT.RooRealVar("s_err", "", 1.))
        with self.assertRaisesRegex(mc.ModelError, "Product"):
            ch.add("a", mc.Pdf(gauss, 1.))
        with self.assertRaisesRegex(mc.ModelError, "TH1"):
            model.channel("bad", (px, tx), data=self.hx)
        with self.assertRaisesRegex(mc.ModelError, "same name"):
            model.channel("bad2", (px, px))
        with self.assertRaisesRegex(mc.ModelError, "one RooRealVar"):
            model.channel("bad3", (px, tx, px))
        with self.assertRaises(ValueError):
            mc.Product(self.bkg, self.hy)            # a TH2 component
        with self.assertRaises(ValueError):
            mc.Product(self.hx, mc.Template(self.bkg))
        one = mc.Model("err1d")
        c1 = one.channel("sr", mc.observable("p_err1", 100, 110, bins=10))
        with self.assertRaisesRegex(mc.ModelError, "2D channel"):
            c1.add("a", mc.Product(self.hx, self.hy))
        with self.assertRaisesRegex(mc.ModelError, "1D channel"):
            c1.add("b", self.bkg)
        with self.assertRaisesRegex(mc.ModelError, "TH2"):
            one.channel("d", mc.observable("p_err1b", 100, 110, bins=10), data=self.bkg)
        with self.assertRaisesRegex(mc.ModelError, "edge"):     # a TH2 that does not have the channel's bin edges
            m3 = mc.Model("err2d_c")
            m3.channel("sr", (mc.observable("p_c", 100.05, 110.05, bins=10), mc.observable("t_c", 0, 4, bins=4))
                       ).add("a", self.bkg)
        ch.add("ok", mc.Product(self.hx, self.hy))             # a shift without an axis acts on a 2D channel
        model.shape_shift("s", shift=0.1)
        with self.assertRaisesRegex(mc.ModelError, "axis"):
            model.build()
        model2 = mc.Model("err2d_b")                            # a PDF component on the other axis
        ch2 = model2.channel("sr", self.axes("err2"))
        pdf = ROOT.RooGaussian("g2_err", "", ch2.axes[1], ROOT.RooRealVar("m2_err", "", 1.),
                               ROOT.RooRealVar("s2_err", "", 1.))
        ch2.add("y", mc.Product(mc.Fit(pdf, self.hx), self.hy))
        with self.assertRaisesRegex(mc.ModelError, "other observable"):
            model2.build()

    def test_narrow_refuses_2d(self):
        from mu2ecombine.narrow import NarrowError
        px, tx = self.axes("narrow")
        model = mc.Model("narrow2d")
        model.channel("sr", (px, tx)).add("bkg", self.bkg)
        out, _ = self.write(model, "narrow2d")
        with self.assertRaisesRegex(NarrowError, "2D"):
            mc.narrow_model(out, out + "_narrow", {"sr": (102., 108.)})

    def _plot_model(self, tag):
        px, tx = self.axes(tag)
        model = mc.Model(f"plot2d_{tag}", poi_range=(0, 10))
        ch = model.channel("sr", (px, tx))
        ch.add("signal", mc.Product(mc.Fit("gauss", self.hx), mc.Fit("exponential", self.hy)), signal=True)
        ch.add("bkg", self.bkg)
        return model, px, tx

    def test_plots_2d(self):
        model, px, tx = self._plot_model("plots")
        out = os.path.join(self.tmp, "plots2d")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model.write(out, plots=True)
        names = set(os.listdir(os.path.join(out, "figures")))
        p, t = px.GetName(), tx.GetName()
        expected = {f"sr_prefit_{a}{s}.png" for a in (p, t) for s in ("", "_log")}
        expected |= {f"sr_signal_{a}_fit{s}.png" for a in (p, t) for s in ("", "_log")}
        expected |= {f"sr_{q}_2d.png" for q in ("signal", "bkg")}
        self.assertEqual(expected - names, set())
        self.assertFalse({n for n in names if n.startswith("sr_prefit.")})     # no unrolled stack

    def test_fit_process_of_a_product(self):
        model, px, tx = self._plot_model("fitproc")
        figdir = os.path.join(self.tmp, "fitproc2d")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = model.fit_process("sr", "signal", figdir=figdir)
        self.assertEqual(set(r.components["signal"]), {px.GetName(), tx.GetName()})
        for a in (px.GetName(), tx.GetName()):
            self.assertTrue(os.path.exists(os.path.join(figdir, f"sr_signal_{a}_fit.png")))
        self.assertTrue(os.path.exists(os.path.join(figdir, "sr_signal_2d.png")))
        r = model.fit_process("sr", "bkg", figdir=figdir)        # a TH2 template: its 2D figure only
        self.assertTrue(os.path.exists(os.path.join(figdir, "sr_bkg_2d.png")))

    def test_project_postfit(self):
        from mu2ecombine.plots import project_postfit
        nx, ny = 3, 2
        edges = {"observable": {"title": "u", "edges": list(range(nx * ny + 1)), "axes": [
            {"name": "x", "title": "x [a]", "bins": nx, "edges": [0., 1., 2., 3.]},
            {"name": "y", "title": "y [b]", "bins": ny, "edges": [0., 1., 2.]}]}, "blind": [1, 2]}
        v = np.arange(1., 7.)                                    # k = iy * nx + ix
        e = np.full(6, 2.)
        ch = {"edges": np.arange(7.), "processes": {"a": (v, e)}, "total": (v, e), "data": (v, v, v), "asimov": False}
        px, info = project_postfit(ch, edges, 0)
        np.testing.assert_allclose(px["processes"]["a"][0], [5., 7., 9.])
        np.testing.assert_allclose(px["processes"]["a"][1], np.full(3, np.sqrt(8.)))
        np.testing.assert_allclose(px["data"][1], np.sqrt([5., 7., 9.]))
        self.assertEqual((info["observable"]["title"], info["blind"]), ("x [a]", None))
        py, _ = project_postfit(ch, edges, 1)
        np.testing.assert_allclose(py["processes"]["a"][0], [6., 15.])
        np.testing.assert_allclose(py["edges"], [0., 1., 2.])

    def test_postfit_figures_2d(self):
        """The post-fit figures of a 2D channel from a synthetic FitDiagnostics-like file."""
        from mu2ecombine.plots import plot_fitdiagnostics
        model, px, tx = self._plot_model("postfit")
        out = os.path.join(self.tmp, "postfit2d")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            summary = model.write(out)
        n = summary["channels"]["sr"]["observable"]["bins"]
        f = ROOT.TFile(os.path.join(out, "fitDiagnostics.test.root"), "RECREATE")
        for fit in ("prefit", "fit_s"):
            d = f.mkdir(f"shapes_{fit}").mkdir("sr")
            d.cd()
            for name, k in (("signal", 1.), ("bkg", 5.), ("total", 6.)):
                h = ROOT.TH1D(name, "", n, 0., float(n))
                for b in range(n):
                    h.SetBinContent(b + 1, k * (1. + b % 4))
                    h.SetBinError(b + 1, 0.1 * k)
                h.Write()
            g = ROOT.TGraphAsymmErrors(n)
            for b in range(n):
                g.SetPoint(b, b + 0.5, 6. * (1. + b % 4))
            g.SetName("data")
            g.Write()
        f.Close()
        plot_fitdiagnostics(os.path.join(out, "fitDiagnostics.test.root"), out)
        names = set(os.listdir(os.path.join(out, "figures", "postfit")))
        for fit in ("prefit", "fit_s"):
            for a in (px.GetName(), tx.GetName()):
                self.assertIn(f"sr_{fit}_{a}.png", names)


    @unittest.skipUnless(shutil.which("combine") and shutil.which("text2workspace.py"), "Combine is not set up")
    def test_combine_limit_on_2d_model(self):
        px, tx = self.axes("limit")
        model = mc.Model("limit2d", poi_range=(0, 20))
        ch = model.channel("sr", (px, tx))
        ch.add("signal", mc.Product(self.hx, mc.Fit("exponential", self.hy), rate=3.), signal=True)
        ch.add("bkg", self.bkg)
        model.lnN("lumi", 1.1, processes="signal")
        model.shape_shift("pscale", shift=0.2, axis=px.GetName())
        out, _ = self.write(model, "limit2d")
        model.text2workspace(out)
        proc = subprocess.run("combine -M AsymptoticLimits model.root -t -1 --setParameters r=0 -n .lim2d", shell=True,
                              cwd=out, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        self.assertEqual(proc.returncode, 0, proc.stdout[-2000:])
        match = re.search(r"Expected 50.0%: r < ([0-9.]+)", proc.stdout)
        self.assertIsNotNone(match, proc.stdout[-2000:])
        self.assertGreater(float(match.group(1)), 0.)


if __name__ == "__main__":
    unittest.main()
