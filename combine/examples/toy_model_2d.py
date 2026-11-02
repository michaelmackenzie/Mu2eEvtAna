"""A self-contained mu2ecombine example of a 2D (p, t) channel on synthetic histograms (no input files needed).

    python3 toy_model_2d.py <output dir> [--plots]

The signal region is a channel on two observables, momentum p and time t. It shows the ways to define a 2D process:
  - a signal that is a Product of a double-sided Crystal Ball fit in p and a template in t
  - a background that is a Product of a template in p and an exponential fit in t
  - a background with a p-t correlation, as a TH2 template (a Product cannot describe it)
  - a momentum-scale shape uncertainty (axis "p") and a time-offset shape uncertainty (axis "t")
Then run Combine in the output directory (commands.sh).
"""

import argparse

import numpy as np
import ROOT

import mu2ecombine as mc

parser = argparse.ArgumentParser()
parser.add_argument("output")
parser.add_argument("--plots", action="store_true")
args = parser.parse_args()

rng = np.random.default_rng(2026)
PLIM, TLIM = (100., 110.), (450., 1650.)


def hist1(name, x, expected, nbins, lim):
    h = ROOT.TH1D(name, "", nbins, *lim)
    h.Sumw2()
    x = np.ascontiguousarray(x, dtype=float)
    h.FillN(len(x), x, np.full(len(x), expected / len(x)))
    return h


def hist2(name, x, y, expected, nx=100, ny=120):
    h = ROOT.TH2D(name, "", nx, *PLIM, ny, *TLIM)
    h.Sumw2()
    x, y = np.ascontiguousarray(x, dtype=float), np.ascontiguousarray(y, dtype=float)
    h.FillN(len(x), x, y, np.full(len(x), expected / len(x)))
    return h


n = 100000
# signal: a peak in p with a low tail, a prompt-like time distribution (shifted exponential)
sig_p = hist1("sig_p", rng.normal(104.97, 0.35, n) - rng.exponential(0.4, n), 1., 100, PLIM)
sig_t = hist1("sig_t", 450. + rng.exponential(450., n), 1., 120, TLIM)
# a background with factorized shapes: falling spectrum in p, exponential in t
bkg_p = hist1("bkg_p", 100. + 10. * rng.beta(1., 3., n), 8., 100, PLIM)
bkg_t = hist1("bkg_t", 450. + rng.exponential(700., n), 8., 120, TLIM)
# a background whose time distribution depends on p (later at higher p): a TH2 template keeps the correlation
cp = 100. + 10. * rng.beta(1.5, 3., n)
corr = hist2("corr", cp, 450. + rng.exponential(300. + 40. * (cp - 100.), n), 5.)

model = mc.Model("toy2d", poi_range=(0, 20))
p = mc.observable("p", *PLIM, bin_width=0.5, title="p [MeV/c]")
t = mc.observable("t", *TLIM, bin_width=100., title="t [ns]")

sr = model.channel("sr", (p, t))
sr.add("signal", mc.Product(mc.Fit("dscb", sig_p, nuisances={"pscale": {"mean": 0.05}}), sig_t), signal=True,
       title="Signal")
sr.add("dio", mc.Product(bkg_p, mc.Fit("exponential", bkg_t)), title="DIO")
sr.add("corr", corr, title="Correlated")

model.lnN("lumi", 1.1)
model.shape_shift("pscale_mc", shift=0.03, processes="corr", axis="p")        # a generated +-0.03 MeV/c shift of p
model.shape_shift("toffset", shift=20., processes=["dio", "corr"], axis="t")  # +-20 ns of time offset

model.write(args.output, plots=args.plots)
