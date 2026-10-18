"""A self-contained mu2ecombine example on synthetic histograms (no input files needed).

    python3 toy_model.py <output dir> [--no-mcstats] [--plots]

It shows the main ways to define a process:
  - a signal peak fit with a double-sided Crystal Ball, with a momentum-scale nuisance on its mean
  - a background template whose normalization floats, constrained by a control region (rateParam)
  - a background whose shape comes from the control region and the signal region combined
  - a sparse background smoothed by a fit (binned=True), with up/down histogram variations
  - RooFit PDFs built by hand: one used as is with a constrained parameter, one fit to MC with a nuisance inside
  - a second signal with its own POI
Then run Combine in the output directory (commands.sh).
"""

import argparse

import numpy as np
import ROOT

import mu2ecombine as mc

parser = argparse.ArgumentParser()
parser.add_argument("output")
parser.add_argument("--no-mcstats", action="store_true")
parser.add_argument("--plots", action="store_true")
args = parser.parse_args()

rng = np.random.default_rng(2026)


def mc_hist(name, x, expected):
    """TH1D of MC values x (0.1 MeV/c bins), normalized to `expected` events."""
    h = ROOT.TH1D(name, "", 300, 90., 120.)
    h.Sumw2()
    x = np.ascontiguousarray(x, dtype=float)
    h.FillN(len(x), x, np.full(len(x), expected / len(x)))
    return h


signal = mc_hist("signal", rng.normal(104.97, 0.35, 50000) - rng.exponential(0.4, 50000), 1.)
dio = mc_hist("dio", 105. - 15. * rng.power(6., 200000), 3.)            # falling spectrum toward an endpoint
cosmic_sr = mc_hist("cosmic_sr", rng.uniform(90, 120, 400), 2.)        # few MC events in the signal region
cosmic_cr = mc_hist("cosmic_cr", rng.uniform(90, 120, 4000), 20.)      # same process, no CRV veto (control region)
rpc = mc_hist("rpc", 90 + rng.exponential(6., 300), 0.5)                # sparse: smoothed by a fit
rpc_up, rpc_down = mc.transformed(rpc, stretch=0.1), mc.transformed(rpc, stretch=-0.1)
other = mc_hist("other", rng.normal(104.97, 0.35, 20000), 1.)           # second signal hypothesis
flat_mc = mc_hist("flat_mc", 100 + rng.exponential(30., 5000), 0.5)

model = mc.Model("toy", poi_range=(0, 20))
model.add_poi("r_other", range=(0, 20))
p = mc.observable("p", 100, 110, bin_width=0.2, title="p [MeV/c]")      # a RooRealVar with 50 bins

sr = model.channel("sr", p, blind=(103, 106))
sr.add("signal", mc.Fit("dscb", signal, nuisances={"pscale": {"mean": 0.05}}), signal=True)
sr.add("other", mc.Fit("gauss", other, nuisances={"pscale": {"mean": 0.05}}), signal=True, poi="r_other")
sr.add("dio", mc.Fit("endpoint", dio, params={"e0": (105.5, 105.0, 106.0)}))
sr.add("cosmic", mc.Fit("chebychev", [cosmic_sr, cosmic_cr], rate=cosmic_sr, order=1))
sr.add("rpc", mc.Fit("exponential", rpc, variations={"rpc_res": (rpc_up, rpc_down)}))

# a RooFit PDF used as is: its parameter 'flat_slope' floats in Combine with a Gaussian constraint
slope = ROOT.RooRealVar("flat_slope", "slope", -0.01, -0.1, 0.1)
flat = ROOT.RooExponential("flat_bkg", "flat", p, slope)
sr.add("flat", mc.Pdf(flat, rate=0.5, variations={"flatN": (1.3, 0.8)}))
model.param("flat_slope", -0.01, 0.005)

# a RooFit PDF fit to MC with a nuisance inside it: 'tail_scale' (declared with model.param) is held at 0 during the
# fit, while 'tail_c' is fit and frozen
tail_scale = ROOT.RooRealVar("tail_scale", "", 0., -5., 5.)
tail_c = ROOT.RooRealVar("tail_c", "", -0.05, -1., 0.)
tail_c_eff = ROOT.RooFormulaVar("tail_c_eff", "@0*(1 + 0.1*@1)", ROOT.RooArgList(tail_c, tail_scale))
sr.add("tail", mc.Fit(ROOT.RooExponential("tail_pdf", "", p, tail_c_eff), flat_mc))
model.param("tail_scale", 0., 1.)

# control region: the cosmic template is its own normalization reference
cr = model.channel("cr", mc.observable("p_cr", 90, 120, bin_width=1., title="p [MeV/c]"))
cr.add("signal", signal, signal=True)
cr.add("cosmic", cosmic_cr)
cr.add("dio", dio)
model.rate_param("k_cosmic", processes="cosmic")                  # one free factor for cosmic in sr and cr

model.lnN("lumi", 1.1, processes="!cosmic", group="rate")
model.lnN("dioN", (1.05, 0.97), processes="dio", group="rate")    # asymmetric (up, down)
model.shape_shift("pscale_mc", shift=0.05, processes="dio")      # generated +-0.05 MeV/c shift, refit
if args.no_mcstats:
    model.mc_stats(enabled=False)

model.write(args.output, plots=args.plots)
