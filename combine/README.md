# mu2ecombine

Build [Combine](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/) models for Mu2e analyses in Python, from
ROOT histograms and RooFit PDFs. You describe channels, processes, and uncertainties in a short script. The package
writes ordinary Combine datacards and shape files, and Combine does all of the statistics.

```python
import mu2ecombine as mc

p = mc.observable("p", 100, 110, bin_width=0.2, title="p [MeV/c]")   # a RooRealVar with its binning
model = mc.Model("bnv", poi_range=(0, 50))
sr = model.channel("sr", p, blind=(102, 106))
sr.add("signal", mc.Fit("dscb", signal_th1), signal=True)
sr.add("dio", dio_th1)                                                # a TH1 is a binned template
sr.add("cosmic", mc.Template(cosmic_sr, shape=[cosmic_sr, cosmic_cr]))
model.lnN("lumi", 1.1, processes="!cosmic")
model.shape_shift("pscale", shift=0.1)                                # momentum scale, every MC shape
model.write("output/bnv", plots=True)
```

Then, in `output/bnv`, run the commands from `commands.sh`, or any other Combine command:

```bash
text2workspace.py combined.txt -o model.root -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO 'map=sr/signal:r[1,0,50]'
combine -M AsymptoticLimits model.root -t -1
```

## Environment

The package needs PyROOT, numpy, and scipy (e.g. `pyenv rootana`). Building envelopes and running Combine also need a
Combine installation on `PATH` and `LD_LIBRARY_PATH`. Add `Mu2eEvtAna/combine` to `PYTHONPATH`.

```bash
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh
pyenv rootana 2.5.0
source <Combine directory>/env_standalone.sh        # or your site's Combine setup
export PYTHONPATH=<work dir>/Mu2eEvtAna/combine:$PYTHONPATH
```

## Inputs

Histograms are ROOT `TH1`s, from any file or built in memory; a list of TH1s is summed. The package copies its inputs
(as detached `TH1D`s), so closing your files afterwards is safe. MC statistical uncertainties are taken from the bin errors.

```python
f = mc.Mu2eEvtAnaFile("ConvAna.cnv_ana.dio00b1s5r0101.m1.root")
h = f.hist("trk_20/obs")                       # Ana/Hist/trk_20/obs, as a TH1D
up, down = f.variation("Scale", set=20)        # Ana/Hist/sys_20/obs_1, obs_2 (names from inc/Systematics.hh)

dio = mc.Sample(f, weight=npot * rate / n_generated, expected_events=9368976)
h = dio.hist("trk_20/obs")                     # weighted, and scaled by expected / seen events (Norm tree)
up, down = dio.variation("Scale", set=20)

h = ROOT.TFile("any.root").Get("dir/name")     # any TH1 works directly
```

Keep the dataset definitions of an analysis in their own module and import it, e.g. a dict of `Sample`s. Shared
constants come from `plotter/physics/Mu2ePhysics.C`: `mc.mu2e_physics.run1a.npot`.

Ordinary TH1 methods (`Add`, `Scale`, `Rebin`) work as usual. A few helpers return new TH1Ds:

| Helper | Result |
|---|---|
| `mc.rebin_to(h, edges)` | Bins merged into `edges`, which must be edges of `h` (an error otherwise, never an approximation). |
| `mc.total([h1, h2])` | Sum. |
| `mc.transformed(h, shift=, stretch=, about=)` | x -> about + (x - about)(1 + stretch) + shift, contents redistributed. |
| `mc.convolved(truth, response, edges)` | A truth spectrum smeared by a response histogram (x_reco - x_true). |
| `mc.replace_tail(h, x, reference)` | Bins above x take the shape of `reference`, keeping h's content there. |

## Channels

```python
p = mc.observable("p", lo, hi, bin_width=w, title="p [MeV/c]")      # or bins=, or edges=
ch = model.channel(name, p, data=None, blind=None)
```

- The observable is a `RooRealVar`; its binning is the channel binning. `mc.observable` is a shortcut for building one
  and calling `setBins`/`setBinning`, so you can pass your own instead. Build your RooFit PDFs on the same variable.
- Channels may share an observable. Observables with the same name must have the same range and binning.
- Input histograms are rebinned to the channel automatically. Their bin edges must include the channel's edges.
- `data`: a TH1 of observed data. By default the data is the Asimov background expectation.
- `blind`: a (lo, hi) window that envelope fits exclude. It is also drawn on the figures.

## Processes: shape and rate

`ch.add(name, shape, signal=False, poi=None, title=None)`, where `shape` is one of:

| Shape | Use |
|---|---|
| `th1` or `[th1s]` | Binned template. The rate is its integral in the channel range. |
| `Template(th1, shape=[cr1, cr2])` | Rate from `th1`, shape from other histograms (control regions, merged selections). |
| `Fit(function, th1, rate=None, order=, params=, fit_range=)` | A PDF fit to the MC histogram, then frozen. `function` is any `RooAbsPdf` built on the channel observable, or a predefined name. The rate defaults to the histogram integral. |
| `Fit(..., binned=True)` | The fitted function used as a binned template, to smooth sparse MC. The rate's MC uncertainty is kept. |
| `Pdf(roo_pdf, rate, norm=None)` | Any `RooAbsPdf`, used as is. Its free parameters float in Combine; constrain them with `model.param`. |
| `Envelope(families=, absorbs=[...])` | Data-driven background: function families fit to the data outside the blind window, then discrete profiling (`RooMultiPdf`). `absorbs` lists the MC processes it replaces. |

Predefined names for `Fit` build ordinary RooFit PDFs: `gauss`, `cb`, `dscb` (`RooCrystalBall`), `exponential`, `power`,
`chebychev`, `bernstein`, `exp_poly`, `endpoint`, and `uniform` (see `mu2ecombine/pdfs.py`). Their starting values come
from the histogram, and `params={name: value | (init, lo, hi)}` overrides them. Gaussian cores are Crystal Balls with far
tails, because Combine treats a top-level `RooGaussian` process PDF as a constraint term.

For your own shape, build it in RooFit, e.g. `ROOT.RooGenericPdf("rmc", "pow(p - x0, a) * exp(-l * p)",
ROOT.RooArgList(p, x0, a, l))`. `Fit` fits its non-constant parameters to the histogram and freezes them; constant
parameters stay as you set them.

A channel whose processes are all binned is written as TH1 shapes. Any other channel becomes a `RooWorkspace`, with its
templates stored as `RooDataHist` on the channel observable. `Model(..., workspace=True)` (or `channel(...,
workspace=True)`) writes template channels as workspaces too. Combine maps TH1 shapes to its own observable
(`CMS_th1x`), so only a workspace model can generate data, e.g. `combine -M GenerateOnly --saveToys`, that a PDF model
of the same channel fits with `--toysFile`.

## Uncertainties

| Call | Effect |
|---|---|
| `model.lnN(name, 1.1)` or `model.lnN(name, (k_up, k_down))` | Log-normal rate uncertainty. |
| `Template(h, variations={name: (up, down)})` | Shape and rate from up/down histograms, e.g. Mu2eEvtAna `sys_<set>` variations. |
| `variations={name: (up, down, "shape_only")}` | Up/down normalized to the nominal rate. |
| `variations={name: (up, down, "rate")}` | Only the change in integral, as an asymmetric log-normal. |
| `Fit(..., variations={name: (up, down)})` | The function is refit to up and down, so its parameters move with the nuisance. The rate change is included unless the mode is `"shape_only"`. |
| `Fit(..., nuisances={name: {"mean": 0.05, "sigma": "10%"}})` | A nuisance shifting fitted parameters (absolute, or relative with `%`). For your own PDF, use its RooRealVar names. |
| A nuisance inside your PDF | Put a `RooRealVar` named after the nuisance into your PDF and declare it with `model.param(name, 0, 1)` (or use the name of another nuisance in the model). `Fit` holds it during the fit, and it floats in Combine. |
| `Pdf(..., variations={name: (k_up, k_down)})` | Rate factors, or up/down histograms used for their integrals. |
| `model.shape_shift(name, shift=, stretch=)` | Generated up/down shapes (x shifted, or stretched about the peak) for every selected Template and Fit process. |
| `model.rate_param(name, processes, channels)` | A free normalization shared across channels, e.g. constrained by a control region. |
| `model.param(name, mean, sigma)` | Gaussian constraint on a parameter of your own PDF. |
| `model.mc_stats(enabled=True)` | MC statistics, on by default: one log-normal per process (the MC uncertainty of its yield). `auto=True` uses `autoMCStats` in binned channels instead; it is off because it stalls fits in some Combine builds. |
| `group=` on `lnN`/`shape_shift`, `model.group(name, [...])` | Nuisance groups, for `--freezeNuisanceGroups`. |

- A nuisance name used in several processes or channels is one correlated parameter. This holds across templates,
  fitted PDFs, and workspaces: the workspace variables carry the nuisance's name.
- Fitted parameters become `p + theta * shift` through RooFit's `RooCustomizer`, so the workspace holds ordinary RooFit
  objects (`<param>_eff` formulas).
- Selections (`processes=`, `channels=`) take a name, a glob, or a list. A `!pattern` entry excludes.
- A name cannot be both `lnN` and a shape uncertainty: Combine forbids it.

## Signals and POIs

`signal=True` scales a process by the default POI (`r`, range `poi_range`). For several POIs, call
`model.add_poi("r_x", range=(lo, hi))` and pass `poi="r_x"` to the signal. `commands.sh` maps every signal to its POI with
`multiSignalModel` and selects one POI per command (`--redefineSignalPOIs`, the others fixed to zero). Normalize signal
samples to a reference rate (e.g. 1e-15); the limit on `r` is then in units of that reference.

## Output

| File | Content |
|---|---|
| `combined.txt`, `<channel>.txt` | Datacards. Paths are relative to the card, so the directory can be moved. |
| `shapes/<channel>.root` | TH1 shapes, or the workspace `w`. |
| `model.json` | Yields, MC uncertainties, nuisances per process, fit results, and envelope choices. |
| `commands.sh` | text2workspace and common Combine commands for this model. |
| `figures/` (`plots=True`) | Expected stack per channel, and each MC fit over its histogram. |

After `combine -M FitDiagnostics ... --saveShapes --saveWithUncertainties`, run
`python3 -m mu2ecombine.plots fitDiagnostics.<name>.root` in the output directory for pre-fit and post-fit figures with
data/prediction ratios (in `figures/postfit`, with `postfit_yields.json`).

To scan a parameter (e.g. a signal mass), loop in Python and write one model per point.

## Examples

- [examples/toy_model.py](examples/toy_model.py): synthetic, self-contained. It covers every process kind, a control
  region with a `rateParam`, hand-built RooFit PDFs (one with a nuisance inside), and two POIs.
- `python3 -m unittest discover -s tests` runs the tests. The text2workspace test is skipped when Combine is not set up.

## Layout

| Module | Contents |
|---|---|
| `hist.py` | TH1 helpers (`rebin_to`, `total`, `transformed`, `convolved`, `replace_tail`) |
| `inputs.py` | `Mu2eEvtAnaFile` and `Sample`: histograms and systematic variations from Mu2eEvtAna files |
| `physics.py` | `mu2e_physics`: constants read from `plotter/physics/Mu2ePhysics.C` |
| `shapes.py` | Process shapes: `Template`, `Fit`, `Pdf` |
| `pdfs.py` | Predefined RooFit PDFs, fitting, and nuisance parameters on fitted PDFs |
| `envelope.py` | Data-driven `Envelope` backgrounds (`RooMultiPdf`) |
| `model.py` | `Model`, channels, uncertainties, and `commands.sh` |
| `datacard.py` | Datacard and shape file writing |
| `plots.py` | Expected stacks, fit QA, and post-fit figures from FitDiagnostics |
| `combine.py` | Loading the Combine library and reading limit outputs |

## Limitations

- Tested with Combine v10.1.0 and ROOT 6.32.
- `autoMCStats` (`model.mc_stats(auto=True)`) is off by default, because it stalls AsymptoticLimits in some Combine
  builds. The default MC statistics term is one log-normal per process on its yield, which ignores bin-to-bin MC
  fluctuations in the shape.
- The AsymptoticLimits results were validated against earlier, independent implementations of the same analyses. The
  toy-based `HybridNew` command in `commands.sh` has not been checked yet.
- Not yet supported: 2D observables and unfolding (truth-bin signals with regularization).
