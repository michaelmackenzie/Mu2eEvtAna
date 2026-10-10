# plotter: signal/background histogram plotting

A ROOT-macro plotting tool for Mu2e histogram files: stacked backgrounds, any number of
signals, data, and a lower pad. It merges three existing plotters:

- `mumep_ana/analysis/plotter/Plotter.C`: the stack/ratio/systematic/ROC plots
- `BNVAna/analysis/plotter/Plotter.C`: several signals, each with its own rate; signal-only plots
- `Run1BAna/scripts/plotSigvsBkgFromNtuple.C`: the signal stacked on the backgrounds and also drawn
  alone in front of the stack, an S/σ(B) lower pad, and shape-normalized comparisons

The plotter holds no physics. Each analysis says, per histogram file, how many events one
MC event stands for, and the plotter multiplies. Constants that several analyses share
(muon capture on Al, DIO fractions, RPC/RMC rates, beam parameters) are kept in one place,
`physics/Mu2ePhysics.C`.

## Files

| file | contents |
|---|---|
| `Plotter.C` | `mu2eplot::Plotter`: loading, normalization, and all plots |
| `PlotTypes.C` | `Process_t`, `Exposure_t`, `Layout_t`, `plot_t`, `mc_norm`, and the enums |
| `PlotUtils.C` | range helpers, `format_sci`, the Mu2e stamp (`Stamp_t`), `handle_canvas` |
| `physics/Mu2ePhysics.C` | shared constants in `namespace mu2e_physics` (plus Run 1A/Run 2 campaign numbers) |
| `physics/Exposures.C` | `Exposure_t` presets: `run1a_exposure()`, `run2_exposure()`, `exposure_1bb(livetime, ...)` |
| `examples/bnv_run1a.C` | BNV signals over the Run 1A backgrounds (replaces `BNVAna/analysis/make_plots.C`) |
| `examples/run1b_rmc.C` | Run 1B calorimeter RMC with the signal stacked (replaces `plotRMCvsBkgFromNtuple.C`) |
| `examples/run1b/` | Run 1B stack plots of the selection variables from the Mu2eEvtAna `Run1BAna` outputs (see its README) |
| `test/test_plotter.C` | toy-file test of the normalization arithmetic and every drawing mode |

Everything is in a namespace (`mu2eplot`, `mu2e_physics`). That way the tool can be loaded next
to an analysis' own `physics.C`/`types.C` without name clashes.

## Quick start

```cpp
#include "Mu2eEvtAna/plotter/Plotter.C"
#include "Mu2eEvtAna/plotter/physics/Exposures.C"
using namespace mu2eplot;
namespace phys = mu2e_physics;

Plotter p;
p.figdir_   = "figures/my_plots";
p.layout_   = Layout_t::evtana();          // Ana/Hist/<type>_<set>/<hist>, Ana/data/Norm
p.exposure_ = phys::run1a_exposure();

// Signal: rate per POT for R = 1, times R through the scale (shown in the legend as [R = ...])
const double rate_capture = phys::muon_capture_fraction_al*phys::run1a::nmuons_per_pot;
p.add_signal("ce", "#mu^{-}#rightarrowe^{-}", kBlue, "ce.root", mc_norm(rate_capture, 1e7), kPOT, 1e-13, "R_{#mue}");

// Backgrounds, in stack order (first at the bottom)
p.add_background("dio", "DIO", kMagenta-10, "dio.root",
                 mc_norm(phys::muon_decay_fraction_al*phys::run1a::nmuons_per_pot*phys::dio_frac_above_95, 25e6))
 .expected(9368976)      // N(events) in the input; corrects for events not processed
 .sys("beam", 0.10);     // 10% normalization uncertainty, source "beam"
p.add_background("cosmic", "Cosmic ray", kAzure-4, "cosmic.root", mc_norm(1., 5.72e6), kLivetime)
 .sys("cosmic", 0.20);
p.add_data("data", "data.root");

if(p.init()) return;
p.print_stack(plot_t("p", "trk", 1, 2, 95., 110., 1., -1., true, false, "p", "MeV/c"));
p.print_yields("p", "trk", 1, 103.85, 105.1);
```

Run it with `root.exe -q -b -l macro.C` from the muse work directory (or use relative
includes, as the examples do). ROOT alone is enough; no libraries need to be loaded.

## Normalization

Each process' histogram is scaled by

```
norm_ × exposure(exposure_) × scale_ × N(expected)/N(seen)
```

- `norm_`: expected events per MC event per unit of exposure. `mc_norm(rate, ngen, emin, emax)`
  gives `rate/ngen`, times `(emax - emin)` for samples generated flat in a range.
- `exposure_`: which number from `Plotter::exposure_` to multiply by. The options are `kPOT`, `kLivetime`
  (cosmic rays), `kNEvents` (pileup-only samples), `kMuonStops`, and `kFixed` (data).
- `scale_`: an extra factor, such as a signal branching ratio. If `scale_symbol_` is set, the legend
  shows `[<symbol> = <scale>]`. `set_signal_scale(R)` changes it for all signals, or for the named ones.
- `N(expected)/N(seen)`: set with `.expected(n)`. `N(seen)` is read from `Layout_t::norm_path_`, either
  by summing a TTree branch or from a TH1's entries or integral.

Data is never rescaled. If a data input is missing events (`N(seen) < N(expected)`), the
exposure is reduced by that fraction instead (`scale_exposure_to_data_`). `exposure_used_`
holds the result. `exposure_` stays as configured, so calling `init()` again doesn't apply the correction twice.

Processes with the same legend label are summed, for example the RMC external 0n and 1n samples, or a
signal split across pileup categories.

Control-region shapes: `.offset(1000, scale, use_nominal_norm)` reads set `selection + 1000`.
If `use_nominal_norm` is true, the plot keeps the nominal set's yield and takes only the shape
from the offset set. This is the same set-offset rule as the mumep plotter. Turn it off with `use_offsets_ = false`.

## Signals

- `signal_mode_ = kOverlay` (default): signals are drawn as lines over the background stack.
- `signal_mode_ = kStacked`: signals go on top of the stack (`signal_stack_fill_`). With
  `signal_front_` (default true), each signal is also drawn alone in front of the stack as a dashed line.
  The line color is `signal_front_color_`, or the signal's color if that is not set. Its legend entry shows both.
  `stacked_signals_ = {"name", ...}` stacks only the named signals and overlays the rest.

## Lower pad (`lower_pad_`)

| mode | with data | without data |
|---|---|---|
| `kRatio` | data / (stacked total), with stat and stat+normalization bands | each signal / background |
| `kDifference` | data − (stacked total) | each signal − background |
| `kSignificance` | per signal: S/√B, or S/√(B + σ²) with `significance_sys_` | same |
| `kNoPad` | single pad | single pad |

`significance_range_` sets which events go into S and B for each bin: `kPerBin` (that bin only),
`kAbove` (everything from the bin up, as for a lower cut), or `kBelow` (everything up to the bin).

The significance is always computed from the absolute yields, even when the plot is
shape-normalized (`plot_t(...).normalized()`).

## Normalization uncertainties

`.sys(source, fraction)` gives a process a fractional normalization uncertainty. Processes
with the same source add linearly (fully correlated). Different sources add in quadrature.
For example, the old mumep ad-hoc band is
`sys("beam", 0.10)` on every beam background, `sys("cosmic", 0.20)` on cosmic rays,
`sys("dio", 0.025)` on DIO, and `sys("rpc_a", 0.093).sys("rpc_b", 0.27)` on RPC.
Uncorrelated per-process uncertainties use one source per process, as in `examples/run1b_rmc.C`.

## File layouts (`layout_`)

| preset | histograms | normalization |
|---|---|---|
| `Layout_t::evtana()` | `Ana/Hist/{type}_{set}/{hist}` | `Ana/data/Norm`, sum of `nseen` |
| `Layout_t::convana()` | `Ana/ConvAna_ConvAna/Hist/{type}_{set}/{hist}` | `Ana/ConvAna_ConvAna/data/Norm` |
| `Layout_t::run1bana()` | `Run1BAna/{type}_{set}/{hist}` | `Run1BAna/data/norm` entries |
| `Layout_t::run1bana_calo()` | `hist_{set}/{hist}` | `norm` integral |

Any other layout is a `hist_format_` string with `{type}`, `{set}`, and `{hist}`.

## Plots

- `plot_stack` / `print_stack(plot_t)`: the main plot. `print_*` writes `figdir_/<prefix>_<hist>_<type>_<set>[_norm][_log].<fmt>` for each of `formats_`.
- `plot_component(plot, name_tag)`: the processes whose names contain the tag, on their own.
- `plot_systematic(plot)`: the model with up/down shifts read from `<hist>_<sys>` in type `sys`.
- `plot_roc(plot, left, eff)`: the ROC curve of each signal against the total background, or the efficiency curves.
- `print_yields(hist, type, set, xmin, xmax)`: yields for each process, the background total with its MC-statistical and normalization uncertainties, and S/√B and S/σ(B) for each signal.
- `get_histograms` / `get_histogram(hist, type, set, role, name_tag)`: the scaled histograms, for custom plots.

Text, stamp, and style options are public members (`legend_columns_`, `stamp_`, `data_label_`, ...).
The stamp's exposure text comes from `exposure_used_` and can be replaced with `stamp_.text_`.

## Validation

- `test/test_plotter.C` builds toy files with known content and checks the normalization, sampling
  correction, offsets, label merging, uncertainty combination, and significance. It also draws every mode.
  `root.exe -q -b -l 'Mu2eEvtAna/plotter/test/test_plotter.C("<scratch dir>")'` returns the number of failures.
- `examples/bnv_run1a.C` reproduces all 48 per-process yields of `BNVAna/analysis/make_plots.C`
  (sets 200, 201, 204, 207, p > 105 MeV/c) exactly.
- `examples/run1b_rmc.C` reproduces all 30 per-process yields of `plotRMCvsBkgFromNtuple.C`
  (sets 70, 74, 75) exactly, and its stacked-signal plot and S/σ(B) curve.

## Differences from the plotters it replaces

- Incomplete data reduces the exposure. The mumep/BNV `init` also scaled the data up by
  `N(expected)/N(seen)`, which applied the correction twice.
- Normalization uncertainties are set per process (`.sys`). The mumep version used a fixed formula that matched process names.
- The stamp's beam power is computed (3.86 kW for Run 1A). The mumep/BNV version rounded it to a hard-coded 3.84 kW.
- Fixed: the ROC background total (`bkg->Add(bkg)`), the ROC efficiency axis title, and a crash when cleaning up
  canvases that have no sub-pads (`Empty_Canvas` also deleted the canvas' own frame).
