# Mu2eEvtAna

Histogramming and selections for Mu2e EventNtuple inputs, built on EventNtuple's `rooutil`, with a ROOT plotting tool
for signal/background stacks.

## Build

From a muse work directory that holds `Mu2eEvtAna` and `EventNtuple`:

```bash
git clone git@github.com:michaelmackenzie/Mu2eEvtAna.git
mu2einit          # source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh
muse setup
muse build -j<N>
```

Run every command below from the work directory, after `mu2einit` and `muse setup`.

## Quick test

```bash
root.exe -q -b "Mu2eEvtAna/test/test.C"        # Mu2eEvtAna on a short file list
root.exe -q -b "Mu2eEvtAna/test/test_run1b.C"  # Run1BAna
```

`test_conv.C`, `test_rmc.C`, and `test_eventntuple.C` cover the other analyzers and the raw `rooutil` reading.

## Making histograms

`scripts/make_histograms.C` runs an analyzer over a dataset, a file list, or a single ntuple:

```bash
# make_histograms(processes, dataset, mode, function, threads, max_entries, name_tag)
root.exe -q -b 'Mu2eEvtAna/scripts/make_histograms.C(1, "cele1b1s5r0101", 1, "cnv_ana", 4, 1e6)'
root.exe -q -b 'Mu2eEvtAna/scripts/make_histograms.C(1, "my_files.txt", 1, "run1b_ana", 1, 1e5, "mytag")'
```

- **Datasets** are defined in `scripts/datasets.C` as a key, the ntuple dataset name, N(events), and N(generated).
  Their file lists are `file_lists/<dataset name>.files`. With an empty dataset, every entry flagged for processing runs.
- **Output** is `<Analyzer>.<function>.<tag>.m<mode>.root` in the current directory, e.g.
  `Run1BAna.run1b_ana.fele0b1s51r0204.m1.root`. Threads write per-thread files that are merged
  (`scripts/merge_hist_files.sh`).
- **Beam intensity:** MC is re-weighted from its simulated N(POT) per microbunch, read from the dataset tag (Run 1B
  `b1r02XX`: 5.92e6, `b1r00XX`/`b1r01XX`: 1.58e7, `b2`: 3.93e7), to the analyzer's `beam_mu_goal_` (Run1BAna: 6.14e6,
  others: off). Mode 2 re-weights every sample but Run 1B to 2BB (3.93e7). Cosmic-ray and no-primary samples follow a
  log-normal intensity distribution, all others x * log-normal. A name without a dataset tag is not re-weighted.
- **Remote files:** file lists are read over xrootd, which needs a token (`getToken`). A single local `/pnfs` file is
  read directly.
- Keep threads plus any other jobs within the node's shared core budget.

### Analyzers

| Function | Class | Purpose |
|---|---|---|
| `mu2e_ana` | `Mu2eEvtAna` | Base event, track, calorimeter, and CRV histograms |
| `cnv_ana` | `ConvAna` | Conversion-electron/positron selections (Run 1A) |
| `rmc_ana` | `RMCAna` | RMC selections |
| `run1b_ana` | `Run1BAna` | Run 1B (field-off) calorimeter selections: RMC, RPC, CE, protons |

Other packages register their own analyzers with `RegisterAnalyzer(function, class, libraries)` (see
`scripts/functions.C`).

`Run1BAna` applies the per-event physics weights and MC vetoes described in its header (`inc/Run1BAna.hh`):
- DIO and Plestid RMC spectrum weights for the flat-spectrum samples,
- the pion survival weight for RPC,
- beam-intensity re-weighting to 1.5 kW (6.14e6 POT per microbunch),
- a veto of non-primary clusters in primary samples, and of the RMC/proton/neutron/DIO tails in pileup.

It stops with an error when an ntuple lacks the information a weight or veto needs.

### Output layout

```
Ana/Hist/<type>_<set>/<hist>   # type: evt, trk, cls, tcs, lns, crv, ...; set: selection index
Ana/data/Norm                  # TTree with ngen, nntuple, nseen, naccept, nneg
```

## Plotting

`plotter/` is a ROOT-macro tool for stacked backgrounds with any number of signals, stacked or overlaid, and ratio,
difference, or significance lower panels. It holds no physics: each analysis gives its samples' normalizations, and
shared constants are in `plotter/physics/Mu2ePhysics.C`. See [plotter/README.md](plotter/README.md).

| Example | What it makes |
|---|---|
| [plotter/examples/run1a/](plotter/examples/run1a/README.md) | Run 1A μ⁻ → e⁻ and μ⁻ → e⁺ stack plots of the track variables |
| [plotter/examples/run1b/](plotter/examples/run1b/README.md) | Run 1B stack plots of every selection's variables (RMC, RPC, CE, protons, neutrons) |

```bash
root.exe -q -b 'Mu2eEvtAna/plotter/examples/run1b/make_plots.C(".", "figures/run1b")'
root.exe -q -b -l 'Mu2eEvtAna/plotter/test/test_plotter.C("<scratch dir>")'   # plotter self-test
```

## Statistical models (Combine)

`combine/` holds `mu2ecombine`, a Python package that turns the histogram files (or any ROOT TH1s and RooFit PDFs) into
[Combine](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/) datacards and shape files. A short script
declares the channels, processes, and uncertainties. The package writes ordinary cards and a `commands.sh`, and you run
Combine yourself for the limits, fits, and post-fit figures. See [combine/README.md](combine/README.md).

```python
import mu2ecombine as mc
p = mc.observable("p", 100, 110, bin_width=0.2, title="p [MeV/c]")
model = mc.Model("my_search", poi_range=(0, 50))
sr = model.channel("sr", p)
sr.add("signal", mc.Fit("dscb", signal_th1), signal=True)
sr.add("dio", dio_th1)
model.lnN("lumi", 1.1)
model.write("output/my_search")        # then run Combine there, e.g. the commands in commands.sh
```

## Python

`python/` produces the same histogram files from EventNtuples with
[`pyevtana`](https://github.com/michaelmackenzie/pyevtana); see [python/README.md](python/README.md).

## Layout

| Directory | Contents |
|---|---|
| `inc/`, `src/` | Analyzer classes and the event, track, cluster, and histogram types |
| `scripts/` | Processing (`make_histograms.C`, `datasets.C`, `functions.C`) and helper scripts |
| `file_lists/` | Ntuple file lists, one per dataset |
| `data/` | MVA weight files |
| `test/` | Test macros |
| `plotter/` | Plotting tool and examples |
| `python/` | Python histogramming |
| `combine/` | Combine model building (`mu2ecombine`) |
