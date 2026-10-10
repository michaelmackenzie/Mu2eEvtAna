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
- beam-intensity re-weighting,
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
| [plotter/examples/run1b/](plotter/examples/run1b/README.md) | Run 1B stack plots of every selection's variables (RMC, RPC, CE, protons, neutrons) |

```bash
root.exe -q -b 'Mu2eEvtAna/plotter/examples/run1b/make_plots.C(".", "figures/run1b")'
root.exe -q -b -l 'Mu2eEvtAna/plotter/test/test_plotter.C("<scratch dir>")'   # plotter self-test
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
