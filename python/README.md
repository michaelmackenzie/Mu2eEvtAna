# Mu2eEvtAna in python, via pyevtana

An example of producing Mu2eEvtAna-style histogram files from EventNtuple data using
[`pyevtana`](../../pyevtana), as a python counterpart to
`scripts/make_histograms.C`. Output follows `Mu2eEvtAna::InitializeOutput`:

```
Ana/
  Hist/
    trk_75/  evt_75/      # signal selection
    trk_0/   evt_0/       # p > 90 MeV/c
  data/
    Norm                  # TTree: ngen, nntuple, nseen, naccept, nneg
```

Histogram names and binning are those of `Mu2eEvtAna::BookTrackHist`, so the existing
plotting macros read it unchanged — `make_plots.C` builds
`<dir>Hist/<type>_<set>/<hist>`, which resolves to `Ana/Hist/trk_75/p_2` and so on.

## Setup

```bash
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh
pyenv rootana 2.5.0
export PYTHONPATH=/exp/mu2e/app/users/mmackenz/main/pyevtana:$PYTHONPATH
```

## Use

```bash
python3 make_histograms.py --list                                  # dataset keys
python3 make_histograms.py cele1b1s5r0100 --jobs 10                # a dataset, 10 workers
python3 make_histograms.py cele1b1s5r0100 --jobs 10 --max-files 20 # a quick subset
python3 make_histograms.py cele1b1s5r0100 --jobs 10 --split-files  # per-worker files + hadd
```

Output goes to `Mu2eEvtAna/hist/<key>.hist.root` unless `--output` says otherwise.

## The four pieces

| file | what it holds |
|---|---|
| `datasets.py` | dataset key -> full name, mirroring `scripts/datasets.C`; the file list is `file_lists/<full_name>.files` |
| `selections.py` | selection sets keyed by set ID, and the per-track quantities they need |
| `histograms.py` | the histogram books, one per set ID, plus the `Norm` tree |
| `make_histograms.py` | the driver: resolve a dataset, run N workers, merge, write, print the cut flow |

### Selection sets

Keyed by the set ID that names the output directory, as `Mu2eEvtAna::BookHistograms` does:

- **75** — the signal selection. The same 20 cuts `Mu2eEvtAna` applies; the per-event
  cumulative cut flow was validated stage by stage against pyfitter's implementation of
  the same selection (`pyevtana/benchmarks/validate_against_pyfitter.py`) and agrees
  exactly.
- **0** — a loose `p > 90 MeV/c` set: a good downstream electron-hypothesis track, nothing
  more, for looking at the spectrum before the ID cuts bite.

Adding a set is one entry in `SELECTIONS`; its histograms and directories follow.

### Histograms

A deliberately small subset to start, chosen from what `make_plots.C` plots most:
`p`, `p_2`, `pt`, `t0`, `d0`, `rMax`, `tanDip`, `cosTheta`, `nActive`, `trkQual`, `trkpid`,
`clusterE`, `ep`, plus event-level `nTracks`. Names and binning are copied from
`Mu2eEvtAna::BookTrackHist`. Adding one is a row in `TRACK_HISTS`.

Momentum and the fit parameters are taken at `TT_Front`, matching `Track_t::PFront()`.

### Normalization

`Ana/data/Norm` carries the branches of `Mu2eEvtAna.cc:843-848` (`Norm_t.hh`), all
`Long64_t`:

| branch | meaning | here |
|---|---|---|
| `ngen` | N(generated events) | from the dataset table, not the ntuple |
| `nntuple` | N(events in the input ntuple) | entries in the files read |
| `nseen` | N(processed events) | equals `nntuple`; every entry is read |
| `naccept` | N(accepted events) | events with >= 1 track surviving set `PRIMARY_SET` (75) |
| `nneg` | N(negative weight events) | 0 — no weights are applied yet |

`ngen` improves on the `FIXME` in `Mu2eEvtAna.cc`, which falls back to the entry count: it
is taken from `datasets.py`'s generated-event count. When only part of the file list is
processed (`--max-files`) it is **scaled** by the fraction of ntuple events seen, and the
scaling is printed so it is not mistaken for an exact count.

With `--split-files` each worker writes its own `Norm` entry and `hadd` concatenates them,
so the merged tree holds one entry per worker and the totals are the sums — the same shape
`Mu2eEvtAna` produces when its jobs are merged.

### Parallel processing and merging

`--jobs N` runs N worker processes over the file list, one partition per file.

By default a worker returns its book as **numpy counts**, which merge by addition, so there
are no temporary files and no `hadd`. `--split-files` instead has each worker write its own
ROOT file and merges them with `scripts/merge_hist_files.sh` (`--clean` removes the parts;
`--keep-parts` keeps them). The two paths were checked to produce bit-identical histograms.

Prefer the default: it avoids the temporary files, and a worker never has to load ROOT.
Both paths were checked to give bit-identical histograms and the same `Norm` totals.

At the end of a run the rate is reported both ways — events per second, and microseconds
per event (in total and per worker), the latter being the figure that stays comparable as
the worker count changes.

## Extending it

- **another selection set** — add to `SELECTIONS` in `selections.py`
- **another histogram** — add a row to `TRACK_HISTS` in `histograms.py`
- **another dataset** — add a `Dataset(...)` to `datasets.py` and drop its file list in
  `file_lists/`
- **another quantity per track** — add a slot to `selections.Track`; it is read once per
  track and shared by the cuts and the histograms

## Note on weights

Every fill currently uses weight 1. `Mu2eEvtAna` applies an event weight
(`evt_.weight_`); `Book.fill` takes a `weight` argument for when that is wanted, but
nothing passes one yet, so these histograms are raw counts.
