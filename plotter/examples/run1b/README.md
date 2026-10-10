# Run 1B stack plots

Stack plots from the Mu2eEvtAna `Run1BAna` (`run1b_ana`) histogram files. Each selection in
`src/Run1BAna.cc` is plotted at every step, showing the variables it cuts on.

```bash
# from the directory holding Run1BAna.run1b_ana.<dataset>.m1.root
root.exe -q -b 'Mu2eEvtAna/plotter/examples/run1b/make_plots.C(".", "figures/run1b")'
root.exe -q -b 'Mu2eEvtAna/plotter/examples/run1b/make_plots.C(".", "figures/run1b", {"rmc"}, 30.)'  # RMC only, 30 days
```

| selection | sets | signal (stacked, drawn alone in front) |
|---|---|---|
| `rmc` | 70 → 71 → 72 → 73 → 74 | RMC (`fgam`) |
| `rpc` | 90 → 93 → 94 | RPC, external and internal (`rpce`, `rpci`) |
| `ce` | 80 | CE, no pileup (`cele0b0`) by default; `ce_sample = "ce"` for the mixed sample |
| `proton` | 40 → 43 → 44 | protons (`prot`) |
| `neutron` | 1, 2, 70 → 71 → 72 | neutrons (`neut`); `Run1BAna` has no neutron selection yet, so this is the shared calorimeter preselection |

All other samples are the backgrounds; CE is left out of the non-CE selections. Every set gets these plots:

- **Calorimeter cluster:** energy, time, radius, N(crystals), E1/E, (E1+E2)/E, E(5×5)/E, time variance, second moment, disk.
- **Electron time cluster:** N(hits) and N(hits with z > 1300 mm), the set-74/94 veto variable.
- **CE and proton selections only:** line cos θ (near 1, where the cut is), N(active hits), t0, and line-seed cos θ.
- **Proton selection only:** the time cluster's average hit E(dep).

Each plot is drawn in linear and log scale, with an S/σ(B) lower pad that assumes a 10% uncorrelated uncertainty on each
background. `print_yields` gives the yields for each set.

Figures go to `<figdir>/<selection>/stack_<hist>_<type>_<set>[_log].png`.

## Normalization

`samples.C` holds the Run 1B definitions:

- **Sample sizes:** N(generated) and N(events) come from the r0204 entries of `scripts/datasets.C`.
- **Rates:** the v40 rates of `Run1BAna/scripts/dataset_info.C`, built from the shared constants in `plotter/physics/Mu2ePhysics.C`.
- **Exposure:** `days` of 1BB running at a 0.322 duty cycle, 6.14×10⁶ POT per microbunch, and 5.066×10⁻⁴ stopped muons per POT.
- **Cosmic rays:** normalized to the on-spill livetime. `wall_time_cosmics = true` uses the wall time instead, as the older `Run1BAna/scripts` macros do.
- **Pileup** (`mnbs`): normalized per microbunch event.
- **Partly processed samples:** corrected by N(events)/N(seen) from the `Norm` tree.

Samples without a histogram file are skipped with a message. To make one:

```bash
root.exe -q -b 'Mu2eEvtAna/scripts/make_histograms.C(1, "fgam0b1s51r0204", 1, "run1b_ana", <threads>, <max entries>)'
```

## Physics weights

The flat electron and photon samples (`fele`, `fgam`, `pgamc`) are generated flat over 50–110 MeV (checked in the
FlateMinus generator configuration for `fele`; the photon ranges are those of `Run1BAna/scripts/dataset_info.C`) and
normalized here as rate × (110 − 50) / N(gen). That normalization needs each event weighted by the physical spectrum's
momentum density at the generated momentum. The `run1b_ana` analyzer (`src/Run1BAna.cc`) applies the physics weights
it has the per-event information for, and never substitutes an average weight:

- **DIO (`fele`):** weighted to Offline's leading-log Al DIO spectrum (`czarnecki_szafron_Al_2016.tbl`), using the
  generated electron from the ntuple's `primary` branch. The stored events are filtered toward high momentum, so each
  event needs its own primary; `run1b_ana` stops with an error on ntuples without the branch (versions before Run1B-010).
- **RMC (`fgam`, COL5 `pgamc`):** not weighted. They should use the Plestid RMC phase-space model, whose weights aren't
  available yet, so their shapes and yields in these plots are not physical (COL5 RMC comes out far too large).
- **RPC (`rpce`, `rpci`):** weighted by the pion survival probability exp(−proper time), the `RPCGun` event weight
  stored in the ntuple as `evtwt.generate` (the samples are generated with pion decay off). The current Run1B-010 RPC
  ntuples can't be used: they were made with an EventNtupleMaker that filled `evtwt` from the previous event's weight
  handles (fixed on EventNtuple branch `fix-evtwt-fill-order`), and `run1b_ana` stops on their invalid weights. The RPC
  ntuples need to be remade with the fixed EventNtuple.

## Datasets

The r0204 entries in `scripts/datasets.C` use the EventNtuple Run1B-010 (pileup: Run1B-011) datasets, which have the
`primary` branch. They are made from the same Run1Baw_best_v1_5 reco files as the earlier ntuples, so the event counts
are unchanged. Added in this update:

- RPC internal (`rpci0b1s51r0204`): 230,559 events from 4,999,841,344 generated.
- No-pileup flat electrons (`fele0b0`) and photons (`fgam0b0`), and unbiased pileup (`mnbs0b1`, Run1B-012).
- No-pileup cosmic rays (`csms0b0`), using the same assumed livetime as the mixed sample (2.29×10⁴ s, not evaluated).

The `file_lists/` for the pileup (Run1B-011) and RPC internal ntuples were incomplete and now list every file in the
dataset.
