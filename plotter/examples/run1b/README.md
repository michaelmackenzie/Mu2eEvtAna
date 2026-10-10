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
| `ce` | 80 | CE with pileup (`cele0b1`), R_μe = 10⁻⁸ by default; `ce_sample = "ce_nomix"` for the sample without pileup |
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

- **Sample sizes:** N(generated) and N(events) come from the r0204 entries of `scripts/datasets.C` (EventNtuple
  Run1B-010; pileup Run1B-011).
- **Rates:** v40 Run 1B rates, built from the shared constants in `plotter/physics/Mu2ePhysics.C`.
- **Exposure:** `days` of 1BB running at a 0.322 duty cycle, 6.14×10⁶ POT per microbunch, and 5.066×10⁻⁴ stopped muons
  per POT. Cosmic rays are normalized to the on-spill livetime, pileup (`mnbs`) per microbunch.
- **Partly processed samples:** corrected by N(events)/N(seen) from the `Norm` tree.

Samples without a histogram file are skipped with a message. To make one:

```bash
root.exe -q -b 'Mu2eEvtAna/scripts/make_histograms.C(1, "fgam0b1s51r0204", 1, "run1b_ana", <threads>, <max entries>)'
```

## Physics weights and MC vetoes

`run1b_ana` (`src/Run1BAna.cc`) applies these per event, and stops with an error when an ntuple lacks the information
one needs; it never substitutes an average weight.

- **Flat-spectrum samples** (`fele`, `fgam`, `pgamc`) are generated flat over 50–110 MeV and normalized as
  rate × (110 − 50) / N(gen), so each event is weighted by the physical spectrum's density at the generated energy,
  taken from the ntuple's `primary` branch:
  - DIO (`fele`): Offline's leading-log Al spectrum (`czarnecki_szafron_Al_2016.tbl`).
  - RMC (`fgam`: Al, `pgamc`: C): the Plestid phase-space model, 0 and 1 knock-out channels.
- **RPC** (`rpce`, `rpci`, generated with pion decay off): the pion survival probability, `evtwt.generate`. This needs
  an EventNtuple version that fills `evtwt` from the current event.
- **Beam intensity:** every MC sample is re-weighted from the simulated to the target intensity (log-normal in N(POT)
  per microbunch).
- **Pileup double counting:**
  - In the primary samples, clusters not from the primary particle or its descendants are vetoed, since the pileup
    sample models them (`CaloCluster_t::IsPrimaryAssociated`). This needs an EventNtuple version that fills
    `calomcsim` `prirel`.
  - In the pileup sample, clusters from the RMC, proton, and neutron tails and DIO electrons above 60 MeV are vetoed,
    since the dedicated samples model them.
