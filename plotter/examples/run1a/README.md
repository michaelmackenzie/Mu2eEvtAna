# Run 1A μ⁻ → e⁻ and μ⁻ → e⁺ stack plots

Stack plots of the track variables from the Mu2eEvtAna `ConvAna` (`cnv_ana`) histogram files, for the main selections in
`src/ConvAna.cc`.

```bash
# from the directory holding ConvAna.cnv_ana.<dataset>.m1.root
root.exe -q -b 'Mu2eEvtAna/plotter/examples/run1a/make_plots.C(".", "figures/run1a")'
root.exe -q -b 'Mu2eEvtAna/plotter/examples/run1a/make_plots.C(".", "figures/run1a", {"mumep"}, 1.e-13, true)'  # e+ only, signal stacked
```

| selection | sets | signal |
|---|---|---|
| `mumem` | 20 (e⁻ ID, 100–110 MeV/c), 60 (Run 1A e⁻ ID), 75 (provided e⁻ cut set) | μ⁻ → e⁻ (`cele1b1`) |
| `mumep` | 45 (loose e⁺ selection), 40 (e⁺ ID, 80–100 MeV/c), 42 (e⁺ ID, any p) | μ⁻ → e⁺ (`cpos1b1`) |

The backgrounds are DIO, cosmic rays, RPC, RMC (external and internal, 0 and 1 knock-out channels), and antiprotons.
Every set gets these track plots, each in linear and log scale with an S/σ(B) lower pad:
- **Kinematics:** momentum, t0, d0, R(max), and cos θ.
- **Fit quality:** N(active hits) and log₁₀ p(χ²).
- **Classifiers:** track quality and PID.
- **Calorimeter and CRV:** cluster energy, E/p, and the track–cluster and track–CRV time differences.

The S/σ(B) pad assumes a 10% uncorrelated normalization uncertainty on each background. `print_yields` gives the
momentum yields for each set.

## Normalization

`samples.C` holds the Run 1A definitions:

- **Sample sizes:** N(generated) and N(events) come from the r0101 entries of `scripts/datasets.C` (MDC2025au EventNtuples).
- **Rates:** built from the shared constants in `plotter/physics/Mu2ePhysics.C`.
- **Signal rate:** R, the rate relative to muon capture. It is set by the `rate` argument (10⁻¹³ by default) and
  shown in the legend.
- **Exposure:** Run 1A, i.e. 28 days of 1BB running (`run1a_exposure()`).

The antiproton weights are applied by `cnv_ana`. Samples without a histogram file are skipped with a message. To make one:

```bash
root.exe -q -b 'Mu2eEvtAna/scripts/make_histograms.C(1, "cele1b1s5r0101", 1, "cnv_ana", <threads>, <max entries>)'
```
