#ifndef __MU2EEVTANA_PLOTTER_EXPOSURES__
#define __MU2EEVTANA_PLOTTER_EXPOSURES__
// Exposure presets built from the shared constants. An analysis can also fill an Exposure_t directly.

#include "../PlotTypes.C"
#include "Mu2ePhysics.C"

namespace mu2e_physics {

  // Exposure for a given N(POT) and livetime
  inline mu2eplot::Exposure_t make_exposure(const double npot, const double livetime, const double nmuons_per_pot,
                                            const double duty_cycle = duty_cycle_1bb, const double npot_per_pulse = npot_per_pulse_1bb,
                                            const TString label = "") {
    mu2eplot::Exposure_t e;
    e.npot       = npot;
    e.livetime   = livetime;
    e.nmuons     = npot*nmuons_per_pot;
    e.nevents    = (npot_per_pulse > 0.) ? npot/npot_per_pulse : -1.;
    e.duty_cycle = duty_cycle;
    e.label      = label;
    return e;
  }

  // Run 1A: 28 days of 1BB running
  inline mu2eplot::Exposure_t run1a_exposure() {
    return make_exposure(run1a::npot, run1a::livetime, run1a::nmuons_per_pot, duty_cycle_1bb, npot_per_pulse_1bb, "Run 1A");
  }

  // Run 2: 3e20 POT of 2BB running
  inline mu2eplot::Exposure_t run2_exposure() {
    return make_exposure(run2::npot, run2::livetime, run2::nmuons_per_pot, duty_cycle_2bb, npot_per_pulse_2bb, "Run 2");
  }

  // 1BB running for a given on-spill livetime
  inline mu2eplot::Exposure_t exposure_1bb(const double livetime, const double nmuons_per_pot, const TString label = "") {
    return make_exposure(livetime*npot_rate_1bb, livetime, nmuons_per_pot, duty_cycle_1bb, npot_per_pulse_1bb, label);
  }
}

#endif
