#ifndef __MU2EEVTANA_PLOTTER_MU2EPHYSICS__
#define __MU2EEVTANA_PLOTTER_MU2EPHYSICS__
// Physics and beam constants shared across Mu2e analyses.
//
// Only quantities that are common to the analyses belong here: muon fates on aluminum, spectrum fractions,
// branching ratios, and beam parameters. Numbers that depend on an analysis' own choices (sample sizes,
// selections, which spectrum model a sample was generated with) stay in that analysis' configuration.
// Campaign normalizations (stopped muons per POT, ...) used by several analyses are grouped by campaign.
//
// Everything is in the mu2e_physics namespace so it can be included alongside an analysis' own physics.C,
// whose global names (muon_capture_fraction_, ...) would otherwise clash.

namespace mu2e_physics {

  //--------------------------------------------------------------------------------------------------
  // Muon fate

  const double muon_capture_fraction_al = 0.609 ; // capture probability for a stopped mu- on aluminum
  const double muon_decay_fraction_al   = 1. - muon_capture_fraction_al; // decay-in-orbit probability on aluminum
  const double muon_capture_fraction_c  = 0.0701; // capture probability on carbon

  //--------------------------------------------------------------------------------------------------
  // Decay-in-orbit spectrum on aluminum: fraction of the DIO spectrum in a momentum range (MeV/c)

  const double dio_frac_above_50 = 8.766e-02;
  const double dio_frac_above_60 = 2.735e-04;
  const double dio_frac_above_80 = 5.671e-08;
  const double dio_frac_above_90 = 7.262e-10;
  const double dio_frac_above_95 = 3.637e-11;
  const double dio_frac_0_60     = 0.999728 ;
  const double dio_frac_60_80    = 2.734e-04;
  const double dio_frac_80_90    = 5.950e-08;

  // Decay-in-orbit on carbon (IPA), fraction above 70 MeV/c
  const double dio_frac_above_70_c = 2.538e-06;

  //--------------------------------------------------------------------------------------------------
  // Radiative pion capture

  const double rpc_br              = 0.0215; // BR(pi- capture --> gamma) on aluminum
  const double rpc_frac_above_50   = 0.9888; // fraction of the RPC photon spectrum above 50 MeV
  const double rpc_internal_ratio  = 0.0069; // internal / external conversion, using hydrogen

  //--------------------------------------------------------------------------------------------------
  // Radiative muon capture

  const double rmc_br_above_57     = 1.40e-5  ; // R(RMC, E > 57 MeV) relative to ordinary muon capture, aluminum
  const double rmc_br_above_57_c   = 1.67e-5  ; // carbon, E > 57 MeV, taken from the TRIUMF oxygen measurement
  const double rmc_frac_above_57   = 0.14800  ; // closure approximation, k_max = 90.1 MeV
  const double rmc_frac_above_85   = 0.0010641; // closure approximation, k_max = 90.1 MeV
  const double rmc_internal_ratio  = 0.0069   ; // internal / external conversion, using hydrogen

  // RMC phase-space model (knock-out channels)
  const double rmc_ps_0n_r_above_57    = 0.099    ; // R(0 knock-outs | E > 57) / R(RMC | E > 57)
  const double rmc_ps_1n_r_above_57    = 0.901    ; // R(1 knock-out  | E > 57) / R(RMC | E > 57)
  const double rmc_ps_0n_frac_above_57 = 0.22887  ;
  const double rmc_ps_1n_frac_above_57 = 0.061620 ;
  const double rmc_ps_0n_frac_above_80 = 0.03319  ;
  const double rmc_ps_1n_frac_above_80 = 0.0013175;

  //--------------------------------------------------------------------------------------------------
  // Beam

  const double proton_kinetic_energy_gev = 8.     ;
  const double microbunch_period_s       = 1.695e-6;
  const double duty_cycle_1bb            = 0.323  ; // taken from Production/JobConfig/ensemble/python/normalizations.py
  const double duty_cycle_2bb            = 0.246  ;
  const double npot_per_pulse_1bb        = 1.58e7 ; // average N(POT) per microbunch
  const double npot_per_pulse_2bb        = 3.93e7 ;
  const double npot_rate_1bb             = npot_per_pulse_1bb / microbunch_period_s; // instantaneous N(POT) per second on-spill
  const double npot_rate_2bb             = npot_per_pulse_2bb / microbunch_period_s;
  const double seconds_per_day           = 24.*60.*60.;

  // Beam power (kW) for a given exposure: proton energy * N(POT) / wall time
  inline double beam_power_kw(const double npot, const double livetime, const double duty_cycle) {
    if(npot <= 0. || livetime <= 0. || duty_cycle <= 0.) return -1.;
    const double energy_j = proton_kinetic_energy_gev*1.602176634e-10;
    return energy_j*npot/(livetime/duty_cycle)/1000.;
  }

  //--------------------------------------------------------------------------------------------------
  // Campaign normalizations shared by several analyses

  namespace run1a { // MDC2025 Run 1A simulation
    const double nmuons_per_pot     = 0.000767114        ; // stopped muons per POT
    const double pion_stops_per_pot = 0.0018801*0.51656  ; // pion stops in the target (infinite lifetime)
    const double pion_survive_frac  = 2393.60487 / 1e10  ; // sum of sampled pion survival weights / N(sampled pions)
    const double pbar_stops_per_pot = 4.7e-18            ; // N(pbar at the stopping target) / POT
    const double ipa_nmuons_per_pot = 2.062e-08          ; // stopped muons in the IPA per POT
    const double livetime           = 28.*seconds_per_day*duty_cycle_1bb; // 28 days of 1BB running at full uptime
    const double npot               = livetime*npot_rate_1bb;
  }

  namespace run2 { // Run 2 (2BB) projections, using the Run 1A stopping rate
    const double nmuons_per_pot = run1a::nmuons_per_pot;
    const double npot           = 3.e20;
    const double livetime       = npot/npot_rate_2bb;
  }
}

#endif
