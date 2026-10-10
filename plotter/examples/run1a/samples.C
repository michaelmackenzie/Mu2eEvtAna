// Run 1A samples and normalization for the Mu2eEvtAna ConvAna (cnv_ana) histograms.
//
// N(generated) and N(events) are those of the r0101 entries in Mu2eEvtAna/scripts/datasets.C (MDC2025au EventNtuples).
// Rates use the shared constants in plotter/physics/Mu2ePhysics.C.

#ifndef __MU2EEVTANA_PLOTTER_EXAMPLES_RUN1A_SAMPLES__
#define __MU2EEVTANA_PLOTTER_EXAMPLES_RUN1A_SAMPLES__

#include "../../Plotter.C"
#include "../../physics/Exposures.C"

namespace run1a_ana {

  using namespace mu2eplot;
  namespace phys = mu2e_physics;

  //----------------------------------------------------------------------------------------------------
  // Samples
  struct Sample_t {
    TString        key;     // short name, also the process name
    TString        label;   // legend label
    int            color;
    TString        dsid;    // histogram file dataset tag (Mu2eEvtAna/scripts/datasets.C)
    double         ngen;    // N(generated), or the generated livetime for cosmic rays
    Long64_t       ndigi;   // N(events) in the input sample
    double         rate;    // per POT (per second for cosmic rays)
    ExposureKind_t exposure;
  };

  inline std::vector<Sample_t> samples() {
    using namespace phys;
    const double nmuons_per_pot = run1a::nmuons_per_pot;
    const double rate_capture   = muon_capture_fraction_al*nmuons_per_pot; // times R, set as the signal scale
    const double rate_dio_95    = muon_decay_fraction_al*nmuons_per_pot*dio_frac_above_95;
    const double rate_rmc       = muon_capture_fraction_al*nmuons_per_pot*rmc_br_above_57;
    const double rate_rmc_0n_80 = rate_rmc*(rmc_ps_0n_r_above_57/rmc_ps_0n_frac_above_57)*rmc_ps_0n_frac_above_80;
    const double rate_rmc_1n_80 = rate_rmc*(rmc_ps_1n_r_above_57/rmc_ps_1n_frac_above_57)*rmc_ps_1n_frac_above_80;
    const double rate_rpc       = run1a::pion_stops_per_pot*rpc_br*run1a::pion_survive_frac;
    const double cosmic_time    = 4437713.0*5.72/4.438; // generated livetime, with a correction to its evaluation
    // In stack order (first at the bottom); each selection's signal is not stacked as a background
    //        key           label                       color        dsid              N(gen)      N(events) rate                               exposure
    return {
      {"mumem"     , "#mu^{-}#rightarrowe^{-}", kBlue      , "cele1b1s5r0101",  10000000, 4141125, rate_capture                     , kPOT     },
      {"mumep"     , "#mu^{-}#rightarrowe^{+}", kBlue      , "cpos1b1s5r0101",  10000000, 3235878, rate_capture                     , kPOT     },
      {"pbar"      , "Antiproton"             , kMagenta-3 , "pbar1b1s5r0101",  30000000, 6461314, run1a::pbar_stops_per_pot        , kPOT     },
      {"rpc_int"   , "RPC"                    , kGreen-6   , "rpci1b1s5r0101", 125000000, 1899806, rate_rpc*rpc_internal_ratio      , kPOT     },
      {"rpc_ext"   , "RPC"                    , kGreen-6   , "rpce1b1s5r0101",    5.e9  ,  458818, rate_rpc                         , kPOT     },
      {"rmc_ext_0n", "RMC (external)"         , kRed-7     , "rmce0b1s5r0101",    7.e9  , 4967393, rate_rmc_0n_80                   , kPOT     },
      {"rmc_ext_1n", "RMC (external)"         , kRed-7     , "rmce1b1s5r0101",    7.e9  , 2974188, rate_rmc_1n_80                   , kPOT     },
      {"rmc_int_0n", "RMC (internal)"         , kRed-9     , "rmci0b1s5r0101",  50000000, 1229006, rate_rmc_0n_80*rmc_internal_ratio, kPOT     },
      {"rmc_int_1n", "RMC (internal)"         , kRed-9     , "rmci1b1s5r0101",  50000000,  517029, rate_rmc_1n_80*rmc_internal_ratio, kPOT     },
      {"cosmic"    , "Cosmic ray"             , kAzure-4   , "cry4ab1s5r0101", cosmic_time, 4155435, 1.                             , kLivetime},
      {"dio"       , "DIO"                    , kMagenta-10, "dio00b1s5r0101",  25000000, 9368976, rate_dio_95                      , kPOT     },
    };
  }

  inline const Sample_t* sample(const TString& key) {
    static const std::vector<Sample_t> all = samples();
    for(const auto& s : all) if(s.key == key) return &s;
    return nullptr;
  }

  inline TString hist_file(const TString& dir, const TString& dsid, const int mode = 1) {
    return Form("%s/ConvAna.cnv_ana.%s.m%i.root", dir.Data(), dsid.Data(), mode);
  }

  // A Process_t for a sample
  inline Process_t make_process(const Sample_t& s, const Role_t role, const TString& hist_dir, const int mode = 1) {
    Process_t p(s.key, s.label, role, s.color, hist_file(hist_dir, s.dsid, mode), mc_norm(s.rate, s.ngen), s.exposure);
    p.expected(s.ndigi);
    if(role == kBackground) p.sys(s.key, 0.10); // 10% on each background, uncorrelated
    return p;
  }
}

#endif
