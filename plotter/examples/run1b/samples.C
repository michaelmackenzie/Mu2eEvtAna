// Run 1B (field-off, straight-line) samples and normalization for the Mu2eEvtAna Run1BAna (run1b_ana) histograms.
//
// N(generated) and N(events) are those of the r0204 entries in Mu2eEvtAna/scripts/datasets.C (EventNtuple Run1B-010/-011).
// Rates are for the v40 Run 1B simulation, using the shared constants in plotter/physics/Mu2ePhysics.C.

#ifndef __MU2EEVTANA_PLOTTER_EXAMPLES_RUN1B_SAMPLES__
#define __MU2EEVTANA_PLOTTER_EXAMPLES_RUN1B_SAMPLES__

#include "../../Plotter.C"
#include "../../physics/Mu2ePhysics.C"

namespace run1b {

  using namespace mu2eplot;
  namespace phys = mu2e_physics;

  //----------------------------------------------------------------------------------------------------
  // Run 1B running conditions (v40 simulation)

  const double nmuons_per_pot      = 5.066e-04; // stopped muons per POT
  const double poly_nmuons_per_pot = 1.967e-03; // muon stops in the polyethylene (carbon) target per POT
  const double duty_cycle          = 0.322    ; // 4 spills * 107.3 ms / 1333 ms
  const double npot_per_event      = 6.14e6   ; // mean N(POT) per microbunch of the 1BB mixing, as used by Run1BAna

  // RPC: pion stops per POT, from the v40 PiBeam/PiTargetStops stages
  const double pibeam_rate         = 11978542. / 1000000000.; // PiBeam efficiency
  const double pi_target_rate      = 8145497. / 250000000.  ; // pion stops in the target per beam pion
  const double pi_select_eff       = 0.999968               ; // pion stop selection efficiency
  const double rpc_per_pot         = pibeam_rate*pi_target_rate*pi_select_eff*phys::rpc_br*phys::rpc_frac_above_50;

  // Exposure: <days> of 1BB running. Cosmic rays are normalized to the on-spill livetime.
  inline Exposure_t exposure(const double days = 7., const double pot_per_event = npot_per_event) {
    Exposure_t e;
    const double walltime = days*phys::seconds_per_day;
    e.nevents    = walltime*duty_cycle/phys::microbunch_period_s;
    e.npot       = e.nevents*pot_per_event;
    e.nmuons     = e.npot*nmuons_per_pot;
    e.livetime   = e.nevents*phys::microbunch_period_s;
    e.duty_cycle = duty_cycle;
    e.label      = Form("Run 1B, %.3g days", days);
    return e;
  }

  //----------------------------------------------------------------------------------------------------
  // Samples
  struct Sample_t {
    TString        key;     // short name, also the process name
    TString        label;   // legend label
    int            color;
    TString        dsid;    // histogram file dataset tag (Mu2eEvtAna/scripts/datasets.C)
    double         ngen;    // N(generated), or the generated livetime for cosmic rays
    Long64_t       ndigi;   // N(events) in the input sample
    double         rate;    // per POT, per event (pileup), or per second (cosmic rays), within [emin, emax]
    double         emin, emax;
    ExposureKind_t exposure;
  };

  inline std::vector<Sample_t> samples() {
    using namespace phys;
    const double rate_ce   = muon_capture_fraction_al*nmuons_per_pot;              // times R(mue), set as the signal scale
    // Flat-spectrum samples: Run1BAna weights fele to the DIO spectrum (from the ntuple's primary branch), so rate * (emax - emin)
    // with the full-spectrum rate is its normalization. The flat photon samples (rmc, poly) are NOT weighted -- the RMC
    // phase-space model weights are not available -- so their shapes and yields are not physical. RPC samples are generated
    // with pion decay off; Run1BAna weights each event by the pion survival probability (evtwt.generate), which turns the
    // pion stop rate below into the physical rate.
    const double rate_rmc  = muon_capture_fraction_al*rmc_br_above_57*nmuons_per_pot; // samples are normalized above 57 MeV
    const double rate_dio  = muon_decay_fraction_al*nmuons_per_pot;                 // all DIO decays; the weights give the spectrum
    const double rate_neut = muon_capture_fraction_al*1.2*5.05038e-06*nmuons_per_pot;  // N(n/capture)*(KE > 60)*(cz restriction)
    const double rate_prot = muon_capture_fraction_al*0.05*1.10588e-05*nmuons_per_pot; // N(p/capture)*(KE > 60)*(cz restriction)
    const double rate_pgam = (1. - 0.99)/2.*rmc_br_above_57_c*muon_capture_fraction_c*poly_nmuons_per_pot; // calo angle restriction
    // In stack order (first at the bottom); the signals (ce, ce_nomix) are not stacked as backgrounds in their own selection
    //        key         label                 color       dsid               N(gen)      N(events)  rate        emin  emax  exposure
    return {
      {"ce"     , "CE"                , kBlue     , "cele0b1s51r0204", 1999000000,   1326786, rate_ce  ,   0.,   1., kPOT     },
      {"ce_nomix", "CE (no pileup)"   , kBlue     , "cele0b0s51r0204", 1900000000,   1184477, rate_ce  ,   0.,   1., kPOT     },
      {"poly"   , "COL5 RMC"          , 26        , "pgamcb1s51r0204",  100000000,   5249814, rate_pgam,  50., 110., kPOT     },
      {"dio"    , "DIO"               , kMagenta-10,"fele0b1s51r0204", 1998000000,    704053, rate_dio ,  50., 110., kPOT     },
      {"rmc"    , "RMC"               , kGray     , "fgam0b1s51r0204", 1999000000,   1039674, rate_rmc ,  50., 110., kPOT     },
      {"rpc"    , "RPC"               , kGreen+2  , "rpce0b1s51r0204", 4999841344,    643155, rpc_per_pot, 0.,   1., kPOT     },
      {"rpc_int", "RPC"               , kGreen+2  , "rpci0b1s51r0204", 4999841344,    230559, rpc_per_pot*rpc_internal_ratio, 0., 1., kPOT},
      {"protons", "Protons"           , kAtlantic , "prot0b1s51r0204",  100000000,    158688, rate_prot,   0.,   1., kPOT     },
      {"neutrons","Neutrons"          , kViolet+6 , "neut0b1s51r0204",  125000000,    154086, rate_neut,   0.,   1., kPOT     },
      {"cosmics", "Cosmics"           , kGreen-6  , "csms0b1s51r0204",     2.29e4,   2351533, 1.       ,   0.,   1., kLivetime},
      {"pileup" , "Pileup"            , kOrange   , "mnbs1b1s51r0204", 5000000000, 344196254, 1.       ,   0.,   1., kNEvents },
    };
  }

  inline const Sample_t* sample(const TString& key) {
    static const std::vector<Sample_t> all = samples();
    for(const auto& s : all) if(s.key == key) return &s;
    return nullptr;
  }

  inline TString hist_file(const TString& dir, const TString& dsid, const int mode = 1) {
    return Form("%s/Run1BAna.run1b_ana.%s.m%i.root", dir.Data(), dsid.Data(), mode);
  }

  // A Process_t for a sample
  inline Process_t make_process(const Sample_t& s, const Role_t role, const TString& hist_dir, const int mode = 1) {
    Process_t p(s.key, s.label, role, s.color, hist_file(hist_dir, s.dsid, mode),
                mc_norm(s.rate, s.ngen, s.emin, s.emax), s.exposure);
    p.expected(s.ndigi);
    if(role == kBackground) p.sys(s.key, 0.10); // 10% on each background, uncorrelated
    return p;
  }
}

#endif
