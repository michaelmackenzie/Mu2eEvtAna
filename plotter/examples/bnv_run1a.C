// BNV muon capture signals over the Run 1A EventNtuple backgrounds, from Mu2eEvtAna BNVAna (bnv_ana) histograms.
// The same plots as BNVAna/analysis/make_plots.C, with the analysis definitions below and the plotting in the Plotter.
//
// Usage (from the muse work directory):
//   root.exe -q -b 'Mu2eEvtAna/plotter/examples/bnv_run1a.C({"bnve3b0", "bnvc2b0", "bnvl1b0"}, {201}, "figures/bnv")'
// stacked = true adds the signals to the background stack, drawn alone in front of it, with an S/sigma(B) lower pad.

#include "../Plotter.C"
#include "../physics/Exposures.C"

#include <fstream>
#include <sstream>

namespace bnv_run1a_config {

  using namespace mu2eplot;
  namespace phys = mu2e_physics;

  //----------------------------------------------------------------------------------------------------
  // Analysis definitions

  struct Signal_t {
    TString family;     // BNVAna dataset family
    TString label;
    int     color;
    double  decay_frac; // fraction of generated chi2 decaying inside the generator's fiducial region
  };

  const std::vector<Signal_t> signals = {
    {"bnve0b0", "e^{-}#chi_{1}, m_{1} = 840"  , kBlue    , 1.    },
    {"bnve1b0", "e^{-}#chi_{1}, m_{1} = 860"  , kAzure+2 , 1.    },
    {"bnve2b0", "e^{-}#chi_{1}, m_{1} = 880"  , kCyan+2  , 1.    },
    {"bnve3b0", "e^{-}#chi_{1}, m_{1} = 900"  , kBlue+2  , 1.    },
    {"bnve4b0", "e^{-}#chi_{1}, m_{1} = 915"  , kViolet+1, 1.    },
    {"bnvc0b0", "#chi_{0}#chi_{2}, m_{1} = 840", kOrange+7, 1.    },
    {"bnvc1b0", "#chi_{0}#chi_{2}, m_{1} = 870", kOrange+2, 1.    },
    {"bnvc2b0", "#chi_{0}#chi_{2}, m_{1} = 900", kRed+1   , 1.    },
    {"bnvr0b0", "e^{-}#chi_{2}, m_{1} = 870"  , kPink+7  , 1.    },
    {"bnvr1b0", "e^{-}#chi_{2}, m_{1} = 900"  , kPink+2  , 1.    },
    {"bnvl0b0", "LL #chi_{2}, #tau = 10 ns"   , kTeal+3  , 0.9837},
    {"bnvl1b0", "LL #chi_{2}, #tau = 30 ns"   , kGreen+3 , 0.7952},
    {"bnvl2b0", "LL #chi_{2}, #tau = 100 ns"  , kSpring+5, 0.4244},
    {"bnvl3b0", "LL #chi_{2}, #tau = 300 ns"  , kYellow+3, 0.1811},
    {"bnvl4b0", "LL #chi_{2}, #tau = 1000 ns" , kGray+2  , 0.0595},
  };

  // N(generated) per signal family, from <hist_dir>/bnv_ngen.txt ("<family> <N(gen)>" lines), else the grid production size
  double signal_ngen(const TString& hist_dir, const TString& family) {
    std::ifstream in((hist_dir + "/bnv_ngen.txt").Data());
    std::string line;
    while(std::getline(in, line)) {
      std::istringstream ss(line);
      std::string name; double ngen;
      if(line.empty() || line[0] == '#' || !(ss >> name >> ngen)) continue;
      if(family == name.c_str()) return ngen;
    }
    return 1.e6;
  }

  struct Background_t {
    TString name, label, dsid;
    int color;
    double ngen;     // N(generated), or the generated livetime for cosmic rays
    Long64_t ndigi;  // N(events) in the input sample
    double rate;     // rate per POT (per second for cosmic rays)
    ExposureKind_t exposure;
  };

  // Run 1A EventNtuple backgrounds, as in mumep_ana/analysis/datasets.C
  std::vector<Background_t> backgrounds() {
    using namespace phys;
    const double nmuons_per_pot = run1a::nmuons_per_pot;
    const double rate_dio_95    = muon_decay_fraction_al*nmuons_per_pot*dio_frac_above_95;
    const double rate_rmc       = muon_capture_fraction_al*nmuons_per_pot*rmc_br_above_57;
    const double rate_rmc_0n_80 = rate_rmc*(rmc_ps_0n_r_above_57/rmc_ps_0n_frac_above_57)*rmc_ps_0n_frac_above_80;
    const double rate_rmc_1n_80 = rate_rmc*(rmc_ps_1n_r_above_57/rmc_ps_1n_frac_above_57)*rmc_ps_1n_frac_above_80;
    const double rate_rpc       = run1a::pion_stops_per_pot*rpc_br*run1a::pion_survive_frac;
    const double t_cry          = 4437713.0*5.72/4.438; // generated livetime, with the livetime calculation correction
    return {
      {"pbar"      , "Antiproton"    , "pbar1b1s5r0100", kMagenta-3,  30000000, 6461314, run1a::pbar_stops_per_pot      , kPOT     },
      {"rpc_int"   , "RPC"           , "rpci1b1s5r0100", kGreen-6  , 125000000, 1899806, rate_rpc*rpc_internal_ratio     , kPOT     },
      {"rpc_ext"   , "RPC"           , "rpce1b1s5r0100", kGreen-6  ,     5.e9 ,  458818, rate_rpc                        , kPOT     },
      {"rmc_ext_0n", "RMC (external)", "rmce0b1s5r0100", kRed-7    ,     7.e9 , 4967393, rate_rmc_0n_80                  , kPOT     },
      {"rmc_ext_1n", "RMC (external)", "rmce1b1s5r0100", kRed-7    ,     7.e9 , 2974188, rate_rmc_1n_80                  , kPOT     },
      {"rmc_int_0n", "RMC (internal)", "rmci0b1s5r0100", kRed-9    ,  50000000, 1229006, rate_rmc_0n_80*rmc_internal_ratio, kPOT    },
      {"rmc_int_1n", "RMC (internal)", "rmci1b1s5r0100", kRed-9    ,  50000000,  517029, rate_rmc_1n_80*rmc_internal_ratio, kPOT    },
      {"cosmic"    , "Cosmic ray"    , "cry4ab1s5r0100", kAzure-4  , t_cry    , 4155435, 1.                              , kLivetime},
      {"dio"       , "DIO"           , "dio00b1s5r0100", kMagenta-10,    25e6 , 9368976, rate_dio_95                     , kPOT     },
    };
  }

  TString hist_file(const TString& dir, const TString& dsid) { return dir + "/BNVAna.bnv_ana." + dsid + ".m1.root"; }

  //----------------------------------------------------------------------------------------------------
  // Configure a Plotter for the given signals
  Plotter* make_plotter(const std::vector<TString>& signal_names, const TString& figdir, const double signal_rate,
                        TString hist_dir, const bool stacked) {
    if(hist_dir == "") hist_dir = Form("/exp/mu2e/data/users/%s/BNVAna/histograms", gSystem->Getenv("USER"));
    Plotter* plotter = new Plotter();
    plotter->figdir_   = figdir;
    plotter->layout_   = Layout_t::evtana();
    plotter->exposure_ = phys::run1a_exposure();
    plotter->legend_columns_ = 2;
    plotter->use_offsets_ = false;

    // Signals: the rate R is relative to ordinary muon capture
    const double rate_capture = phys::muon_capture_fraction_al*phys::run1a::nmuons_per_pot;
    for(const auto& name : signal_names) {
      auto sig = std::find_if(signals.begin(), signals.end(), [&](const Signal_t& s) { return s.family == name; });
      if(sig == signals.end()) {
        printf("bnv_run1a: Unknown signal %s\n", name.Data());
        delete plotter;
        return nullptr;
      }
      const double ngen = signal_ngen(hist_dir, sig->family);
      plotter->add_signal(sig->family, sig->label, sig->color, hist_file(hist_dir, sig->family + "s61r0100"),
                          mc_norm(rate_capture*sig->decay_frac, ngen), kPOT, signal_rate, "R");
      printf("  %-8s N(generated) = %.3g\n", sig->family.Data(), ngen);
    }

    // Backgrounds
    for(const auto& bkg : backgrounds()) {
      plotter->add_background(bkg.name, bkg.label, bkg.color, hist_file(hist_dir, bkg.dsid), mc_norm(bkg.rate, bkg.ngen), bkg.exposure)
        .expected(bkg.ndigi);
    }

    if(stacked) {
      plotter->signal_mode_ = kStacked;
      plotter->lower_pad_   = kSignificance;
      plotter->significance_sys_ = false;
    } else {
      plotter->lower_pad_ = kNoPad;
    }
    return plotter;
  }
}

//----------------------------------------------------------------------------------------------------
int bnv_run1a(std::vector<TString> signals = {"bnve3b0", "bnvc2b0", "bnvl1b0"}, std::vector<int> sets = {201},
              TString figdir = "figures/bnv_run1a", double signal_rate = 1.e-13, TString hist_dir = "", bool stacked = false) {
  using namespace mu2eplot;
  Plotter* plotter = bnv_run1a_config::make_plotter(signals, figdir, signal_rate, hist_dir, stacked);
  if(!plotter || plotter->init()) return 1;

  int status = 0;
  for(int set : sets) {
    for(int logy = 0; logy < 2; ++logy) {
      handle_canvas(plotter->print_stack(plot_t("p_wide"     , "trk", set, 2,  75.,  250., 1., -1., logy, false, "p", "MeV/c")), status);
      handle_canvas(plotter->print_stack(plot_t("cosTheta"   , "trk", set, 5,  -1.,    1., 1., -1., logy, false, "cos(#theta)", "")), status);
      handle_canvas(plotter->print_stack(plot_t("t0"         , "trk", set, 5, 400., 1700., 1., -1., logy, false, "t_{0}", "ns")), status);
      handle_canvas(plotter->print_stack(plot_t("d0"         , "trk", set, 5,-200.,  200., 1., -1., logy, false, "D_{0}", "mm")), status);
      handle_canvas(plotter->print_stack(plot_t("rMax"       , "trk", set, 5, 300.,  800., 1., -1., logy, false, "R_{max}", "mm")), status);
      handle_canvas(plotter->print_stack(plot_t("clusterE"   , "trk", set, 4,   0.,  200., 1., -1., logy, false, "Cluster energy", "MeV")), status);
    }
    plotter->print_yields("p_wide", "trk", set, 105., 1.e9);
  }
  printf("Figures in %s\n", plotter->figdir_.Data());
  delete plotter;
  return status;
}
