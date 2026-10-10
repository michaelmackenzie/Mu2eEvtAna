// Run 1B calorimeter RMC selection: the RMC photon signal stacked on the backgrounds and drawn alone in front, with an
// S/sigma(B) lower pad. The same configuration as Run1BAna/scripts/plotRMCvsBkgFromNtuple.C (v40 datasets), with the
// analysis definitions below and the plotting in the Plotter.
//
// Usage (from the directory holding the Run1BAna.<dataset>.hist files):
//   root.exe -q -b 'Mu2eEvtAna/plotter/examples/run1b_rmc.C(".", "figures/run1b_rmc")'

#include "../Plotter.C"
#include "../physics/Mu2ePhysics.C"

namespace run1b_rmc_config {

  using namespace mu2eplot;
  namespace phys = mu2e_physics;

  //----------------------------------------------------------------------------------------------------
  // Analysis definitions (Run1BAna/analysis/physics.C and Run1BAna/scripts/dataset_info.C, v40)

  const double nmuons_per_pot      = 5.066e-04; // Run 1B v40 stopped muons per POT
  const double poly_nmuons_per_pot = 1.967e-03; // muon stops in the polyethylene (carbon) target per POT
  const double duty_cycle          = 0.322    ; // 4 spills * 107.3 ms / 1333 ms
  const double walltime_week       = 603420.  ;
  const double cosmic_livetime     = 2.29e4   ; // livetime equivalent of the simulated cosmic-ray sample

  struct Sample_t {
    TString dsid;
    double ngen;
    double ndigi;
    double rate; // per POT, per event (pileup), or per second (cosmic rays)
    double emin, emax;
  };

  std::map<TString, Sample_t> samples() {
    using namespace phys;
    const double rate_rmc  = muon_capture_fraction_al*rmc_br_above_57; // samples are normalized above 57 MeV
    const double rate_neut = muon_capture_fraction_al*1.2*5.05038e-06; // N(neutrons per capture)*(KE > 60)*(cz restriction)
    const double rate_prot = muon_capture_fraction_al*0.05*1.10588e-05; // N(protons per capture)*(KE > 60)*(cz restriction)
    const double rate_pgam = (1. - 0.99)/2.*rmc_br_above_57_c*muon_capture_fraction_c; // calorimeter angle restriction
    return {
      {"mnbs", {"mnbs1b1s51r0004", 5000000000, 344196254, 1.                            ,  0.,   1.}},
      {"fgam", {"fgam0b1s51r0004", 1999000000,   1039674, rate_rmc*nmuons_per_pot       , 50., 110.}},
      {"csms", {"csms0b1s51r0004", cosmic_livetime, 2351533, 1.                         ,  0.,   1.}},
      {"pgam", {"pgamcb1s51r0004",  100000000,   5249814, rate_pgam*poly_nmuons_per_pot , 50., 110.}},
      {"neut", {"neut0b1s51r0004",  125000000,    154086, rate_neut*nmuons_per_pot      ,  0.,   1.}},
      {"prot", {"prot0b1s51r0004",  100000000,    158688, rate_prot*nmuons_per_pot      ,  0.,   1.}},
    };
  }

  TString hist_file(const TString& dir, const TString& dsid) { return dir + "/Run1BAna." + dsid + ".hist"; }

  //----------------------------------------------------------------------------------------------------
  Plotter* make_plotter(const TString& hist_dir, const TString& figdir) {
    auto s = samples();

    // Exposure: one week of 1BB running, with N(POT) per event from the signal sample
    TFile* f_sig = TFile::Open(hist_file(hist_dir, s["fgam"].dsid), "READ");
    TH1* h_npot = (f_sig) ? dynamic_cast<TH1*>(f_sig->Get("npot")) : nullptr;
    if(!h_npot) {
      printf("run1b_rmc: N(POT) histogram not found in %s\n", hist_file(hist_dir, s["fgam"].dsid).Data());
      return nullptr;
    }
    const double npot_per_event = h_npot->Integral();
    f_sig->Close();
    Exposure_t exposure;
    exposure.nevents    = walltime_week*duty_cycle/phys::microbunch_period_s;
    exposure.npot       = exposure.nevents*npot_per_event;
    exposure.nmuons     = exposure.npot*nmuons_per_pot;
    exposure.livetime   = walltime_week; // as in plotRMCvsBkgFromNtuple.C, the cosmic rays are scaled by the week's wall time
    exposure.duty_cycle = 1.;            // so the stamp's beam power is N(POT) * E / wall time, as in the original

    Plotter* plotter = new Plotter();
    plotter->figdir_   = figdir;
    plotter->layout_   = Layout_t::run1bana_calo();
    plotter->exposure_ = exposure;
    plotter->stamp_.show_muons_ = false;

    // Processes: signal and pileup shapes are split by set offsets (0: base, 100/200: pileup categories)
    auto add = [&](TString name, TString label, TString key, int offset, bool signal, int color, ExposureKind_t kind) {
      const Sample_t& sample = s[key];
      Process_t p(name, label, (signal) ? kSignal : kBackground, color, hist_file(hist_dir, sample.dsid),
                  mc_norm(sample.rate, sample.ngen, sample.emin, sample.emax), kind);
      p.expected((Long64_t) sample.ndigi).offset(offset);
      if(!signal) p.sys(name, 0.10); // 10% on each background, uncorrelated
      plotter->add_process(p);
    };
    add("rmc"      , "RMC"                , "fgam",   0, true , kGray     , kPOT     );
    add("rmc_pu"   , "RMC"                , "fgam", 100, true , kGray     , kPOT     );
    add("rmc_cpu"  , "RMC"                , "fgam", 200, true , kGray     , kPOT     );
    add("poly"     , "COL5 RMC"           , "pgam",   0, false, 26        , kPOT     );
    add("protons"  , "Protons"            , "prot",   0, false, kAtlantic , kPOT     );
    add("neutrons" , "Neutrons"           , "neut",   0, false, kViolet+6 , kPOT     );
    add("cosmics"  , "Cosmics"            , "csms",   0, false, kGreen-6  , kLivetime);
    add("pileup_lo", "Low pileup clusters", "mnbs",   0, false, kPink     , kNEvents );
    add("pileup_ot", "Other pileup"       , "mnbs", 100, false, kViolet   , kNEvents );
    add("calomu"   , "Calo muon stops"    , "mnbs", 200, false, kOrange   , kNEvents );

    plotter->signal_mode_        = kStacked;
    plotter->signal_front_color_ = kBlue;
    plotter->lower_pad_          = kSignificance;
    plotter->significance_sys_   = true;
    return plotter;
  }
}

//----------------------------------------------------------------------------------------------------
int run1b_rmc(TString hist_dir = ".", TString figdir = "figures/run1b_rmc", std::vector<int> sets = {70, 74, 75}) {
  using namespace mu2eplot;
  Plotter* plotter = run1b_rmc_config::make_plotter(hist_dir, figdir);
  if(!plotter || plotter->init()) return 1;
  int status = 0;
  for(int set : sets) {
    for(bool logy : {false, true}) {
      handle_canvas(plotter->print_stack(plot_t("cluster_energy", "", set, 2,  60.,  100., 1., -1., logy, false, "E", "MeV")), status);
      handle_canvas(plotter->print_stack(plot_t("cluster_time"  , "", set, 5, 600., 1650., 1., -1., logy, false, "t", "ns")), status);
      handle_canvas(plotter->print_stack(plot_t("cluster_radius", "", set, 1, 300.,  700., 1., -1., logy, false, "R", "mm")), status);
    }
    handle_canvas(plotter->print_stack(plot_t("cluster_energy", "", set, 2, 60., 100., 1., -1., false, false, "E", "MeV").normalized()), status);
    plotter->print_yields("cluster_energy", "", set);
  }
  printf("Figures in %s\n", plotter->figdir_.Data());
  delete plotter;
  return status;
}
