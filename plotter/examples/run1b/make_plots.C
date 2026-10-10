// Run 1B stack plots from Mu2eEvtAna Run1BAna (run1b_ana) histograms: the variables each selection cuts on, at each
// step of the RMC, RPC, CE, and proton selections of Mu2eEvtAna/src/Run1BAna.cc, and neutrons on the shared calorimeter
// preselection.
//
// Usage (from the directory holding the Run1BAna.run1b_ana.<dataset>.m1.root files):
//   root.exe -q -b 'Mu2eEvtAna/plotter/examples/run1b/make_plots.C(".", "figures/run1b")'
//   root.exe -q -b 'Mu2eEvtAna/plotter/examples/run1b/make_plots.C(".", "figures/run1b", {"rmc"})'
//
// Each selection's signal is stacked on the other samples and drawn alone in front, with an S/sigma(B) lower pad
// (10% normalization uncertainty on each background). Samples without a histogram file are skipped; make one with
//   root.exe -q -b 'Mu2eEvtAna/scripts/make_histograms.C(1, "<dataset>", 1, "run1b_ana", <threads>, <max entries>)'

#include "samples.C"

namespace run1b {

  //----------------------------------------------------------------------------------------------------
  // Selections: histogram sets (Run1BAna::InitHistSelections), signal, and the samples left out
  struct Selection_t {
    TString              name;
    std::vector<int>     sets;
    std::vector<TString> signals; // stacked together; samples with the same label are summed
    std::vector<TString> exclude; // samples not shown (e.g. the other CE sample)
  };

  inline std::vector<Selection_t> selections(const TString& ce_sample) {
    const TString other_ce = (ce_sample == "ce") ? "ce_nomix" : "ce";
    return {
      // 70: 60 < E < 120 MeV, 500 < t < 1650 ns; 71: + pileup veto (N(cr), E1/E, (E1+E2)/E, t variance, 2nd moment, disk 0);
      // 72: + 500 < r < 580 mm; 73: + no electron line; 74: + no electron line seed, N(electron time cluster hits, z > 1300 mm) < 3
      {"rmc"   , {70, 71, 72, 73, 74}, {"rmc"}             , {"ce", "ce_nomix"}},
      // 90: 60 < E < 140 MeV, 300 < t < 550 ns, pileup and radius vetoes; 93: + no electron line; 94: + no line seed, high-z hits
      {"rpc"   , {90, 93, 94}        , {"rpc", "rpc_int"}  , {"ce", "ce_nomix"}},
      // 80: set 72 with a matched electron line and time cluster, 10 < N(hits) < 40, N(active) > 0, |cos(theta)| > 0.985
      {"ce"    , {80}                , {ce_sample}         , {other_ce}},
      // 40: base selection with a proton line and time cluster; 43: + average hit E(dep) > 2.8 keV; 44: + disk 0, N(cr) < 4, |cos(theta)| > 0.985
      {"proton", {40, 43, 44}        , {"protons"}         , {"ce", "ce_nomix"}},
      // Run1BAna has no neutron selection yet: the calorimeter preselection the other selections share,
      // 1: E > 50 MeV; 2: E > 70 MeV; 70-72 as for the RMC selection
      {"neutron", {1, 2, 70, 71, 72} , {"neutrons"}        , {"ce", "ce_nomix"}},
    };
  }

  //----------------------------------------------------------------------------------------------------
  // Plots for a set: the calorimeter cluster variables cut on, the electron time cluster, and the line for the
  // selections that require one
  inline void plot_set(Plotter& plotter, const int set, const bool has_line, const bool proton, int& status) {
    std::vector<plot_t> plots = {
      //     hist              type   set  rebin  xmin   xmax  ymin ymax logy  logx  xtitle                        unit
      plot_t("energy"        , "cls", set,  4,   50.,  150.,  1., -1., false, false, "Cluster energy"            , "MeV"),
      plot_t("t0"            , "cls", set,  2,  250., 1750.,  1., -1., false, false, "Cluster time"              , "ns" ),
      plot_t("r"             , "cls", set,  1,  350.,  700.,  1., -1., false, false, "Cluster radius"            , "mm" ),
      plot_t("ncr0"          , "cls", set,  1,    0.,   15.,  1., -1., false, false, "N(crystals)"               , ""   ),
      plot_t("fre1"          , "cls", set,  4,  0.3 , 1.05 ,  1., -1., false, false, "E_{1}/E"                   , ""   ),
      plot_t("fre2"          , "cls", set,  4,  0.5 , 1.05 ,  1., -1., false, false, "(E_{1}+E_{2})/E"           , ""   ),
      plot_t("e25_over_e"    , "cls", set,  4,  0.5 , 1.05 ,  1., -1., false, false, "E(5#times5)/E"             , ""   ),
      plot_t("time_var"      , "cls", set,  2,   0. ,   5. ,  1., -1., false, false, "Cluster time variance"     , "ns^{2}"),
      plot_t("second_moment" , "cls", set,  1,   0. , 3000.,  1., -1., false, false, "Cluster second moment"     , ""   ),
      plot_t("disk_id"       , "cls", set,  1,   0. ,   2. ,  1., -1., false, false, "Disk"                      , ""   ),
      plot_t("nhitsabovez"   , "tcs", set,  1,   0. ,  20. ,  1., -1., false, false, "N(time cluster hits, z > 1300 mm)", ""),
      plot_t("nhits"         , "tcs", set,  2,   0. , 100. ,  1., -1., false, false, "N(time cluster hits)"      , ""   ),
    };
    if(proton) plots.push_back(plot_t("avgedep", "tcs", set, 2, 0., 0.005, 1., -1., false, false, "Average hit E(dep)", "MeV"));
    if(has_line) {
      plots.push_back(plot_t("cosTheta", "trk", set, 1, 0.9,  1., 1., -1., false, false, "Line cos(#theta)", ""  )); // the cut is |cos| > 0.985
      plots.push_back(plot_t("nActive" , "trk", set, 1,  0., 80., 1., -1., false, false, "Line N(active hits)", ""));
      plots.push_back(plot_t("t0"      , "trk", set, 4, 250., 1750., 1., -1., false, false, "Line t_{0}", "ns"));
      plots.push_back(plot_t("cos"     , "lns", set, 1, 0.9,  1., 1., -1., false, false, "Line seed cos(#theta)", ""));
    }
    for(auto plot : plots) {
      for(bool logy : {false, true}) {
        plot.logy_ = logy;
        handle_canvas(plotter.print_stack(plot), status);
      }
    }
  }
}

//----------------------------------------------------------------------------------------------------
// hist_dir     : directory of the Run1BAna.run1b_ana.<dataset>.m<mode>.root files
// figdir       : figure directory; each selection goes in <figdir>/<selection>
// names        : selections to plot (rmc, rpc, ce, proton, neutron)
// days         : 1BB running time to normalize to
// ce_sample    : "ce" (mixed with pileup) or "ce_nomix"
// ce_rate      : R(mue) for the CE signal
// stacked      : stack the signal on the other samples (and draw it in front), otherwise overlay it
// wall_time_cosmics: normalize cosmic rays to the wall time, as the older Run1BAna/scripts macros do
int make_plots(TString hist_dir = ".", TString figdir = "figures/run1b",
               std::vector<TString> names = {"rmc", "rpc", "ce", "proton", "neutron"},
               double days = 7., TString ce_sample = "ce_nomix", double ce_rate = 1.e-9,
               bool stacked = true, bool wall_time_cosmics = false, int mode = 1) {
  using namespace mu2eplot;
  int status = 0;
  for(const auto& sel : run1b::selections(ce_sample)) {
    if(std::find(names.begin(), names.end(), sel.name) == names.end()) continue;
    printf("=== Run 1B %s selection: sets", sel.name.Data());
    for(int set : sel.sets) printf(" %i", set);
    printf(", signal");
    for(const auto& sig : sel.signals) printf(" %s", sig.Data());
    printf("\n");

    Plotter plotter;
    plotter.figdir_          = figdir + "/" + sel.name;
    plotter.layout_          = Layout_t::evtana();
    plotter.exposure_        = run1b::exposure(days, run1b::npot_per_event, wall_time_cosmics);
    plotter.require_signals_ = false; // a missing signal sample still gives the background plots
    plotter.use_offsets_     = false;
    plotter.legend_columns_  = 3;
    plotter.signal_mode_     = (stacked) ? kStacked : kOverlay;
    plotter.signal_front_color_ = (stacked) ? kBlue : -1;
    plotter.lower_pad_       = kSignificance;
    plotter.print_missing_hist_summary_ = false;

    // Signals first, then the backgrounds in stack order (the largest backgrounds at the top)
    for(const auto& key : sel.signals) {
      Process_t& p_sig = plotter.add_process(run1b::make_process(*run1b::sample(key), kSignal, hist_dir, mode));
      if(key.BeginsWith("ce")) p_sig.rate_scale(ce_rate, "R_{#mue}");
    }
    for(const auto& s : run1b::samples()) {
      if(std::find(sel.signals.begin(), sel.signals.end(), s.key) != sel.signals.end()) continue;
      if(std::find(sel.exclude.begin(), sel.exclude.end(), s.key) != sel.exclude.end()) continue;
      plotter.add_process(run1b::make_process(s, kBackground, hist_dir, mode));
    }
    if(plotter.init()) { ++status; continue; }

    const bool has_line = sel.name == "ce" || sel.name == "proton";
    for(int set : sel.sets) {
      run1b::plot_set(plotter, set, has_line, sel.name == "proton", status);
      plotter.print_yields("energy", "cls", set);
    }
    printf("Figures in %s\n", plotter.figdir_.Data());
  }
  return status;
}
