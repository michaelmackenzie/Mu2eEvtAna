// Run 1A mu- --> e- and mu- --> e+ stack plots from Mu2eEvtAna ConvAna (cnv_ana) histograms: the track variables of the
// main selections of Mu2eEvtAna/src/ConvAna.cc.
//
// Usage (from the directory holding the ConvAna.cnv_ana.<dataset>.m1.root files):
//   root.exe -q -b 'Mu2eEvtAna/plotter/examples/run1a/make_plots.C(".", "figures/run1a")'
//   root.exe -q -b 'Mu2eEvtAna/plotter/examples/run1a/make_plots.C(".", "figures/run1a", {"mumep"})'
//
// Each selection's signal is drawn over the background stack (stacked = true stacks it and draws it alone in front),
// with an S/sigma(B) lower pad (10% normalization uncertainty on each background). Samples without a histogram file are
// skipped; make one with
//   root.exe -q -b 'Mu2eEvtAna/scripts/make_histograms.C(1, "<dataset>", 1, "cnv_ana", <threads>, <max entries>)'

#include "samples.C"

namespace run1a_ana {

  //----------------------------------------------------------------------------------------------------
  // Selections: histogram sets (ConvAna::InitHistSelections), signal, momentum window, and the samples left out
  struct Selection_t {
    TString              name;
    std::vector<int>     sets;
    TString              signal;
    double               p_min, p_max; // momentum plot range (MeV/c)
    std::vector<TString> exclude;
  };

  inline std::vector<Selection_t> selections() {
    return {
      // 20: e- ID with the upstream and multiple-track vetoes, 100 < p < 110 MeV/c; 60: Run 1A e- ID, t0 > 640 ns, one
      // downstream track; 75: provided e- cut set, 94 < p < 110 MeV/c, t0 > 540 ns
      {"mumem", {20, 60, 75}, "mumem",  95., 110., {"mumep"}},
      // 45: loose e+ selection, 87 < p < 97 MeV/c, triggered; 40: e+ ID, 80 < p < 100 MeV/c, triggered; 42: e+ ID, any p
      {"mumep", {45, 40, 42}, "mumep",  80., 100., {"mumem"}},
    };
  }

  //----------------------------------------------------------------------------------------------------
  // Track variables for a set
  inline void plot_set(Plotter& plotter, const Selection_t& sel, const int set, int& status) {
    std::vector<plot_t> plots = {
      //     hist              type   set  rebin  xmin        xmax      ymin ymax logy  logx  xtitle                       unit
      plot_t("p_2"           , "trk", set,  4,   sel.p_min, sel.p_max,  1., -1., false, false, "p"                         , "MeV/c"),
      plot_t("t0"            , "trk", set,  5,   400.     , 1700.    ,  1., -1., false, false, "t_{0}"                     , "ns"   ),
      plot_t("d0"            , "trk", set,  4,  -150.     ,  150.    ,  1., -1., false, false, "d_{0}"                     , "mm"   ),
      plot_t("rMax"          , "trk", set, 10,   300.     ,  800.    ,  1., -1., false, false, "R_{max}"                   , "mm"   ),
      plot_t("cosTheta"      , "trk", set,  4,     0.     ,    1.    ,  1., -1., false, false, "cos(#theta)"               , ""     ),
      plot_t("nActive"       , "trk", set,  2,     0.     ,  100.    ,  1., -1., false, false, "N(active hits)"            , ""     ),
      plot_t("fitCons_log"   , "trk", set,  4,    -6.     ,    0.    ,  1., -1., false, false, "log_{10} p(#chi^{2})"      , ""     ),
      plot_t("trkQual"       , "trk", set,  4,     0.     ,    1.    ,  1., -1., false, false, "Track quality"             , ""     ),
      plot_t("pid"           , "trk", set,  4,     0.     ,    1.    ,  1., -1., false, false, "PID"                       , ""     ),
      plot_t("clusterE"      , "trk", set,  4,     0.     ,  150.    ,  1., -1., false, false, "Cluster energy"            , "MeV"  ),
      plot_t("ep"            , "trk", set,  2,     0.     ,    1.5   ,  1., -1., false, false, "E/p"                       , ""     ),
      plot_t("dt"            , "trk", set,  2,   -10.     ,   10.    ,  1., -1., false, false, "#Deltat(track - cluster)"  , "ns"   ),
      plot_t("crv_min_deltat", "trk", set,  4,  -250.     ,  250.    ,  1., -1., false, false, "#Deltat(track - CRV)"      , "ns"   ),
    };
    for(auto plot : plots) {
      for(bool logy : {false, true}) {
        plot.logy_ = logy;
        handle_canvas(plotter.print_stack(plot), status);
      }
    }
  }
}

//----------------------------------------------------------------------------------------------------
// hist_dir   : directory of the ConvAna.cnv_ana.<dataset>.m<mode>.root files
// figdir     : figure directory; each selection goes in <figdir>/<selection>
// names      : selections to plot (mumem, mumep)
// rate       : R = Gamma(mu- --> e-+) / Gamma(capture) for the signals
// stacked    : stack the signal on the backgrounds (and draw it in front), otherwise overlay it
int make_plots(TString hist_dir = "/exp/mu2e/data/projects/run1a/mumep_ana/histograms",
               TString figdir = "figures/run1a", std::vector<TString> names = {"mumem", "mumep"},
               double rate = 1.e-13, bool stacked = false, int mode = 1) {
  using namespace mu2eplot;
  int status = 0;
  for(const auto& sel : run1a_ana::selections()) {
    if(std::find(names.begin(), names.end(), sel.name) == names.end()) continue;
    printf("=== Run 1A %s selection: sets", sel.name.Data());
    for(int set : sel.sets) printf(" %i", set);
    printf("\n");

    Plotter plotter;
    plotter.figdir_          = figdir + "/" + sel.name;
    plotter.layout_          = Layout_t::evtana();
    plotter.exposure_        = mu2e_physics::run1a_exposure();
    plotter.require_signals_ = false; // a missing signal sample still gives the background plots
    plotter.use_offsets_     = false;
    plotter.legend_columns_  = 3;
    plotter.signal_mode_     = (stacked) ? kStacked : kOverlay;
    plotter.lower_pad_       = kSignificance;
    plotter.print_missing_hist_summary_ = false;

    // Signal, then the backgrounds in stack order
    plotter.add_process(run1a_ana::make_process(*run1a_ana::sample(sel.signal), kSignal, hist_dir, mode))
      .rate_scale(rate, (sel.signal == "mumep") ? "R_{#mue^{+}}" : "R_{#mue^{-}}");
    for(const auto& s : run1a_ana::samples()) {
      if(s.key == sel.signal) continue;
      if(std::find(sel.exclude.begin(), sel.exclude.end(), s.key) != sel.exclude.end()) continue;
      plotter.add_process(run1a_ana::make_process(s, kBackground, hist_dir, mode));
    }
    if(plotter.init()) { ++status; continue; }

    for(int set : sel.sets) {
      run1a_ana::plot_set(plotter, sel, set, status);
      plotter.print_yields("p_2", "trk", set);
    }
    printf("Figures in %s\n", plotter.figdir_.Data());
  }
  return status;
}
