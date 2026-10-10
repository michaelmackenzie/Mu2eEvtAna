// Test the Plotter with toy histogram files of known content: normalizations, sampling corrections, control-region
// offsets, label merging, normalization uncertainties, significance, and every drawing mode.
//
// Usage: root.exe -q -b -l 'Mu2eEvtAna/plotter/test/test_plotter.C("<output directory>")'
// Returns the number of failed checks.

#include "../Plotter.C"

using namespace mu2eplot;

namespace plotter_test {

  int nfail_ = 0;

  void check(const bool pass, const TString& what) {
    printf("  [%s] %s\n", (pass) ? "PASS" : "FAIL", what.Data());
    if(!pass) ++nfail_;
  }

  bool close_to(const double a, const double b, const double tol = 1.e-6) {
    return std::fabs(a - b) <= tol*std::max(1., std::max(std::fabs(a), std::fabs(b)));
  }

  // Toy histogram file in the Mu2eEvtAna layout: trk_<set>/p with the given bin contents, and a Norm tree
  void make_file(const TString& name, const std::map<int, std::function<double(double)>>& sets, const std::vector<Long64_t>& nseen) {
    TFile f(name, "RECREATE");
    TDirectory* hist_dir = f.mkdir("Ana")->mkdir("Hist");
    for(const auto& set : sets) {
      hist_dir->mkdir(Form("trk_%i", set.first))->cd();
      TH1F* h = new TH1F("p", "Track momentum", 100, 0., 200.);
      for(int bin = 1; bin <= 100; ++bin) {
        const double val = set.second(h->GetBinCenter(bin));
        h->SetBinContent(bin, val);
        h->SetBinError(bin, std::sqrt(val));
      }
      h->Write();
    }
    f.cd("Ana");
    TDirectory* data_dir = f.GetDirectory("Ana")->mkdir("data");
    data_dir->cd();
    Long64_t n;
    TTree* t = new TTree("Norm", "Normalization information");
    t->Branch("nseen", &n);
    for(auto val : nseen) { n = val; t->Fill(); }
    t->Write();
    f.Close();
  }
}

int test_plotter(TString outdir = "plotter_test") {
  using namespace plotter_test;
  gSystem->mkdir(outdir, true);
  const TString files = outdir + "/files";
  gSystem->mkdir(files, true);

  // Toy inputs
  auto flat  = [](double)   { return 10.; };
  auto falls = [](double x) { return 1000.*std::exp(-x/20.); };
  auto peak  = [](double c) { return [c](double x) { return 100.*std::exp(-0.5*std::pow((x - c)/4., 2)); }; };
  make_file(files + "/bkg_flat.root"  , {{1, flat}}                      , {60, 40}  ); // 100 seen
  make_file(files + "/bkg_falls.root" , {{1, falls}}                     , {50}      );
  make_file(files + "/bkg_falls2.root", {{1, falls}}                     , {50}      );
  make_file(files + "/cosmic.root"    , {{1, flat}, {1001, peak(150.)}}  , {10}      );
  make_file(files + "/sig_a.root"     , {{1, peak(105.)}}                , {10}      );
  make_file(files + "/sig_b.root"     , {{1, peak(130.)}}                , {10}      );
  // Data: the model at half the exposure (half of the data events are seen), with the cosmic shape of the offset set
  make_file(files + "/data.root"      , {{1, [](double x) {
    return std::round(20. + 500.*std::exp(-x/20.) + 2500./(100.*std::sqrt(2.*M_PI)*2.)*100.*std::exp(-0.5*std::pow((x - 150.)/4., 2)));
  }}}, {100});

  // Configuration
  Exposure_t exposure;
  exposure.npot = 1000.; exposure.livetime = 50.; exposure.nmuons = 10.; exposure.duty_cycle = 0.323;
  Plotter plotter;
  plotter.figdir_   = outdir + "/figures";
  plotter.exposure_ = exposure;
  plotter.add_background("flat"  , "Flat"   , kGreen-6 , files + "/bkg_flat.root"  , 2.e-3).expected(200).sys("beam", 0.1);   // x1000 POT x2 sampling = 4
  plotter.add_background("falls" , "Falling", kRed-7   , files + "/bkg_falls.root" , 5.e-4).sys("beam", 0.1);                  // 0.5
  plotter.add_background("falls2", "Falling", kRed-7   , files + "/bkg_falls2.root", 5.e-4).sys("beam", 0.1);                  // 0.5, merged
  plotter.add_background("cosmic", "Cosmic" , kAzure-4 , files + "/cosmic.root"    , 0.1, kLivetime).offset(1000, 1., true).sys("cosmic", 0.2); // x50 = 5
  plotter.add_signal("sig_a", "Signal A", kBlue   , files + "/sig_a.root", 1.e-3, kPOT, 1.e-2, "R");                           // 0.01
  plotter.add_signal("sig_b", "Signal B", kOrange+7, files + "/sig_b.root", 1.e-3, kPOT, 2.e-2, "R");                          // 0.02
  plotter.add_background("missing", "Missing", kGray, files + "/does_not_exist.root", 1.);

  printf("Initialization:\n");
  check(plotter.init() == 0, "init succeeds with a missing background");
  check(plotter.process("missing")->f_ == nullptr, "missing background is skipped");
  check(close_to(plotter.process("flat")->sample_corr_, 2.), "sampling correction N(expected)/N(seen) = 200/100");
  check(plotter.display_label(*plotter.process("sig_a")) == "Signal A [R = 10^{-2}]", "signal label includes the scale");
  check(plotter.display_label(*plotter.process("sig_b")) == "Signal B [R = 2.0 #times 10^{-2}]", "signal label with a mantissa");

  // Expected yields
  printf("Yields:\n");
  const double n_flat  = 100*10.;
  double n_falls = 0.; for(int bin = 1; bin <= 100; ++bin) n_falls += 1000.*std::exp(-(2.*bin - 1.)/20.);
  auto integral = [](TH1* h) { return (h) ? h->Integral(0, h->GetNbinsX()+1) : -1.; };
  auto bkgs = plotter.get_histograms("p", "trk", 1, kBackground);
  check(bkgs.size() == 3, "backgrounds merged by label: Flat, Falling, Cosmic");
  check(close_to(integral(bkgs[0]), 4.*n_flat), "flat background = norm * N(POT) * sampling correction");
  check(close_to(integral(bkgs[1]), 2.*0.5*n_falls), "same-label backgrounds are summed");
  check(close_to(integral(bkgs[2]), 5.*n_flat), "offset set keeps the nominal yield (use_nominal_norm)");
  check(std::fabs(bkgs[2]->GetBinCenter(bkgs[2]->GetMaximumBin()) - 150.) < 2., "offset set gives the shape");
  auto sigs = plotter.get_histograms("p", "trk", 1, kSignal);
  double n_peak = 0.; for(int bin = 1; bin <= 100; ++bin) n_peak += 100.*std::exp(-0.5*std::pow((2.*bin - 1. - 105.)/4., 2));
  check(sigs.size() == 2 && close_to(integral(sigs[0]), 0.01*n_peak, 1.e-4), "signal = norm * N(POT) * scale");
  plotter.use_offsets_ = false;
  TH1* h_cosmic = plotter.get_histogram("p", "trk", 1, kAnyRole, "cosmic");
  check(h_cosmic && h_cosmic->GetMaximumBin() != h_cosmic->FindBin(150.), "use_offsets_ = false reads the nominal set");
  plotter.use_offsets_ = true;

  // Normalization uncertainties: linear within a source, quadrature across sources
  printf("Normalization uncertainties:\n");
  auto sources = plotter.sys_sources("p", "trk", 1, 1, [](const Process_t& p) { return p.is_background(); });
  const double b_beam = 4.*n_flat + 2.*0.5*n_falls, b_cosmic = 5.*n_flat;
  check(sources.size() == 2, "two sources");
  check(close_to(Plotter::sys_in_range(sources, 0, 101), std::sqrt(std::pow(0.1*b_beam, 2) + std::pow(0.2*b_cosmic, 2))),
        "total = sqrt((0.1 beam)^2 + (0.2 cosmic)^2)");
  // Significance per bin and integrated
  TH1* b_total = plotter.get_histogram("p", "trk", 1, kBackground);
  plotter.significance_range_ = kPerBin;
  TH1* z = plotter.significance_hist(sigs[0], b_total, sources);
  const int bin = z->FindBin(105.);
  const double s = sigs[0]->GetBinContent(bin), b = b_total->GetBinContent(bin), sys = Plotter::sys_in_range(sources, bin, bin);
  check(close_to(z->GetBinContent(bin), s/std::sqrt(b + sys*sys)), "per-bin S/sigma(B)");
  plotter.significance_range_ = kAbove;
  z = plotter.significance_hist(sigs[0], b_total, sources);
  const double s_above = sigs[0]->Integral(bin, 101), b_above = b_total->Integral(bin, 101), sys_above = Plotter::sys_in_range(sources, bin, 101);
  check(close_to(z->GetBinContent(bin), s_above/std::sqrt(b_above + sys_above*sys_above)), "integrated (above) S/sigma(B)");
  plotter.significance_range_ = kPerBin;
  Plotter::delete_sources(sources);

  // Data: incomplete data reduces the exposure rather than scaling the data
  printf("Data:\n");
  plotter.add_data("data", files + "/data.root").expected(200);
  check(plotter.init() == 0, "init with data");
  check(close_to(plotter.exposure_used_.npot, 500.) && close_to(plotter.exposure_used_.livetime, 25.), "exposure halved for half the data");
  check(close_to(plotter.exposure_.npot, 1000.), "configured exposure unchanged");
  TH1* h_data = plotter.get_histogram("p", "trk", 1, kData);
  check(close_to(integral(h_data), integral(plotter.get_histogram("p", "trk", 1, kData))) && h_data->GetBinContent(1) > 0., "data is not rescaled");
  check(close_to(integral(plotter.get_histograms("p", "trk", 1, kBackground)[0]), 2.*n_flat), "MC follows the reduced exposure");
  plotter.print_yields("p", "trk", 1, 90., 120.);

  // Drawing modes
  printf("Drawing:\n");
  int status = 0;
  plot_t plot("p", "trk", 1, 2, 0., 200., 1., -1., false, false, "p", "MeV/c");
  plot_t plot_log = plot; plot_log.logy_ = true;

  plotter.signal_mode_ = kOverlay;
  plotter.lower_pad_ = kRatio;
  handle_canvas(plotter.print_stack(plot), status);
  handle_canvas(plotter.print_stack(plot_log), status);
  plotter.lower_pad_ = kDifference;
  plotter.figdir_ = outdir + "/figures/difference"; gSystem->mkdir(plotter.figdir_, true);
  handle_canvas(plotter.print_stack(plot), status);

  plotter.signal_mode_ = kStacked;
  plotter.lower_pad_ = kSignificance;
  plotter.figdir_ = outdir + "/figures/stacked"; gSystem->mkdir(plotter.figdir_, true);
  handle_canvas(plotter.print_stack(plot), status);
  handle_canvas(plotter.print_stack(plot_log), status);
  plotter.stacked_signals_ = {"sig_a"};
  plotter.figdir_ = outdir + "/figures/stacked_one"; gSystem->mkdir(plotter.figdir_, true);
  handle_canvas(plotter.print_stack(plot), status);
  plotter.significance_range_ = kAbove;
  plotter.figdir_ = outdir + "/figures/stacked_one_above"; gSystem->mkdir(plotter.figdir_, true);
  handle_canvas(plotter.print_stack(plot), status);
  plotter.significance_range_ = kPerBin;
  plotter.stacked_signals_.clear();

  plotter.signal_mode_ = kOverlay;
  plotter.lower_pad_ = kSignificance;
  plotter.figdir_ = outdir + "/figures/normalized"; gSystem->mkdir(plotter.figdir_, true);
  handle_canvas(plotter.print_stack(plot_t(plot).normalized()), status);
  plotter.lower_pad_ = kNoPad;
  plotter.figdir_ = outdir + "/figures/nopad"; gSystem->mkdir(plotter.figdir_, true);
  handle_canvas(plotter.print_stack(plot), status);
  handle_canvas(plotter.print_component(plot, "cosmic"), status);
  handle_canvas(plotter.print_roc(plot), status);
  handle_canvas(plotter.print_roc(plot, true, true), status);

  // Signals only (no backgrounds), as for new selection sets
  plotter.lower_pad_ = kRatio;
  plotter.processes_.erase(std::remove_if(plotter.processes_.begin(), plotter.processes_.end(),
                                          [](const Process_t& p) { return p.is_data(); }), plotter.processes_.end());
  for(auto& p : plotter.processes_) if(p.is_background()) p.file_ = "";
  check(plotter.init() == 0, "init with signals only");
  plotter.figdir_ = outdir + "/figures/signals_only"; gSystem->mkdir(plotter.figdir_, true);
  handle_canvas(plotter.print_stack(plot), status);
  check(status == 0, "all canvases drawn");

  // Figures written
  const std::vector<TString> expected = {
    "figures/stack_p_trk_1.png", "figures/stack_p_trk_1_log.png", "figures/difference/stack_p_trk_1.png",
    "figures/stacked/stack_p_trk_1.png", "figures/stacked/stack_p_trk_1_log.png", "figures/stacked_one/stack_p_trk_1.png",
    "figures/stacked_one_above/stack_p_trk_1.png", "figures/normalized/stack_p_trk_1_norm.png", "figures/nopad/stack_p_trk_1.png",
    "figures/nopad/comp_cosmic_p_trk_1.png", "figures/nopad/roc_p_trk_1.png", "figures/nopad/eff_p_trk_1.png",
    "figures/signals_only/stack_p_trk_1.png"
  };
  int nmissing = 0;
  for(const auto& fig : expected) if(gSystem->AccessPathName(outdir + "/" + fig)) { printf("  missing %s\n", fig.Data()); ++nmissing; }
  check(nmissing == 0, "all figures written");

  printf("%i check(s) failed\n", nfail_);
  return nfail_;
}
