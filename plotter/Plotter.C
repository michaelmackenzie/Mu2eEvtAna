#ifndef __MU2EEVTANA_PLOTTER_PLOTTER__
#define __MU2EEVTANA_PLOTTER_PLOTTER__
// General signal/background/data histogram plotter.
//
// Stacked backgrounds with any number of signals, each with its own rate, overlaid or added to the stack and drawn alone
// in front of it; data; and a ratio, difference, or S/sqrt(B) lower pad. Also shape-normalized comparisons, component,
// systematic-shift, and ROC plots.
//
// The Plotter holds no physics: each Process_t gives its histogram file, its expected events per MC event per unit
// exposure, which exposure that is (POT, livetime, ...), and any extra scale (e.g. a signal rate). The analysis computes
// these from its own dataset and rate definitions, using the shared constants in physics/Mu2ePhysics.C where they apply.
// See README.md and examples/.

#include "PlotTypes.C"
#include "PlotUtils.C"

#include <algorithm>
#include <functional>

namespace mu2eplot {

class Plotter {
public:

  //-------------------------------------------------------------------------------------------------------
  // Inputs
  std::vector<Process_t> processes_                 ; // processes to plot, in stack order (first at the bottom)
  Layout_t   layout_                                ; // histogram and normalization paths in the input files
  Exposure_t exposure_                              ; // exposure to normalize to
  Exposure_t exposure_used_                         ; // exposure after the data sampling correction, set by init()
  bool       scale_exposure_to_data_ = true         ; // reduce the exposure if the data inputs are incomplete
  bool       require_signals_        = true         ; // missing signal or data files fail init(); missing backgrounds are skipped

  // Output
  TString              figdir_  = "figures/plots"   ; // figure directory
  std::vector<TString> formats_ = {"png"}           ; // figure file formats
  int                  debug_   = 0                 ; // printout level

  // Signals
  SignalMode_t         signal_mode_        = kOverlay; // draw signals over the stack, or stack them
  std::vector<TString> stacked_signals_              ; // in kStacked mode, the signals to stack; empty stacks all signals
  int     signal_fill_        = 0                    ; // fill style of overlaid signals
  int     signal_width_       = 3                    ; // line width of overlaid signals
  int     signal_stack_fill_  = 1001                 ; // fill style of stacked signals
  bool    signal_front_       = true                 ; // draw stacked signals alone in front of the stack
  int     signal_front_style_ = 7                    ; // line style of the signals in front
  int     signal_front_width_ = 4                    ; // line width of the signals in front
  int     signal_front_color_ = -1                   ; // line color of the signals in front; < 0 uses the signal's color

  // Lower pad
  LowerPad_t          lower_pad_          = kRatio   ; // lower pad content
  double              min_ratio_          = 0.5      ; // data / model range
  double              max_ratio_          = 1.5      ;
  bool                draw_sys_band_      = true     ; // add the process normalization uncertainties to the model band
  bool                significance_sys_   = true     ; // include the process normalization uncertainties in sigma(B)
  SignificanceRange_t significance_range_ = kPerBin  ; // per-bin or integrated significance
  double              sys_shift_min_      = 0.       ; // y-range for the systematic shift lower pad
  double              sys_shift_max_      = 2.       ;

  // Style
  bool    use_offsets_       = true                  ; // use control-region set offsets
  int     legend_columns_    = 3                     ;
  double  legend_text_size_  = 0.045                 ;
  int     bkg_line_color_    = kBlack                ; // outline of stacked histograms
  bool    draw_stat_band_    = false                 ; // draw the stack's MC statistical uncertainty band
  int     stat_band_color_   = kGray+1               ; // model statistical uncertainty band
  double  log_span_          = 5.                    ; // maximum orders of magnitude shown on log-y plots
  TString data_label_        = "Data"                ; // data legend label
  Stamp_t stamp_                                     ; // experiment/exposure text
  bool    print_missing_hist_summary_ = true         ; // print missing histograms per plot
  int     missing_hist_path_debug_    = 2            ; // debug level above which example missing paths are printed

  //-------------------------------------------------------------------------------------------------------
  Plotter() {}
  ~Plotter() { close(); }

  //-------------------------------------------------------------------------------------------------------
  // Process definitions. The returned reference is only valid until the next process is added.
  Process_t& add_process(const Process_t& process) {
    for(auto& p : processes_) {
      if(p.name_ != process.name_) continue;
      printf("Plotter::%s: Replacing process %s\n", __func__, process.name_.Data());
      p = process;
      return p;
    }
    processes_.push_back(process);
    return processes_.back();
  }

  Process_t& add_signal(TString name, TString label, int color, TString file, double norm,
                        ExposureKind_t exposure = kPOT, double scale = 1., TString scale_symbol = "") {
    Process_t p(name, label, kSignal, color, file, norm, exposure, scale);
    p.scale_symbol_ = scale_symbol;
    return add_process(p);
  }

  Process_t& add_background(TString name, TString label, int color, TString file, double norm, ExposureKind_t exposure = kPOT) {
    return add_process(Process_t(name, label, kBackground, color, file, norm, exposure));
  }

  Process_t& add_data(TString name, TString file, TString label = "") {
    return add_process(Process_t(name, (label == "") ? data_label_ : label, kData, kBlack, file, 1., kFixed));
  }

  Process_t* process(const TString& name) {
    for(auto& p : processes_) if(p.name_ == name) return &p;
    return nullptr;
  }

  // Set the scale (e.g. the branching fraction) of the given signals, or of all signals
  void set_signal_scale(const double scale, const std::vector<TString>& names = {}) {
    for(auto& p : processes_) {
      if(!p.is_signal()) continue;
      if(!names.empty() && std::find(names.begin(), names.end(), p.name_) == names.end()) continue;
      p.scale_ = scale;
    }
  }

  // Legend label, including the scale if a symbol is given
  TString display_label(const Process_t& p) const {
    if(p.scale_symbol_ == "" || p.scale_ == 1.) return p.label_;
    return p.label_ + " [" + p.scale_symbol_ + " = " + format_sci(p.scale_) + "]";
  }

  // Whether a process is part of the stack
  bool in_stack(const Process_t& p) const {
    if(p.is_background()) return true;
    if(!p.is_signal() || signal_mode_ != kStacked) return false;
    if(stacked_signals_.empty()) return true;
    return std::find(stacked_signals_.begin(), stacked_signals_.end(), p.name_) != stacked_signals_.end();
  }

  // Total scale applied to a process' histograms
  double process_scale(const Process_t& p) const {
    return p.norm_*p.scale_*p.sample_corr_*exposure_used_.value(p.exposure_);
  }

  //-------------------------------------------------------------------------------------------------------
  void close() {
    for(auto& p : processes_) {
      if(p.f_) p.f_->Close();
      p.f_ = nullptr;
    }
  }

  //-------------------------------------------------------------------------------------------------------
  // N(events seen) from the normalization object of a file, or 0 if not found
  double read_norm(TFile* f) const {
    if(!f || layout_.norm_path_ == "") return 0.;
    TObject* o = f->Get(layout_.norm_path_);
    if(!o) {
      printf("Plotter::%s: Normalization object %s not found in %s\n", __func__, layout_.norm_path_.Data(), f->GetName());
      return 0.;
    }
    if(auto t = dynamic_cast<TTree*>(o)) {
      TLeaf* leaf = t->GetLeaf(layout_.norm_branch_);
      if(!leaf) {
        printf("Plotter::%s: Branch %s not found in %s\n", __func__, layout_.norm_branch_.Data(), layout_.norm_path_.Data());
        return 0.;
      }
      double sum = 0.;
      for(Long64_t entry = 0; entry < t->GetEntries(); ++entry) {
        leaf->GetBranch()->GetEntry(entry);
        sum += leaf->GetValue();
      }
      return sum;
    }
    printf("Plotter::%s: %s is not a TTree\n", __func__, layout_.norm_path_.Data());
    return 0.;
  }

  //-------------------------------------------------------------------------------------------------------
  // Open the inputs and evaluate the sampling corrections. Returns 0 on success.
  int init() {
    close();
    gSystem->mkdir(figdir_, true);
    TGaxis::SetMaxDigits(3);
    TGaxis::SetExponentOffset(-0.06, 0.008, "Y");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    exposure_used_ = exposure_;
    double data_fraction = 1.;
    int nopen = 0;
    for(auto& p : processes_) {
      p.sample_corr_ = 1.;
      p.n_seen_ = 0.;
      if(p.exposure_ != kFixed && exposure_.value(p.exposure_) <= 0.) {
        printf("Plotter::%s: Process %s is normalized per unit of an exposure that is not set\n", __func__, p.name_.Data());
        return 1;
      }
      const bool remote = p.file_.Contains("://");
      p.f_ = (p.file_ != "" && (remote || !gSystem->AccessPathName(p.file_))) ? TFile::Open(p.file_, "READ") : nullptr;
      if(p.f_ && p.f_->IsZombie()) { delete p.f_; p.f_ = nullptr; }
      if(!p.f_) {
        printf("Plotter::%s: File not found for %s: %s\n", __func__, p.name_.Data(), p.file_.Data());
        if(p.is_background() || !require_signals_) {
          printf("  --> skipping %s\n", p.name_.Data());
          continue;
        }
        return 1;
      }
      ++nopen;
      if(debug_ > 0) printf("Plotter::%s: Opened %s for %s\n", __func__, p.file_.Data(), p.name_.Data());

      // Correct for unprocessed events
      if(p.n_expected_ <= 0) continue;
      p.n_seen_ = read_norm(p.f_);
      if(p.n_seen_ <= 0.) {
        printf("Plotter::%s: No events seen for %s, no sampling correction applied\n", __func__, p.name_.Data());
        continue;
      }
      const double ratio = p.n_expected_/p.n_seen_;
      if(std::fabs(ratio - 1.) < 1.e-9) continue;
      if(p.is_data()) { // data is not rescaled: the exposure is reduced instead
        if(ratio > 1.) data_fraction = std::min(data_fraction, 1./ratio);
        printf("Plotter::%s: Data %s has %.0f of %lld expected events\n", __func__, p.name_.Data(), p.n_seen_, p.n_expected_);
      } else {
        p.sample_corr_ = ratio;
        printf("Plotter::%s: %s has %.0f of %lld expected events --> scaling by %.4g\n",
               __func__, p.name_.Data(), p.n_seen_, p.n_expected_, ratio);
      }
    }
    if(nopen == 0) {
      printf("Plotter::%s: No input files opened\n", __func__);
      return 1;
    }
    if(data_fraction < 1. && scale_exposure_to_data_) {
      exposure_used_.scale(data_fraction);
      printf("Plotter::%s: Data inputs are incomplete; scaling the exposure by %.4g\n", __func__, data_fraction);
    }
    return 0;
  }

  //-------------------------------------------------------------------------------------------------------
  // A process' histogram for a set, scaled to its expected yield (unless raw). Returns a new, detached histogram.
  TH1* load_hist(const Process_t& p, const TString& hist, const TString& type, const int selection,
                 TString* missing_path = nullptr, const bool raw = false) {
    if(!p.f_) return nullptr;
    int set_offset = (use_offsets_) ? p.set_offset_ : 0;
    if(set_offset > 0 && selection > set_offset && (selection / set_offset) % 2 == 1) set_offset = 0; // already in the control region
    const TString path = layout_.hist_path(hist, type, selection + set_offset);
    TH1* h_in = dynamic_cast<TH1*>(p.f_->Get(path));
    if(!h_in) {
      if(missing_path) *missing_path = path;
      return nullptr;
    }
    TH1* h = (TH1*) h_in->Clone(Form("h_%s_%i", p.name_.Data(), ++uid_));
    h->SetDirectory(0);
    if(raw) return h;

    double scale = process_scale(p);
    if(set_offset > 0) {
      scale *= p.offset_scale_;
      if(p.offset_use_nominal_norm_) {
        TH1* h_nom = dynamic_cast<TH1*>(p.f_->Get(layout_.hist_path(hist, type, selection)));
        const double n_offset  = h->Integral(0, h->GetNbinsX()+1);
        const double n_nominal = (h_nom) ? h_nom->Integral(0, h_nom->GetNbinsX()+1) : 0.;
        if(!h_nom) printf("Plotter::%s: Nominal set %i of %s not found for %s\n", __func__, selection, hist.Data(), p.name_.Data());
        scale *= (n_offset > 0.) ? n_nominal/n_offset : 0.;
      }
    }
    if(debug_ > 6) printf("Plotter::%s: Scaling %s by %.4g\n", __func__, p.name_.Data(), scale);
    h->Scale(scale);
    return h;
  }

  //-------------------------------------------------------------------------------------------------------
  // Apply the role style to a histogram
  void style_hist(TH1* h, const Process_t& p) const {
    h->SetLineColor(p.color_);
    h->SetFillColor(p.color_);
    h->SetMarkerColor(p.color_);
    h->SetLineWidth(2);
    h->SetFillStyle(0);
    if(p.is_background() || (p.is_signal() && in_stack(p))) {
      h->SetFillStyle((p.is_background()) ? 1001 : signal_stack_fill_);
      h->SetLineColor(bkg_line_color_);
      h->SetLineWidth(1);
    } else if(p.is_signal()) {
      h->SetFillStyle(signal_fill_);
      h->SetLineWidth(signal_width_);
    } else if(p.is_data()) {
      h->SetLineColor(kBlack);
      h->SetMarkerColor(kBlack);
      h->SetMarkerStyle(20);
      h->SetMarkerSize(1.);
    }
    if(p.line_color_ >= 0) h->SetLineColor(p.line_color_);
    if(p.line_style_ >= 0) h->SetLineStyle(p.line_style_);
    if(p.line_width_ >= 0) h->SetLineWidth(p.line_width_);
    if(p.fill_style_ >= 0) h->SetFillStyle(p.fill_style_);
  }

  //-------------------------------------------------------------------------------------------------------
  // Scaled histograms of the selected processes, summed by legend label, in process order
  std::vector<TH1*> collect(const TString& hist, const TString& type, const int selection,
                            std::function<bool(const Process_t&)> select) {
    std::map<TString, TH1*> hist_map;
    std::vector<TString> labels;
    std::map<TString, int> missing_counts;
    std::map<TString, TString> missing_examples;
    for(auto& p : processes_) {
      if(!p.f_ || !select(p)) continue;
      TString missing;
      TH1* h = load_hist(p, hist, type, selection, &missing);
      if(!h) {
        ++missing_counts[p.name_];
        if(!missing_examples.count(p.name_)) missing_examples[p.name_] = missing;
        continue;
      }
      const TString label = display_label(p);
      if(hist_map.count(label)) {
        hist_map[label]->Add(h);
        delete h;
        continue;
      }
      labels.push_back(label);
      h->SetName(Form("h_%s_%s_%i_%s_%i", hist.Data(), type.Data(), selection, p.name_.Data(), ++uid_));
      h->SetTitle(label);
      style_hist(h, p);
      hist_map[label] = h;
    }

    if(print_missing_hist_summary_ && !missing_counts.empty()) {
      printf("Plotter::%s: Missing histogram(s) for %s/%s/%i -->", __func__, hist.Data(), type.Data(), selection);
      for(const auto& entry : missing_counts) printf(" %s", entry.first.Data());
      printf("\n");
      if(debug_ > missing_hist_path_debug_) {
        for(const auto& entry : missing_examples) printf("  %s: %s\n", entry.first.Data(), entry.second.Data());
      }
    }

    std::vector<TH1*> histograms;
    for(const auto& label : labels) histograms.push_back(hist_map[label]);
    return histograms;
  }

  //-------------------------------------------------------------------------------------------------------
  // Histograms by role (kAnyRole for all) and optional process-name substring, summed by legend label
  std::vector<TH1*> get_histograms(TString hist, TString type, int selection, int role = kAnyRole, TString name_tag = "") {
    return collect(hist, type, selection, [&](const Process_t& p) {
      return (role == kAnyRole || p.role_ == role) && (name_tag == "" || p.name_.Contains(name_tag));
    });
  }

  //-------------------------------------------------------------------------------------------------------
  // Sum of the histograms by role and optional process-name substring
  TH1* get_histogram(TString hist, TString type, int selection, int role = kAnyRole, TString name_tag = "") {
    auto hists = get_histograms(hist, type, selection, role, name_tag);
    TH1* h = sum_hists(hists, Form("h_%s_%s_%i_%i%s", hist.Data(), type.Data(), selection, role,
                                   (name_tag == "") ? "" : ("_" + name_tag).Data()));
    for(auto hh : hists) delete hh;
    return h;
  }

  //-------------------------------------------------------------------------------------------------------
  // Normalization uncertainty histograms by source (absolute), for the selected processes
  std::map<TString, TH1*> sys_sources(const TString& hist, const TString& type, const int selection, const int rebin,
                                      std::function<bool(const Process_t&)> select) {
    std::map<TString, TH1*> sources;
    for(auto& p : processes_) {
      if(!p.f_ || p.sys_.empty() || !select(p)) continue;
      TH1* h = load_hist(p, hist, type, selection);
      if(!h) continue;
      if(rebin > 1) h->Rebin(rebin);
      for(const auto& sys : p.sys_) {
        if(!sources.count(sys.first)) {
          TH1* hs = (TH1*) h->Clone(Form("sys_%s_%i", sys.first.Data(), ++uid_));
          hs->SetDirectory(0);
          hs->Reset();
          sources[sys.first] = hs;
        }
        sources[sys.first]->Add(h, sys.second);
      }
      delete h;
    }
    return sources;
  }

  // Combined uncertainty over bins [bin_low, bin_high]: sources summed linearly over bins, added in quadrature
  static double sys_in_range(const std::map<TString, TH1*>& sources, const int bin_low, const int bin_high) {
    double var = 0.;
    for(const auto& source : sources) {
      const double val = source.second->Integral(bin_low, bin_high);
      var += val*val;
    }
    return std::sqrt(var);
  }

  static void delete_sources(std::map<TString, TH1*>& sources) {
    for(auto& source : sources) delete source.second;
    sources.clear();
  }

  //-------------------------------------------------------------------------------------------------------
  // Significance S/sqrt(B) or S/sqrt(B + sigma_sys(B)^2), per bin or integrated above/below each bin
  TH1* significance_hist(TH1* h_sig, TH1* h_bkg, const std::map<TString, TH1*>& bkg_sys) {
    if(!h_sig || !h_bkg) return nullptr;
    TH1* h = (TH1*) h_sig->Clone(Form("%s_significance_%i", h_sig->GetName(), ++uid_));
    h->SetDirectory(0);
    h->Reset();
    h->SetFillStyle(0);
    const int nbins = h->GetNbinsX();
    for(int bin = 1; bin <= nbins; ++bin) {
      const int bin_low  = (significance_range_ == kAbove) ? bin : (significance_range_ == kBelow) ? 0         : bin;
      const int bin_high = (significance_range_ == kAbove) ? nbins + 1 : (significance_range_ == kBelow) ? bin : bin;
      const double s = h_sig->Integral(bin_low, bin_high);
      double var = h_bkg->Integral(bin_low, bin_high);
      if(significance_sys_) var += std::pow(sys_in_range(bkg_sys, bin_low, bin_high), 2);
      h->SetBinContent(bin, (var > 0.) ? s/std::sqrt(var) : 0.);
      h->SetBinError  (bin, 0.);
    }
    h->SetLineWidth(3);
    return h;
  }

  //-------------------------------------------------------------------------------------------------------
  // Sum of histograms (nullptr if none)
  TH1* sum_hists(const std::vector<TH1*>& hists, const char* name) {
    TH1* h = nullptr;
    for(auto hh : hists) {
      if(!hh) continue;
      if(!h) {
        h = (TH1*) hh->Clone(Form("%s_%i", name, ++uid_));
        h->SetDirectory(0);
      } else h->Add(hh);
    }
    return h;
  }

  //-------------------------------------------------------------------------------------------------------
  // Canvas layout and axis styles
  void configure_split_canvas(TPad*& pad1, TPad*& pad2, const bool split) {
    pad1 = new TPad(Form("pad1_%i", ++uid_), "pad1", 0., (split) ? 0.3 : 0., 1., 1.0);
    pad2 = new TPad(Form("pad2_%i", ++uid_), "pad2", 0., 0.0, 1., 0.3);
    pad1->Draw();
    if(split) pad2->Draw();
    pad1->SetBottomMargin((split) ? 0.03 : 0.10);
    pad1->SetTopMargin(0.10); pad1->SetLeftMargin(0.12); pad1->SetRightMargin(0.09);
    pad2->SetBottomMargin(0.35); pad2->SetTopMargin(0.04); pad2->SetLeftMargin(pad1->GetLeftMargin()); pad2->SetRightMargin(pad1->GetRightMargin());
    pad1->SetFillColor(0); pad1->SetTickx(1); pad1->SetTicky(1);
    pad2->SetFillColor(0); pad2->SetTickx(1); pad2->SetTicky(1);
  }

  static TString axis_title(const TString& xtitle, const TString& unit) {
    return (unit == "") ? xtitle : xtitle + " (" + unit + ")";
  }

  void style_main_axis(TH1* haxis, const TString& unit, const TString& xtitle, const TString& ytitle,
                       const double scale, const bool split, const bool normalized = false) {
    haxis->SetTitle("");
    haxis->SetXTitle((split) ? TString("") : axis_title(xtitle, unit));
    if(ytitle != "")    haxis->SetYTitle(ytitle);
    else if(normalized) haxis->SetYTitle(Form("Normalized / %.2g %s", haxis->GetBinWidth(1), unit.Data()));
    else                haxis->SetYTitle(Form("Entries / %.2g %s"   , haxis->GetBinWidth(1), unit.Data()));
    haxis->GetYaxis()->SetLabelSize(0.05*scale);
    haxis->GetXaxis()->SetLabelSize((split) ? 0. : 0.05*scale);
    haxis->GetYaxis()->SetTitleSize(0.06*scale);
    haxis->GetYaxis()->SetTitleOffset(0.9);
    haxis->GetXaxis()->SetTitleSize(0.06*scale);
    haxis->GetXaxis()->SetTitleOffset(0.9);
    haxis->GetXaxis()->SetTitleFont(132);
    haxis->GetYaxis()->SetTitleFont(132);
    haxis->GetXaxis()->SetLabelFont(132);
    haxis->GetYaxis()->SetLabelFont(132);
  }

  void style_ratio_axis(TH1* haxis_r, const TString& xtitle, const TString& unit) {
    haxis_r->SetXTitle(axis_title(xtitle, unit));
    haxis_r->GetYaxis()->SetNdivisions(507);
    haxis_r->GetYaxis()->SetLabelSize(0.125);
    haxis_r->GetYaxis()->SetLabelOffset(0.01);
    haxis_r->GetYaxis()->SetTitleSize(0.145);
    haxis_r->GetYaxis()->SetTitleOffset(0.35);
    haxis_r->GetXaxis()->SetLabelSize(0.15);
    haxis_r->GetXaxis()->SetLabelOffset(0.008);
    haxis_r->GetXaxis()->SetTitleSize(0.18);
    haxis_r->GetXaxis()->SetTitleOffset(0.8);
    haxis_r->GetXaxis()->SetTitleFont(132);
    haxis_r->GetYaxis()->SetTitleFont(132);
    haxis_r->GetXaxis()->SetLabelFont(132);
    haxis_r->GetYaxis()->SetLabelFont(132);
  }

  static void resolve_x_range(const TH1* h, double& xmin, double& xmax) {
    if(!h) return;
    const double abs_xmin = h->GetXaxis()->GetXmin();
    const double abs_xmax = h->GetXaxis()->GetXmax();
    if(xmin >= xmax || xmax < abs_xmin || xmin > abs_xmax) {
      xmin = abs_xmin;
      xmax = abs_xmax;
    }
    xmin = std::max(xmin, abs_xmin);
    xmax = std::min(xmax, abs_xmax);
  }

  //-------------------------------------------------------------------------------------------------------
  // Legends
  struct LegendEntry_t {
    TObject* obj_;
    TString  label_;
    TString  option_;
  };

  // Approximate rendered length of a TLatex label, in characters. Sub/superscripts count as smaller characters.
  static double latex_length(const TString& label) {
    double length = 0.;
    int script_depth = 0; // brace depth inside a sub/superscript
    bool script_next = false;
    for(int index = 0; index < label.Length(); ++index) {
      const char ch = label[index];
      const double weight = (script_depth > 0 || script_next) ? 0.6 : 1.;
      if(ch == '^' || ch == '_') { script_next = true; continue; }
      if(ch == '{') { if(script_next || script_depth > 0) ++script_depth; script_next = false; continue; }
      if(ch == '}') { if(script_depth > 0) --script_depth; continue; }
      if(ch == '#') { // a symbol: one character
        while(index + 1 < label.Length() && std::isalpha(label[index+1])) ++index;
      }
      length += weight;
      script_next = false;
    }
    return length;
  }

  // Legend across the top of the pad, with the column count (at most legend_columns_) and text size chosen so the
  // labels fit. y1 is set to the legend's lower edge.
  TLegend* make_legend(TPad* pad, const std::vector<LegendEntry_t>& entries, const double text_size, double& y1) {
    y1 = 1. - pad->GetTopMargin();
    if(entries.empty()) return nullptr;
    const double x1 = pad->GetLeftMargin() + 0.03, x2 = 1. - pad->GetRightMargin() - 0.02;
    const double aspect = (pad->GetWh()*pad->GetAbsHNDC())/(pad->GetWw()*pad->GetAbsWNDC()); // pad height / width in pixels
    double max_length = 1.;
    for(const auto& entry : entries) max_length = std::max(max_length, latex_length(entry.label_));
    // Text size at which the longest label fills the text part (75%) of a column, for characters ~0.42 of the text size wide
    auto fit_size = [&](int ncol) { return 0.75*(x2 - x1)/ncol/(0.42*aspect*max_length); };
    int ncol = std::max(1, std::min(legend_columns_, (int) entries.size()));
    for(; ncol > 1; --ncol) if(fit_size(ncol) >= 0.75*text_size) break;
    const double size = std::min(text_size, fit_size(ncol));
    const int nrows = (entries.size() + ncol - 1)/ncol;
    const double y2 = 1. - pad->GetTopMargin() - 0.03;
    y1 = y2 - nrows*size*1.35 - 0.01;
    TLegend* leg = new TLegend(x1, y1, x2, y2);
    leg->SetNColumns(ncol);
    leg->SetLineWidth(0); leg->SetLineColor(0); leg->SetFillColor(0); leg->SetFillStyle(0);
    leg->SetTextSize(size);
    leg->SetTextFont(132);
    for(const auto& entry : entries) leg->AddEntry(entry.obj_, entry.label_, entry.option_);
    leg->Draw();
    return leg;
  }

  // y-range leaving room for a legend covering the top legend_frac of the frame
  void resolve_y_range(double& ymin, double& ymax, double max_val, double min_val, const bool logy, const double legend_frac) {
    if(ymin < ymax) return;
    const double room = std::max(0.2, 1. - legend_frac);
    if(max_val <= 0.) max_val = 1.;
    if(!logy) {
      ymin = 0.;
      ymax = max_val/room;
      return;
    }
    if(min_val <= 0. || min_val >= max_val) min_val = 1.e-3*max_val;
    min_val = std::max(min_val, std::pow(10., -log_span_)*max_val);
    const double orders = std::log10(max_val/min_val);
    ymin = std::pow(10., -0.05*orders)*min_val; // ~5% buffer on the bottom
    ymax = ymin*std::pow(max_val/ymin, 1./room);
  }

  //-------------------------------------------------------------------------------------------------------
  // Stacked backgrounds (and signals, in kStacked mode) with signals and data overlaid
  TCanvas* plot_stack(plot_t plot) {
    const TString hist = plot.hist_, type = plot.type_;
    const int selection = plot.selection_, rebin = plot.rebin_;
    double xmin = plot.xmin_, xmax = plot.xmax_, ymin = plot.ymin_, ymax = plot.ymax_;
    if(debug_ > 0) printf("Plotter::%s: Plotting %s/%s/%i\n", __func__, hist.Data(), type.Data(), selection);

    // Retrieve the inputs
    auto bkgs   = collect(hist, type, selection, [&](const Process_t& p) { return p.is_background(); });
    auto s_sigs = collect(hist, type, selection, [&](const Process_t& p) { return p.is_signal() &&  in_stack(p); });
    auto o_sigs = collect(hist, type, selection, [&](const Process_t& p) { return p.is_signal() && !in_stack(p); });
    auto datas  = collect(hist, type, selection, [&](const Process_t& p) { return p.is_data(); });
    if(bkgs.empty() && s_sigs.empty() && o_sigs.empty() && datas.empty()) {
      printf("Plotter::%s: No histograms found for %s/%s/%i\n", __func__, hist.Data(), type.Data(), selection);
      return nullptr;
    }
    TH1* data = sum_hists(datas, Form("data_%s_%s_%i", hist.Data(), type.Data(), selection));
    for(auto h : datas) delete h;
    if(rebin > 1) {
      for(auto h : bkgs  ) h->Rebin(rebin);
      for(auto h : s_sigs) h->Rebin(rebin);
      for(auto h : o_sigs) h->Rebin(rebin);
      if(data) data->Rebin(rebin);
    }
    std::vector<TH1*> sigs(s_sigs);
    sigs.insert(sigs.end(), o_sigs.begin(), o_sigs.end());
    TH1* bkg_total = sum_hists(bkgs, Form("bkg_%s_%s_%i", hist.Data(), type.Data(), selection));
    std::vector<TH1*> stacked(bkgs);
    stacked.insert(stacked.end(), s_sigs.begin(), s_sigs.end());
    TH1* stack_total = sum_hists(stacked, Form("stack_%s_%s_%i", hist.Data(), type.Data(), selection));

    // Decide on the lower pad
    const bool has_ratio = (lower_pad_ == kRatio || lower_pad_ == kDifference) && ((data && stack_total) || (!sigs.empty() && bkg_total));
    const bool has_sig   = lower_pad_ == kSignificance && !sigs.empty() && bkg_total;
    const bool split     = has_ratio || has_sig;
    const float scale    = (split) ? 1. : 0.75; // text scale

    // Lower pad inputs, evaluated before any shape normalization
    auto model_sys = (draw_sys_band_ && has_ratio)
      ? sys_sources(hist, type, selection, rebin, [&](const Process_t& p) { return in_stack(p); }) : std::map<TString, TH1*>{};
    std::vector<TH1*> significances;
    if(has_sig) {
      auto bkg_sys = (significance_sys_) ? sys_sources(hist, type, selection, rebin, [&](const Process_t& p) { return p.is_background(); })
        : std::map<TString, TH1*>{};
      for(auto sig : sigs) {
        TH1* h = significance_hist(sig, bkg_total, bkg_sys);
        const bool stacked_sig = std::find(s_sigs.begin(), s_sigs.end(), sig) != s_sigs.end();
        h->SetLineColor((stacked_sig && signal_front_color_ >= 0) ? signal_front_color_ : sig->GetFillColor());
        significances.push_back(h);
      }
      delete_sources(bkg_sys);
    }

    // Shape comparison: unit-area signals, background total, and data
    resolve_x_range((data) ? data : (bkg_total) ? bkg_total : sigs[0], xmin, xmax);
    if(plot.normalize_) {
      const double n_bkg = integral_in_range(bkg_total, xmin, xmax);
      if(n_bkg > 0.) {
        for(auto h : bkgs) h->Scale(1./n_bkg);
        bkg_total->Scale(1./n_bkg);
        for(auto& source : model_sys) source.second->Scale(1./n_bkg);
      }
      for(auto h : sigs) {
        const double n = integral_in_range(h, xmin, xmax);
        if(n > 0.) h->Scale(1./n);
      }
      if(data) {
        const double n = integral_in_range(data, xmin, xmax);
        if(n > 0.) data->Scale(1./n);
      }
      delete stack_total;
      stack_total = sum_hists(stacked, Form("stack_%s_%s_%i", hist.Data(), type.Data(), selection));
    }

    // Build the stack
    THStack* stack = new THStack(Form("s_%s_%s_%i_%i", hist.Data(), type.Data(), selection, ++uid_), "Stack");
    double min_val = 1.e30;
    for(auto h : stacked) {
      stack->Add(h);
      min_val = std::min(min_val, min_in_range(h, xmin, xmax, false, 1.e-20));
    }
    for(auto h : o_sigs) min_val = std::min(min_val, min_in_range(h, xmin, xmax, false, 1.e-20));
    if(data) min_val = std::min(min_val, min_in_range(data, xmin, xmax, false, 1.e-20));

    // Canvas and axis
    TCanvas* c = new TCanvas(Form("c_%s_%s_%i_%i", hist.Data(), type.Data(), selection, ++uid_), "Canvas", 1200, 1000);
    TPad *pad1(nullptr), *pad2(nullptr);
    configure_split_canvas(pad1, pad2, split);
    pad1->cd();
    TH1* first = (data) ? data : (stack_total) ? stack_total : sigs[0];
    TH1* haxis = (TH1*) first->Clone(Form("axis_%s_%s_%i_%i", hist.Data(), type.Data(), selection, ++uid_));
    haxis->SetDirectory(0);
    haxis->Reset();
    haxis->SetLineWidth(0);
    haxis->SetFillStyle(0);
    style_main_axis(haxis, plot.unit_, plot.xtitle_, plot.ytitle_, scale, split, plot.normalize_);
    haxis->Draw("hist");
    haxis->GetXaxis()->SetRangeUser(xmin, xmax);

    // Draw the stack, and optionally its statistical uncertainty
    double max_val = 0.;
    TH1* stat_band = nullptr;
    if(stack_total) {
      stack->Draw("hist noclear same");
      max_val = max_in_range(stack_total, xmin, xmax);
    }
    if(stack_total && draw_stat_band_) {
      stat_band = (TH1*) stack_total->Clone(Form("stat_%s_%s_%i_%i", hist.Data(), type.Data(), selection, ++uid_));
      stat_band->SetDirectory(0);
      stat_band->SetFillStyle(3001);
      stat_band->SetFillColor(stat_band_color_);
      stat_band->SetLineColor(0);
      stat_band->SetLineWidth(0);
      stat_band->SetMarkerSize(0);
      stat_band->Draw("same E2");
    }

    // Draw the signals over the stack
    for(auto h : o_sigs) {
      h->Draw("hist same");
      max_val = std::max(max_val, max_in_range(h, xmin, xmax));
    }
    std::vector<TH1*> fronts;
    if(signal_front_) {
      for(auto h : s_sigs) {
        TH1* front = (TH1*) h->Clone(Form("%s_front", h->GetName()));
        front->SetDirectory(0);
        front->SetFillStyle(0);
        front->SetLineColor((signal_front_color_ >= 0) ? signal_front_color_ : h->GetFillColor());
        front->SetLineStyle(signal_front_style_);
        front->SetLineWidth(signal_front_width_);
        front->Draw("hist same");
        fronts.push_back(front);
        max_val = std::max(max_val, max_in_range(front, xmin, xmax));
      }
    }

    // Draw the data
    if(data) {
      data->Draw("EX0 same");
      max_val = std::max(max_val, max_in_range(data, xmin, xmax, true));
    }

    // Legend: data, overlaid signals, then the stack from the top down. A stacked signal drawn in front gets one entry
    // showing both its stack fill and its front outline.
    std::vector<LegendEntry_t> entries;
    if(data) entries.push_back({data, data_label_, "PL"});
    for(auto h : o_sigs) entries.push_back({h, h->GetTitle(), (h->GetFillStyle() == 0) ? "L" : "F"});
    for(int index = (int) stacked.size() - 1; index >= 0; --index) {
      TH1* h = stacked[index];
      const auto front = std::find_if(fronts.begin(), fronts.end(), [&](TH1* f) { return TString(f->GetName()) == TString(h->GetName()) + "_front"; });
      if(front == fronts.end()) {
        entries.push_back({h, h->GetTitle(), "F"});
        continue;
      }
      TH1* proxy = (TH1*) h->Clone(Form("%s_legend", h->GetName())); // legend-only copy, owned by the legend entry list
      proxy->SetDirectory(0);
      proxy->SetLineColor((*front)->GetLineColor());
      proxy->SetLineStyle((*front)->GetLineStyle());
      proxy->SetLineWidth(std::max(2, (*front)->GetLineWidth()/2));
      entries.push_back({proxy, h->GetTitle(), "FL"});
    }
    double leg_y1 = 1. - pad1->GetTopMargin();
    TLegend* leg = make_legend(pad1, entries, legend_text_size_*scale, leg_y1);

    // Axis ranges, leaving room for the legend
    const double frame_height = 1. - pad1->GetTopMargin() - pad1->GetBottomMargin();
    const double legend_frac  = (entries.empty()) ? 0. : (1. - pad1->GetTopMargin() - leg_y1)/frame_height + 0.02;
    resolve_y_range(ymin, ymax, max_val, min_val, plot.logy_, legend_frac);
    if(debug_ > 1) printf("Plotter::%s: x = [%g, %g], y = [%g, %g]\n", __func__, xmin, xmax, ymin, ymax);
    haxis->GetYaxis()->SetRangeUser(ymin, ymax);
    haxis->GetXaxis()->SetRangeUser(xmin, xmax);
    if(plot.logy_) pad1->SetLogy();
    if(plot.logx_) pad1->SetLogx();

    // Lower pad
    if(split) {
      pad2->cd();
      TH1* haxis_r = (TH1*) haxis->Clone(Form("axis_r_%s_%s_%i_%i", hist.Data(), type.Data(), selection, ++uid_));
      haxis_r->SetDirectory(0);
      haxis_r->Draw("hist");
      style_ratio_axis(haxis_r, plot.xtitle_, plot.unit_);
      haxis_r->GetXaxis()->SetRangeUser(xmin, xmax);
      if(plot.logx_) pad2->SetLogx();
      if(has_sig) draw_significance(haxis_r, significances, xmin, xmax);
      else        draw_comparison(haxis_r, data, sigs, (data) ? stack_total : bkg_total, model_sys, !s_sigs.empty(), xmin, xmax);
    }
    delete_sources(model_sys);

    pad1->cd();
    stamp_.draw(exposure_used_, scale);
    pad1->RedrawAxis();
    if(split) pad2->RedrawAxis();
    return c;
  }

  //-------------------------------------------------------------------------------------------------------
  // Lower pad: significance of each signal
  void draw_significance(TH1* haxis_r, const std::vector<TH1*>& significances, const double xmin, const double xmax) {
    double max_val = 0.;
    for(auto h : significances) {
      h->Draw("hist same");
      max_val = std::max(max_val, max_in_range(h, xmin, xmax));
    }
    haxis_r->GetYaxis()->SetRangeUser(0., (max_val > 0.) ? 1.3*max_val : 1.);
    haxis_r->GetYaxis()->SetNoExponent(true); // the exponent label would be clipped by the pad edge
    haxis_r->SetYTitle((significance_sys_) ? "S/#sigma_{B}" : "S/#sqrt{B}");
  }

  //-------------------------------------------------------------------------------------------------------
  // Lower pad: data vs. the model with its uncertainty band, or each signal vs. the background without data
  void draw_comparison(TH1* haxis_r, TH1* data, const std::vector<TH1*>& sigs, TH1* model,
                       const std::map<TString, TH1*>& model_sys, const bool stacked_signal, const double xmin, const double xmax) {
    const bool diff = lower_pad_ == kDifference;
    std::vector<TH1*> nums;
    if(data) {
      TH1* h = (TH1*) data->Clone(Form("%s_r", data->GetName()));
      h->SetDirectory(0);
      nums.push_back(h);
    } else {
      for(auto sig : sigs) {
        TH1* h = (TH1*) sig->Clone(Form("%s_r", sig->GetName()));
        h->SetDirectory(0);
        h->SetFillStyle(0);
        h->SetLineColor(sig->GetFillColor());
        h->SetLineStyle(kSolid);
        h->SetLineWidth(2);
        nums.push_back(h);
      }
    }
    if(nums.empty() || !model) return;

    // Uncertainty band on the model, statistical and statistical + normalization
    TGraphErrors *g_stat(nullptr), *g_tot(nullptr);
    if(data) {
      int bin_low, bin_high; bin_range(model, xmin, xmax, bin_low, bin_high);
      const int n = bin_high - bin_low + 1;
      std::vector<double> x(n), y(n), xerr(n), stat(n), tot(n);
      for(int bin = bin_low; bin <= bin_high; ++bin) {
        const int index = bin - bin_low;
        const double val = model->GetBinContent(bin);
        x[index] = model->GetBinCenter(bin);
        xerr[index] = model->GetBinWidth(bin)/2.;
        y[index] = (diff) ? 0. : 1.;
        if(val <= 0.) continue;
        const double norm = (diff) ? 1. : val;
        stat[index] = model->GetBinError(bin)/norm;
        tot [index] = std::sqrt(std::pow(stat[index], 2) + std::pow(sys_in_range(model_sys, bin, bin)/norm, 2));
      }
      g_tot = new TGraphErrors(n, x.data(), y.data(), xerr.data(), tot.data());
      g_tot->SetFillStyle(3004);
      g_tot->SetFillColor(stat_band_color_);
      g_tot->SetLineWidth(0);
      g_stat = new TGraphErrors(n, x.data(), y.data(), xerr.data(), stat.data());
      g_stat->SetFillStyle(3001);
      g_stat->SetFillColor(stat_band_color_);
      g_stat->SetLineWidth(0);
    }

    double max_val(-1.e30), min_val(1.e30);
    for(auto h : nums) {
      if(diff) h->Add(model, -1.);
      else     h->Divide(model);
      max_val = std::max(max_val, max_in_range(h, xmin, xmax, (bool) data));
      min_val = std::min(min_val, min_in_range(h, xmin, xmax, (bool) data, -1.e30));
    }
    const TString num_title = (data) ? data_label_ : TString("Signal");
    const TString den_title = (data && stacked_signal) ? "Total" : "Bkg";
    haxis_r->SetYTitle(num_title + ((diff) ? " - " : " / ") + den_title);
    if(diff) {
      const double buffer = 0.05*(max_val - min_val);
      haxis_r->GetYaxis()->SetRangeUser(min_val - buffer, max_val + buffer);
    } else if(data) {
      haxis_r->GetYaxis()->SetRangeUser(min_ratio_, max_ratio_);
    } else {
      haxis_r->GetYaxis()->SetRangeUser(0., (max_val > 0.) ? 1.1*max_val : 1.);
    }

    if(g_tot ) g_tot ->Draw("E2");
    if(g_stat) g_stat->Draw("E2");
    if(data) {
      TLine* line = new TLine(xmin, (diff) ? 0. : 1., xmax, (diff) ? 0. : 1.);
      line->SetLineWidth(2);
      line->SetLineColor(kBlack);
      line->SetLineStyle(kDashed);
      line->Draw("same");
    }
    for(auto h : nums) h->Draw((data) ? "EX0 same" : "hist same");
  }

  //-------------------------------------------------------------------------------------------------------
  // One process (or group of processes matching the name tag) alone
  TCanvas* plot_component(plot_t plot, TString process) {
    const TString hist = plot.hist_, type = plot.type_;
    const int selection = plot.selection_;
    double xmin = plot.xmin_, xmax = plot.xmax_, ymin = plot.ymin_, ymax = plot.ymax_;
    if(debug_ > 0) printf("Plotter::%s: Plotting %s/%s/%i for %s\n", __func__, hist.Data(), type.Data(), selection, process.Data());

    TH1* h = get_histogram(hist, type, selection, kAnyRole, process);
    if(!h) return nullptr;
    if(plot.rebin_ > 1) h->Rebin(plot.rebin_);
    h->SetName(Form("h_%s_%s_%s_%i_%i", process.Data(), hist.Data(), type.Data(), selection, ++uid_));

    TCanvas* c = new TCanvas(Form("c_%s_%s_%s_%i_%i", process.Data(), hist.Data(), type.Data(), selection, ++uid_), "Canvas", 1200, 700);
    c->SetBottomMargin(0.13); c->SetTopMargin(0.10); c->SetLeftMargin(0.12); c->SetRightMargin(0.09);
    c->SetFillColor(0); c->SetTickx(1); c->SetTicky(1);
    resolve_x_range(h, xmin, xmax);
    style_main_axis(h, plot.unit_, plot.xtitle_, plot.ytitle_, 1., false);
    h->GetXaxis()->SetTitleSize(0.06);
    h->Draw("hist");
    resolve_y_range(ymin, ymax, max_in_range(h, xmin, xmax), min_in_range(h, xmin, xmax, false, 1.e-20), plot.logy_, 0.);
    h->GetYaxis()->SetRangeUser(ymin, ymax);
    h->GetXaxis()->SetRangeUser(xmin, xmax);
    if(plot.logy_) c->SetLogy();
    if(plot.logx_) c->SetLogx();
    stamp_.draw(exposure_used_);
    c->RedrawAxis();
    return c;
  }

  //-------------------------------------------------------------------------------------------------------
  // Model with a systematic shift: histograms <hist>_<sys_up> (and <hist>_<sys_down>) in the "sys" type
  TCanvas* plot_systematic(plot_t plot) {
    const TString hist = plot.hist_, type = plot.type_;
    const int selection = plot.selection_, rebin = plot.rebin_, sys_up = plot.sys_up_, sys_down = plot.sys_down_;
    double xmin = plot.xmin_, xmax = plot.xmax_, ymin = plot.ymin_, ymax = plot.ymax_;
    if(debug_ > 0) printf("Plotter::%s: Plotting %s/%s/%i: up = %i, down = %i\n", __func__, hist.Data(), type.Data(), selection, sys_up, sys_down);
    if(sys_up < 0) {
      printf("Plotter::%s: Systematic up index %i is undefined\n", __func__, sys_up);
      return nullptr;
    }

    auto model_of = [&](const TString& h_name, const TString& h_type, const char* tag) {
      auto hists = collect(h_name, h_type, selection, [&](const Process_t& p) { return in_stack(p); });
      TH1* h = sum_hists(hists, Form("%s_%s_%s_%i", tag, hist.Data(), type.Data(), selection));
      for(auto hh : hists) delete hh;
      if(h && rebin > 1) h->Rebin(rebin);
      return h;
    };
    TH1* nominal = model_of(hist, type, "nom");
    TH1* up      = model_of(hist + Form("_%i", sys_up), "sys", "up");
    TH1* down    = (sys_down >= 0) ? model_of(hist + Form("_%i", sys_down), "sys", "down") : nullptr;
    if(!nominal || !up) {
      printf("Plotter::%s: Model not found for %s/%s/%i: up = %i, down = %i\n", __func__, hist.Data(), type.Data(), selection, sys_up, sys_down);
      return nullptr;
    }
    auto sigs = collect(hist, type, selection, [&](const Process_t& p) { return p.is_signal() && !in_stack(p); });
    auto datas = collect(hist, type, selection, [&](const Process_t& p) { return p.is_data(); });
    TH1* data = sum_hists(datas, Form("data_%s_%s_%i", hist.Data(), type.Data(), selection));
    if(rebin > 1) {
      for(auto h : sigs) h->Rebin(rebin);
      if(data) data->Rebin(rebin);
    }

    // Band from the shifted models
    const int nbins = nominal->GetNbinsX();
    std::vector<double> x(nbins), y(nbins), xerr(nbins), eup(nbins), edown(nbins), ry(nbins, 1.), rup(nbins), rdown(nbins);
    double max_r(1.), min_r(1.);
    for(int bin = 1; bin <= nbins; ++bin) {
      const int index = bin - 1;
      const double val = nominal->GetBinContent(bin);
      const double d_up   = up->GetBinContent(bin) - val;
      const double d_down = (down) ? down->GetBinContent(bin) - val : 0.;
      x[index] = nominal->GetBinCenter(bin);
      xerr[index] = nominal->GetBinWidth(bin)/2.;
      y[index] = val;
      eup  [index] = std::max(0., std::max(d_up, d_down));
      edown[index] = std::max(0., -std::min(d_up, d_down));
      if(val <= 0.) { ry[index] = -10.; continue; }
      rup  [index] = eup  [index]/val;
      rdown[index] = edown[index]/val;
      max_r = std::max(max_r, 1. + rup  [index]);
      min_r = std::min(min_r, 1. - rdown[index]);
    }
    auto graph   = new TGraphAsymmErrors(nbins, x.data(), y.data() , xerr.data(), xerr.data(), edown.data(), eup.data());
    auto graph_r = new TGraphAsymmErrors(nbins, x.data(), ry.data(), xerr.data(), xerr.data(), rdown.data(), rup.data());
    for(auto g : {graph, graph_r}) { g->SetFillColor(kRed); g->SetFillStyle(3001); }

    TCanvas* c = new TCanvas(Form("c_sys_%s_%s_%i_%i_%i", hist.Data(), type.Data(), selection, sys_up, ++uid_), "Canvas", 1200, 1000);
    TPad *pad1(nullptr), *pad2(nullptr);
    configure_split_canvas(pad1, pad2, true);
    pad1->cd();
    resolve_x_range(nominal, xmin, xmax);
    style_main_axis(nominal, plot.unit_, plot.xtitle_, plot.ytitle_, 1., true);
    nominal->SetLineWidth(2);
    nominal->SetLineColor(kRed);
    nominal->SetFillStyle(0);
    nominal->Draw("hist");
    graph->Draw("E2");
    double max_val = std::max(max_in_range(nominal, xmin, xmax), max_in_range(up, xmin, xmax));
    if(down) max_val = std::max(max_val, max_in_range(down, xmin, xmax));
    for(auto h : sigs) {
      h->Draw("hist same");
      max_val = std::max(max_val, max_in_range(h, xmin, xmax));
    }
    if(data) {
      data->Draw("EX0 same");
      max_val = std::max(max_val, max_in_range(data, xmin, xmax, true));
    }

    TLegend* leg = new TLegend(pad1->GetLeftMargin() + 0.03, 1. - pad1->GetTopMargin() - 0.25, 1. - pad1->GetRightMargin() - 0.02, 1. - pad1->GetTopMargin() - 0.06);
    leg->SetNColumns(legend_columns_);
    leg->SetLineWidth(0); leg->SetLineColor(0); leg->SetFillColor(0); leg->SetFillStyle(0);
    leg->SetTextSize(legend_text_size_);
    leg->SetTextFont(132);
    if(data) leg->AddEntry(data, data_label_, "PL");
    for(auto h : sigs) leg->AddEntry(h, h->GetTitle(), (h->GetFillStyle() == 0) ? "L" : "F");
    leg->AddEntry(nominal, (signal_mode_ == kStacked) ? "Total" : "Background", "L");
    leg->AddEntry(graph, "Systematic", "F");
    leg->Draw();

    resolve_y_range(ymin, ymax, max_val, min_in_range(nominal, xmin, xmax, false, 1.e-20), plot.logy_, 0.3);
    nominal->GetYaxis()->SetRangeUser(ymin, ymax);
    nominal->GetXaxis()->SetRangeUser(xmin, xmax);
    if(plot.logy_) pad1->SetLogy();
    if(plot.logx_) pad1->SetLogx();

    pad2->cd();
    TH1* haxis_r = (TH1*) nominal->Clone(Form("axis_r_sys_%i", ++uid_));
    haxis_r->SetDirectory(0);
    haxis_r->Reset();
    haxis_r->SetLineWidth(0);
    haxis_r->Draw("hist");
    graph_r->Draw("E2");
    double shift_ymin = min_r - 0.05*(max_r - min_r), shift_ymax = max_r + 0.05*(max_r - min_r);
    if(sys_shift_min_ < sys_shift_max_) {
      shift_ymin = std::max(shift_ymin, sys_shift_min_);
      shift_ymax = std::min(shift_ymax, sys_shift_max_);
      if(shift_ymax <= shift_ymin) { shift_ymin = sys_shift_min_; shift_ymax = sys_shift_max_; }
    }
    haxis_r->GetYaxis()->SetRangeUser(shift_ymin, shift_ymax);
    haxis_r->GetXaxis()->SetRangeUser(xmin, xmax);
    TLine* line = new TLine(xmin, 1., xmax, 1.);
    line->SetLineWidth(2); line->SetLineColor(kBlack); line->SetLineStyle(kDashed);
    line->Draw("same");
    haxis_r->SetYTitle("Shift / Nominal");
    style_ratio_axis(haxis_r, plot.xtitle_, plot.unit_);
    if(plot.logx_) pad2->SetLogx();

    pad1->cd();
    stamp_.draw(exposure_used_);
    pad1->RedrawAxis();
    pad2->RedrawAxis();
    return c;
  }

  //-------------------------------------------------------------------------------------------------------
  // Signal efficiency vs. background rejection for a cut on the variable, for each signal.
  // left: keep events above the cut; eff: draw the efficiency curves instead of the ROC
  TCanvas* plot_roc(plot_t plot, const bool left = true, const bool eff = false) {
    const TString hist = plot.hist_, type = plot.type_;
    const int selection = plot.selection_;
    double xmin = plot.xmin_, xmax = plot.xmax_;
    auto bkgs = collect(hist, type, selection, [&](const Process_t& p) { return p.is_background(); });
    auto sigs = collect(hist, type, selection, [&](const Process_t& p) { return p.is_signal(); });
    if(bkgs.empty() || sigs.empty()) return nullptr;
    TH1* bkg = sum_hists(bkgs, Form("bkg_%s_%s_%i", hist.Data(), type.Data(), selection));
    for(auto h : bkgs) delete h;
    if(plot.rebin_ > 1) {
      bkg->Rebin(plot.rebin_);
      for(auto h : sigs) h->Rebin(plot.rebin_);
    }

    auto to_eff = [&](TH1* h) { // fraction of the total kept by the cut at each bin edge
      const int nbins = h->GetNbinsX();
      const double sum = h->Integral(0, nbins+1);
      for(int bin = 0; bin <= nbins + 1; ++bin) {
        const double above = h->Integral(bin, nbins+1);
        h->SetBinContent(bin, (sum > 0.) ? ((left) ? above/sum : 1. - above/sum) : 0.);
        h->SetBinError(bin, 0.);
      }
    };
    to_eff(bkg);
    std::vector<TGraph*> graphs;
    for(auto h : sigs) {
      to_eff(h);
      auto g = new TGraph();
      g->SetName(Form("%s_roc", h->GetName()));
      for(int bin = 0; bin <= h->GetNbinsX(); ++bin) g->AddPoint(h->GetBinContent(bin), 1. - bkg->GetBinContent(bin));
      g->SetLineColor(h->GetLineColor());
      g->SetLineWidth(std::max(2, (int) h->GetLineWidth()));
      graphs.push_back(g);
    }

    TCanvas* c = new TCanvas(Form("c_roc_%s_%s_%i_%i", hist.Data(), type.Data(), selection, ++uid_), "Canvas", 1000, 750);
    c->SetRightMargin(0.05);
    c->SetTickx(1); c->SetTicky(1);
    TLegend* leg = new TLegend(0.15, 0.15, 0.6, 0.15 + 0.05*(sigs.size() + 1));
    leg->SetLineWidth(0); leg->SetFillStyle(0); leg->SetTextFont(132);
    if(eff) {
      bkg->SetFillStyle(0);
      bkg->SetLineColor(kRed);
      bkg->SetLineWidth(2);
      bkg->SetTitle("");
      bkg->SetXTitle(axis_title(plot.xtitle_, plot.unit_));
      bkg->SetYTitle("Efficiency");
      bkg->Draw("hist");
      bkg->GetYaxis()->SetRangeUser(0., 1.1);
      resolve_x_range(bkg, xmin, xmax);
      bkg->GetXaxis()->SetRangeUser(xmin, xmax);
      leg->AddEntry(bkg, "Background", "L");
      for(auto h : sigs) {
        h->SetFillStyle(0);
        h->Draw("hist same");
        leg->AddEntry(h, h->GetTitle(), "L");
      }
    } else {
      auto haxis = graphs[0];
      haxis->SetTitle("");
      haxis->GetXaxis()->SetTitle("Signal efficiency");
      haxis->GetYaxis()->SetTitle("Background rejection");
      haxis->Draw("AL");
      for(size_t index = 0; index < graphs.size(); ++index) {
        if(index > 0) graphs[index]->Draw("L"); // the first is drawn with the axis
        leg->AddEntry(graphs[index], sigs[index]->GetTitle(), "L");
      }
      haxis->GetXaxis()->SetLimits(0., 1.05);
      haxis->GetYaxis()->SetRangeUser(0., 1.05);
    }
    leg->Draw();
    stamp_.draw(exposure_used_, 0.7);
    c->RedrawAxis();
    return c;
  }

  //-------------------------------------------------------------------------------------------------------
  // Expected yields of each process for [xmin, xmax] (the full histogram including overflows if xmin >= xmax),
  // and S/sqrt(B) and S/sigma(B) for each signal
  int print_yields(TString hist, TString type, int selection, double xmin = 1., double xmax = -1.) {
    printf("Expected yields for %s", layout_.hist_path(hist, type, selection).Data());
    if(xmin < xmax) printf(" in [%g, %g]", xmin, xmax);
    printf(" (N(POT) = %.3g, livetime = %.3g s, N(muons) = %.3g)\n", exposure_used_.npot, exposure_used_.livetime, exposure_used_.nmuons);
    double n_bkg(0.), var_bkg(0.);
    std::vector<std::pair<TString, double>> signals;
    for(auto& p : processes_) {
      if(!p.f_) continue;
      TH1* h = load_hist(p, hist, type, selection);
      TH1* h_raw = load_hist(p, hist, type, selection, nullptr, true);
      if(!h || !h_raw) {
        printf("  %-45s (%-12s): not found\n", display_label(p).Data(), p.name_.Data());
        delete h; delete h_raw;
        continue;
      }
      double err = 0.;
      const double n = integral_in_range(h, xmin, xmax, &err);
      const double n_mc = integral_in_range(h_raw, xmin, xmax);
      printf("  %-45s (%-12s): %11.4g +- %-10.3g (%9.0f MC)\n", display_label(p).Data(), p.name_.Data(), n, err, n_mc);
      if(p.is_background()) { n_bkg += n; var_bkg += err*err; }
      if(p.is_signal()) { // summed by label, as in the plots
        const TString label = display_label(p);
        auto sig = std::find_if(signals.begin(), signals.end(), [&](const std::pair<TString, double>& s) { return s.first == label; });
        if(sig == signals.end()) signals.push_back({label, n});
        else sig->second += n;
      }
      delete h; delete h_raw;
    }
    std::map<TString, TH1*> bkg_sys;
    double sys = 0.;
    if(significance_sys_) {
      bkg_sys = sys_sources(hist, type, selection, 1, [&](const Process_t& p) { return p.is_background(); });
      if(!bkg_sys.empty()) {
        TH1* ref = bkg_sys.begin()->second;
        const int nbins = ref->GetNbinsX();
        int bin_low = 0, bin_high = nbins + 1;
        if(xmin < xmax) {
          bin_low  = (xmin <= ref->GetXaxis()->GetXmin()) ? 0         : ref->GetXaxis()->FindFixBin(xmin + 1.e-6);
          bin_high = (xmax >= ref->GetXaxis()->GetXmax()) ? nbins + 1 : ref->GetXaxis()->FindFixBin(xmax - 1.e-6);
        }
        sys = sys_in_range(bkg_sys, bin_low, bin_high);
      }
      delete_sources(bkg_sys);
    }
    printf("  Total background: %.4g +- %.3g (MC stat) +- %.3g (norm.)\n", n_bkg, std::sqrt(var_bkg), sys);
    for(const auto& sig : signals) {
      printf("  %-45s: S = %.4g, S/sqrt(B) = %.3g", sig.first.Data(), sig.second, (n_bkg > 0.) ? sig.second/std::sqrt(n_bkg) : 0.);
      if(sys > 0.) printf(", S/sigma(B) = %.3g", sig.second/std::sqrt(n_bkg + sys*sys));
      printf("\n");
    }
    return 0;
  }

  //-------------------------------------------------------------------------------------------------------
  // Figure saving
  TString fig_name(const TString& prefix, const plot_t& plot, const TString& suffix = "") const {
    const TString type = (plot.type_ == "") ? TString("") : "_" + plot.type_;
    return Form("%s/%s_%s%s_%i%s%s%s", figdir_.Data(), prefix.Data(), plot.hist_.Data(), type.Data(), plot.selection_,
                suffix.Data(), (plot.normalize_) ? "_norm" : "", (plot.logy_) ? "_log" : "");
  }

  void save(TCanvas* c, const TString& name) const {
    if(!c) return;
    for(const auto& format : formats_) c->SaveAs(name + "." + format);
  }

  TCanvas* print_stack(plot_t plot) {
    TCanvas* c = plot_stack(plot);
    save(c, fig_name("stack", plot));
    return c;
  }

  TCanvas* print_component(plot_t plot, TString process) {
    TCanvas* c = plot_component(plot, process);
    save(c, fig_name("comp_" + process, plot));
    return c;
  }

  TCanvas* print_systematic(plot_t plot) {
    TCanvas* c = plot_systematic(plot);
    save(c, fig_name("sys", plot, Form("_sys_%i", plot.sys_up_)));
    return c;
  }

  TCanvas* print_systematic(plot_t plot, int sys_up, int sys_down = -1) {
    plot.sys_up_ = sys_up; plot.sys_down_ = sys_down;
    return print_systematic(plot);
  }

  TCanvas* print_roc(plot_t plot, const bool left = true, const bool eff = false) {
    TCanvas* c = plot_roc(plot, left, eff);
    save(c, fig_name((eff) ? "eff" : "roc", plot));
    return c;
  }

  // Positional-argument versions, as in the analysis Plotters
  TCanvas* plot_stack(TString hist, TString type, int selection, int rebin = 1, double xmin = 1., double xmax = -1.,
                      double ymin = 1., double ymax = -1., bool logy = false, bool logx = false) {
    return plot_stack(plot_t(hist, type, selection, rebin, xmin, xmax, ymin, ymax, logy, logx));
  }
  TCanvas* print_stack(TString hist, TString type, int selection, int rebin = 1, double xmin = 1., double xmax = -1.,
                       double ymin = 1., double ymax = -1., bool logy = false, bool logx = false) {
    return print_stack(plot_t(hist, type, selection, rebin, xmin, xmax, ymin, ymax, logy, logx));
  }

private:
  int uid_ = 0; // for unique object names
};

}

#endif
