#ifndef __MU2EEVTANA_PLOTTER_PLOTTYPES__
#define __MU2EEVTANA_PLOTTER_PLOTTYPES__
// Types describing what to plot. None of these carry physics: an analysis fills them from its own
// dataset and rate definitions (see examples/), and the Plotter only multiplies what it is given.

#include <map>
#include <vector>

namespace mu2eplot {

  //--------------------------------------------------------------------------------------------------
  // Enums

  // Process roles. The values match the type codes of the analysis Plotters this replaces.
  enum Role_t { kAnyRole = -10, kSignal = -1, kData = 0, kBackground = 1 };

  // Which exposure a process' normalization is per unit of
  enum ExposureKind_t {
    kPOT,       // beam processes: expected events per MC event per POT
    kLivetime,  // cosmic rays: per second of livetime
    kNEvents,   // pileup-only samples: per event (microbunch)
    kMuonStops, // per stopped muon
    kFixed      // already normalized (e.g. data)
  };

  // Lower pad content
  enum LowerPad_t {
    kNoPad,        // single pad
    kRatio,        // data / model, or signal / background without data
    kDifference,   // data - model, or signal - background without data
    kSignificance  // S/sqrt(B) per signal, or S/sigma(B) including the process normalization uncertainties
  };

  // How signals are drawn with the background stack
  enum SignalMode_t {
    kOverlay, // signals drawn as lines on top of the background stack
    kStacked  // signals added on top of the background stack, and also drawn alone in front of it
  };

  // Significance range, for kSignificance lower pads
  enum SignificanceRange_t {
    kPerBin, // S and B in each bin
    kAbove,  // S and B integrated from the bin to the upper edge (a lower cut at the bin's low edge)
    kBelow   // S and B integrated from the lower edge to the bin (an upper cut at the bin's high edge)
  };

  //--------------------------------------------------------------------------------------------------
  // Integrated exposure the expected yields are normalized to
  struct Exposure_t {
    double npot          = -1.; // N(protons on target)
    double livetime      = -1.; // on-spill livetime (s)
    double nmuons        = -1.; // N(stopped muons)
    double nevents       = -1.; // N(microbunch events)
    double duty_cycle    = -1.; // for the beam power in the stamp, if beam_power_kw is not given
    double beam_power_kw = -1.; // beam power for the stamp; computed from npot/livetime/duty_cycle if <= 0
    TString label        = "" ; // free-form name, e.g. "Run 1A"

    double value(const ExposureKind_t kind) const {
      switch(kind) {
      case kPOT      : return npot    ;
      case kLivetime : return livetime;
      case kNEvents  : return nevents ;
      case kMuonStops: return nmuons  ;
      default        : return 1.      ;
      }
    }

    void scale(const double factor) {
      if(npot     > 0.) npot     *= factor;
      if(livetime > 0.) livetime *= factor;
      if(nmuons   > 0.) nmuons   *= factor;
      if(nevents  > 0.) nevents  *= factor;
    }
  };

  //--------------------------------------------------------------------------------------------------
  // Expected events per MC event for a sample generated with a given rate
  //   rate      : expected rate of the process per unit exposure (e.g. per POT), within the generated range
  //   ngen      : N(generated events) in the sample (or the generated livetime, for cosmic rays)
  //   emin/emax : flat generation range, if the sample was generated flat and the rate is a density in it
  inline double mc_norm(const double rate, const double ngen, const double emin = 0., const double emax = 1.) {
    if(ngen <= 0.) {
      printf("mu2eplot::mc_norm: Non-positive N(generated) = %g, returning 0\n", ngen);
      return 0.;
    }
    double norm = rate/ngen;
    if(emin < emax) norm *= (emax - emin);
    return norm;
  }

  //--------------------------------------------------------------------------------------------------
  // A process (one histogram file) to include in the plots
  struct Process_t {
    TString        name_    = ""         ; // unique key
    TString        label_   = ""         ; // legend label; processes with the same label are summed
    Role_t         role_    = kBackground;
    int            color_   = kRed       ;
    TString        file_    = ""         ; // histogram file path

    // Normalization: histogram * norm_ * exposure(exposure_) * scale_ * (N(expected) / N(seen))
    double         norm_     = 1.  ; // expected events per MC event per unit exposure (see mc_norm)
    ExposureKind_t exposure_ = kPOT; // which exposure norm_ is per unit of
    double         scale_    = 1.  ; // additional factor, e.g. a signal branching fraction
    TString        scale_symbol_ = ""; // if set and scale_ != 1, the label gains " [<symbol> = <scale_>]"
    Long64_t       n_expected_   = 0 ; // N(events) expected in the input sample; if > 0, corrects for unprocessed events

    // Control-region shapes: read the histogram from set (selection + set_offset_) instead
    int            set_offset_        = 0    ;
    double         offset_scale_      = 1.   ; // applied only when the offset set is used
    bool           offset_use_nominal_norm_ = false; // keep the nominal set's yield, take only the shape from the offset set

    // Fractional normalization uncertainties by source. Contributions of the same source are summed linearly
    // (fully correlated) and different sources are added in quadrature
    std::map<TString, double> sys_;

    // Style overrides (< 0 uses the Plotter's defaults for the role)
    int            line_color_ = -1;
    int            line_style_ = -1;
    int            line_width_ = -1;
    int            fill_style_ = -1;

    // Set by Plotter::init
    TFile*         f_           = nullptr;
    double         sample_corr_ = 1.     ; // N(expected) / N(seen)
    double         n_seen_      = 0.     ;

    Process_t() {}
    Process_t(TString name, TString label, Role_t role, int color, TString file,
              double norm = 1., ExposureKind_t exposure = kPOT, double scale = 1.) :
      name_(name), label_(label), role_(role), color_(color), file_(file),
      norm_(norm), exposure_(exposure), scale_(scale) {}

    // Chainable setters
    Process_t& expected(Long64_t n)                    { n_expected_ = n; return *this; }
    Process_t& rate_scale(double s, TString sym = "")  { scale_ = s; if(sym != "") scale_symbol_ = sym; return *this; }
    Process_t& offset(int off, double s = 1., bool use_nominal_norm = false) {
      set_offset_ = off; offset_scale_ = s; offset_use_nominal_norm_ = use_nominal_norm; return *this;
    }
    Process_t& sys(TString source, double frac)       { sys_[source] = frac; return *this; }
    Process_t& style(int line_style, int line_width = -1, int fill_style = -1, int line_color = -1) {
      line_style_ = line_style; line_width_ = line_width; fill_style_ = fill_style; line_color_ = line_color; return *this;
    }

    bool is_signal    () const { return role_ == kSignal    ; }
    bool is_background() const { return role_ == kBackground; }
    bool is_data      () const { return role_ == kData      ; }
  };

  //--------------------------------------------------------------------------------------------------
  // Where histograms and normalization information live in the histogram files
  struct Layout_t {
    // Histogram path, with {type}, {set}, and {hist} replaced
    TString hist_format_ = "Ana/Hist/{type}_{set}/{hist}";
    // Normalization TTree, whose norm_branch_ is summed over entries to give N(events seen)
    TString norm_path_   = "Ana/data/Norm";
    TString norm_branch_ = "nseen";

    TString hist_path(const TString& hist, const TString& type, const int set) const {
      TString path(hist_format_);
      path.ReplaceAll("{type}", type);
      path.ReplaceAll("{set}" , Form("%i", set));
      path.ReplaceAll("{hist}", hist);
      return path;
    }

    // Mu2eEvtAna analyzer outputs (ConvAna, RMCAna, Run1BAna, ...)
    static Layout_t evtana() { return Layout_t(); }
  };

  //--------------------------------------------------------------------------------------------------
  // A plot
  struct plot_t {
    TString hist_      = ""   ;
    TString type_      = ""   ;
    int     selection_ = 0    ;
    int     rebin_     = 1    ;
    double  xmin_      =  1.  ; // xmin >= xmax uses the full axis
    double  xmax_      = -1.  ;
    double  ymin_      =  1.  ; // ymin >= ymax chooses the range automatically
    double  ymax_      = -1.  ;
    bool    logy_      = false;
    bool    logx_      = false;
    TString xtitle_    = ""   ;
    TString unit_      = ""   ;
    TString ytitle_    = ""   ;
    TString title_     = ""   ;
    bool    normalize_ = false; // shape comparison: each signal, the background total, and the data scaled to unit area
    int     sys_up_    = -1   ;
    int     sys_down_  = -1   ;
    plot_t(TString hist = "", TString type = "", int selection = 0,
           int rebin = 1, double xmin = 1., double xmax = -1., double ymin = 1., double ymax = -1.,
           bool logy = false, bool logx = false,
           TString xtitle = "", TString unit = "", TString ytitle = "", TString title = "") :
      hist_(hist), type_(type), selection_(selection), rebin_(rebin), xmin_(xmin), xmax_(xmax),
      ymin_(ymin), ymax_(ymax), logy_(logy), logx_(logx), xtitle_(xtitle), unit_(unit), ytitle_(ytitle), title_(title)
    {}
    plot_t& normalized(bool val = true) { normalize_ = val; return *this; }
  };
}

#endif
