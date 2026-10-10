#ifndef __MU2EEVTANA_PLOTTER_PLOTUTILS__
#define __MU2EEVTANA_PLOTTER_PLOTUTILS__
// Histogram range, text formatting, and stamp helpers for the Plotter

#include "PlotTypes.C"
#include "physics/Mu2ePhysics.C"

namespace mu2eplot {

  //--------------------------------------------------------------------------------------------------
  // Bin range covering [xmin, xmax], or the full axis if xmin >= xmax
  inline void bin_range(const TH1* h, const double xmin, const double xmax, int& bin_low, int& bin_high) {
    const int nbins = h->GetNbinsX();
    bin_low  = (xmin < xmax) ? std::max(1    , h->GetXaxis()->FindFixBin(xmin + 1.e-6)) : 1;
    bin_high = (xmin < xmax) ? std::min(nbins, h->GetXaxis()->FindFixBin(xmax - 1.e-6)) : nbins;
  }

  //--------------------------------------------------------------------------------------------------
  inline double min_in_range(TH1* h, double xmin, double xmax, bool use_errors = false, double cut_off = 0.) {
    if(!h) return cut_off;
    int bin_low, bin_high; bin_range(h, xmin, xmax, bin_low, bin_high);
    double min_val = std::max(h->GetMaximum(), cut_off);
    for(int ibin = bin_low; ibin <= bin_high; ++ibin) {
      double binc = h->GetBinContent(ibin);
      if(binc <= cut_off) continue;
      if(use_errors) binc -= h->GetBinError(ibin);
      if(binc <= cut_off) continue;
      min_val = std::min(min_val, binc);
    }
    return min_val;
  }

  //--------------------------------------------------------------------------------------------------
  inline double max_in_range(TH1* h, double xmin, double xmax, bool use_errors = false) {
    if(!h) return 0.;
    int bin_low, bin_high; bin_range(h, xmin, xmax, bin_low, bin_high);
    double max_val = 0.;
    for(int ibin = bin_low; ibin <= bin_high; ++ibin) {
      double binc = h->GetBinContent(ibin);
      if(use_errors) binc += h->GetBinError(ibin);
      max_val = std::max(max_val, binc);
    }
    return max_val;
  }

  //--------------------------------------------------------------------------------------------------
  // Integral over [xmin, xmax], including the underflow/overflow if the range reaches past the axis
  inline double integral_in_range(TH1* h, double xmin, double xmax, double* error = nullptr) {
    if(!h) return 0.;
    const int nbins = h->GetNbinsX();
    int bin_low  = 0, bin_high = nbins + 1;
    if(xmin < xmax) {
      bin_low  = (xmin <= h->GetXaxis()->GetXmin()) ? 0         : h->GetXaxis()->FindFixBin(xmin + 1.e-6);
      bin_high = (xmax >= h->GetXaxis()->GetXmax()) ? nbins + 1 : h->GetXaxis()->FindFixBin(xmax - 1.e-6);
    }
    double err = 0.;
    const double val = h->IntegralAndError(bin_low, bin_high, err);
    if(error) *error = err;
    return val;
  }

  //--------------------------------------------------------------------------------------------------
  // "1.5 #times 10^{-13}", or "10^{-13}" for a unit mantissa
  inline TString format_sci(const double value, const int digits = 1) {
    if(value == 0.) return "0";
    const double log_val = std::log10(std::fabs(value)) + 1.e-9;
    const int power = (int) std::floor(log_val);
    const double mantissa = value/std::pow(10., power);
    if(power == 0) return Form("%.*f", digits, mantissa);
    if(std::fabs(mantissa - 1.) < 0.5*std::pow(10., -digits)) return Form("10^{%i}", power);
    return Form("%.*f #times 10^{%i}", digits, mantissa, power);
  }

  //--------------------------------------------------------------------------------------------------
  // Text drawn above the top pad
  struct Stamp_t {
    bool    draw_         = true        ;
    TString experiment_   = "Mu2e"      ;
    TString status_       = "Simulation";
    TString text_         = ""          ; // if set, replaces the exposure text
    TString extra_        = ""          ; // appended to the exposure text, e.g. "R_{#mue} = 10^{-13}"
    bool    show_power_   = true        ;
    bool    show_pot_     = false       ;
    bool    show_running_ = true        ; // running time (livetime / duty cycle) in days, or years above a year
    bool    show_livetime_= false       ; // on-spill livetime in seconds
    bool    show_muons_   = true        ;

    TString exposure_text(const Exposure_t& exp) const {
      if(text_ != "") return text_;
      std::vector<TString> parts;
      double power = exp.beam_power_kw;
      if(power <= 0.) power = mu2e_physics::beam_power_kw(exp.npot, exp.livetime, exp.duty_cycle);
      if(show_power_    && power        > 0.) parts.push_back(Form("%.3g kW beam", power));
      if(show_pot_      && exp.npot     > 0.) parts.push_back(format_sci(exp.npot) + " POT");
      if(show_running_  && exp.livetime > 0. && exp.duty_cycle > 0.) {
        const double days = exp.livetime/exp.duty_cycle/(24.*60.*60.);
        parts.push_back((days > 365.) ? Form("%.3g years running", days/365.) : Form("%.3g days running", days));
      }
      if(show_livetime_ && exp.livetime > 0.) parts.push_back(format_sci(exp.livetime) + " s On-Spill");
      if(show_muons_    && exp.nmuons   > 0.) parts.push_back(format_sci(exp.nmuons) + " muon stops");
      if(extra_ != "") parts.push_back(extra_);
      TString text;
      for(size_t i = 0; i < parts.size(); ++i) text += (i > 0) ? "; " + parts[i] : parts[i];
      return text;
    }

    // Draw on the current pad
    TLatex* draw(const Exposure_t& exp, const double scale = 1.) const {
      if(!draw_ || !gPad) return nullptr;
      TLatex* logo = new TLatex();
      logo->SetNDC();
      const float text_size = 0.042*1.25*scale;
      const float x0(gPad->GetLeftMargin() + 0.015), y0(1. - (gPad->GetTopMargin() - 0.017));
      logo->SetTextAlign(11);
      logo->SetTextSize(text_size);
      logo->SetTextFont(61);
      logo->DrawLatex(x0, y0, experiment_);
      logo->SetTextSize(0.042*scale);
      logo->SetTextFont(52);
      const double r = (gPad->GetWh()*gPad->GetAbsHNDC())/(gPad->GetWw()*gPad->GetAbsWNDC()); // pad height / width in pixels
      logo->DrawLatex(x0 + 0.72*experiment_.Length()*text_size*r, y0, status_); // after the experiment name
      const TString text = exposure_text(exp);
      if(text != "") {
        logo->SetTextSize(0.76*text_size);
        logo->SetTextFont(132);
        logo->SetTextAlign(31);
        logo->DrawLatex(1. - gPad->GetRightMargin(), y0, text);
      }
      return logo;
    }
  };

  //--------------------------------------------------------------------------------------------------
  // Delete a canvas and everything drawn in it. The canvas' own frame is left to the canvas, which deletes it on
  // close (deleting it here too crashes for canvases drawn without sub-pads).
  inline int empty_canvas(TCanvas* c) {
    if(!c) return 0;
    TList* list = c->GetListOfPrimitives();
    if(list) {
      while(TObject* frame = list->FindObject("TFrame")) list->Remove(frame);
      list->Delete();
    }
    delete c;
    return 0;
  }

  //--------------------------------------------------------------------------------------------------
  // Count a missing canvas as a failure, and clean up an existing one
  inline bool handle_canvas(TCanvas* c, int& status, const bool cleanup = true) {
    if(!c) {
      ++status;
      return false;
    }
    if(cleanup) empty_canvas(c);
    return true;
  }
}

#endif
