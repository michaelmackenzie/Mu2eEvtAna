//
// SpectrumWeight: physics weight for a sample generated flat in momentum, from a tabulated spectrum
// Michael MacKenzie (2026)

#ifndef MU2EEVTANA_SPECTRUMWEIGHT_HH
#define MU2EEVTANA_SPECTRUMWEIGHT_HH

// standard includes
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ROOT includes
#include "TString.h"
#include "TSystem.h"

namespace Mu2eEvtAna {
  // Reads an Offline spectrum table (two columns: total energy [MeV], relative rate, as used by the
  // StoppedParticleReactionGun "tabulated" spectra) and normalizes it to a probability density.
  //
  // Weight(p) is the density in momentum, dN/dp = f(E) dE/dp = f(E) p/E. For a sample generated flat in
  // momentum over [pmin, pmax] and normalized as rate * (pmax - pmin) / N(gen) -- the usual flat-sample
  // normalization -- weighting each event by Weight(p) gives the physical spectrum.
  class SpectrumWeight {
  public:
    // file: path relative to MU2E_SEARCH_PATH (e.g. "Offline/EventGenerator/data/czarnecki_szafron_Al_2016.tbl"), or absolute
    SpectrumWeight(const TString& file, const double mass) : mass_(mass) {
      const TString path = Resolve(file);
      std::ifstream in(path.Data());
      if(!in.is_open()) throw std::runtime_error(Form("SpectrumWeight: Unable to open spectrum table %s", file.Data()));
      std::string line;
      while(std::getline(in, line)) {
        std::istringstream ss(line);
        double e, val;
        if(line.empty() || line[0] == '#' || !(ss >> e >> val)) continue;
        energy_.push_back(e);
        value_.push_back(val);
      }
      if(energy_.size() < 2) throw std::runtime_error(Form("SpectrumWeight: No spectrum read from %s", path.Data()));
      // Tables are given in either energy order: sort ascending
      std::vector<size_t> order(energy_.size());
      for(size_t i = 0; i < order.size(); ++i) order[i] = i;
      std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return energy_[a] < energy_[b]; });
      std::vector<double> energy, value;
      for(auto i : order) { energy.push_back(energy_[i]); value.push_back(value_[i]); }
      energy_ = energy;
      value_  = value;
      // Normalize to unit area (trapezoidal integration over the table)
      double norm = 0.;
      for(size_t i = 1; i < energy_.size(); ++i) norm += 0.5*(value_[i] + value_[i-1])*(energy_[i] - energy_[i-1]);
      if(norm <= 0.) throw std::runtime_error(Form("SpectrumWeight: Non-positive spectrum integral in %s", path.Data()));
      for(auto& val : value_) val /= norm;
      file_ = path;
    }

    // Probability density in total energy (per MeV), linearly interpolated
    double DensityE(const double e) const {
      if(e <= energy_.front() || e >= energy_.back()) return 0.;
      size_t lo = 0, hi = energy_.size() - 1;
      while(hi - lo > 1) {
        const size_t mid = (lo + hi)/2;
        if(energy_[mid] <= e) lo = mid; else hi = mid;
      }
      const double frac = (e - energy_[lo])/(energy_[hi] - energy_[lo]);
      return value_[lo] + frac*(value_[hi] - value_[lo]);
    }

    // Probability density in momentum (per MeV/c)
    double Weight(const double p) const {
      if(p <= 0.) return 0.;
      const double e = std::sqrt(p*p + mass_*mass_);
      return DensityE(e)*p/e;
    }

    // Fraction of the spectrum with momentum in [pmin, pmax]
    double Fraction(const double pmin, const double pmax, const int nsteps = 20000) const {
      double sum = 0.;
      const double dp = (pmax - pmin)/nsteps;
      for(int i = 0; i < nsteps; ++i) sum += Weight(pmin + (i + 0.5)*dp)*dp;
      return sum;
    }

    const TString& File() const { return file_; }

  private:
    static TString Resolve(const TString& file) {
      if(file.BeginsWith("/")) return file;
      const char* search = gSystem->Getenv("MU2E_SEARCH_PATH");
      if(search) {
        std::stringstream ss(search);
        std::string dir;
        while(std::getline(ss, dir, ':')) {
          const TString path = TString(dir.c_str()) + "/" + file;
          if(!gSystem->AccessPathName(path)) return path;
        }
      }
      return file;
    }

    double mass_;
    TString file_;
    std::vector<double> energy_;
    std::vector<double> value_;
  };
}

#endif
