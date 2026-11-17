//
// Help track N-1 efficiencies of a set of CutID bits
// Michael MacKenzie (2026)

#ifndef MU2EEVTANA_NMINUSONE_HH
#define MU2EEVTANA_NMINUSONE_HH

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include "Mu2eEvtAna/inc/CutID.hh"

namespace Mu2eEvtAna {
  // For each registered ID bit, count candidates that pass all registered cuts except that one.
  // Counts are per event (at most once each), unweighted. Only registered bits are considered;
  // any other bits set in the ID are ignored.
  class NMinusOne {
  public:

    NMinusOne() { ResetEvent(); }

    // Register a cut bit; registration order is the print order
    void AddCut(int bit, const std::string& name) {
      bits_.push_back(bit);
      names_.push_back(name);
      n_minus_.push_back(0);
      counted_minus_.push_back(false);
      mask_ |= (1 << bit);
    }

    void ResetEvent() {
      counted_presel_ = false;
      counted_all_ = false;
      for(size_t i = 0; i < counted_minus_.size(); ++i) counted_minus_[i] = false;
    }

    // Call for a candidate that passed the non-ID preselection
    void Fill(const CutID& ID) {
      if(!counted_presel_) { ++n_presel_; counted_presel_ = true; }
      if(!counted_all_ && ID.Passes(mask_)) { ++n_all_; counted_all_ = true; }
      for(size_t i = 0; i < bits_.size(); ++i) {
        if(counted_minus_[i]) continue;
        if(ID.Passes(mask_ & ~(1 << bits_[i]))) { ++n_minus_[i]; counted_minus_[i] = true; }
      }
    }

    long NPresel() const { return n_presel_; }
    long NAll()    const { return n_all_; }
    long NMinus(size_t i) const { return n_minus_.at(i); }
    const auto& Names() const { return names_; }
    const auto& Bits()  const { return bits_; }

    void Print(std::ostream& out = std::cout) const {
      if(n_presel_ == 0) {
        out << "  N-1 (no candidates passed the preselection)\n";
        return;
      }
      const auto flags = out.flags();
      const auto prec = out.precision();
      out << "  N-1 summary: N(presel) = " << n_presel_
          << ", N(all cuts) = " << n_all_
          << ", acceptance = " << std::fixed << std::setprecision(3)
          << (100.*n_all_)/n_presel_ << "%\n";
      out << "  N-1 " << std::left << std::setw(22) << "removed cut" << std::right
          << std::setw(12) << "N(all but)"
          << std::setw(12) << "N-1 eff"
          << std::setw(14) << "acc w/o cut" << "\n";
      for(size_t i = 0; i < bits_.size(); ++i) {
        const long n = n_minus_[i];
        out << "  N-1 " << std::left << std::setw(22) << names_[i] << std::right
            << std::setw(12) << n;
        if(n > 0) out << std::setw(11) << std::fixed << std::setprecision(3) << (100.*n_all_)/n << "%";
        else      out << std::setw(12) << "-";
        out << std::setw(13) << std::fixed << std::setprecision(3) << (100.*n)/n_presel_ << "%\n";
      }
      out.flags(flags);
      out.precision(prec);
    }

  private:
    std::vector<int> bits_;
    std::vector<std::string> names_;
    std::vector<long> n_minus_;
    std::vector<bool> counted_minus_;
    int mask_ = 0;
    long n_presel_ = 0;
    long n_all_ = 0;
    bool counted_presel_ = false;
    bool counted_all_ = false;
  };
}
#endif
