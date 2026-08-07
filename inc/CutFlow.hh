//
// Help track a cut-flow
// Michael MacKenzie (2026)

#ifndef MU2EEVTANA_CUTFLOW_HH
#define MU2EEVTANA_CUTFLOW_HH

#include <map>
#include <vector>
#include <string>
#include <iostream>

namespace Mu2eEvtAna {
  class CutFlow {
  public:

    CutFlow() {
      ResetEvent();
      cuts_.reserve(1000);
      indices_.reserve(1000);
      counts_.reserve(1000);
    }

    long Increment(std::string cut, const bool check_if_applied = true) {
      if(!counts_.contains(cut)) {
        cuts_.push_back(cut);
        counts_[cut] = 1;
        indices_[cut] = cuts_.size() - 1;
        event_index_ = indices_[cut];
      } else if(!check_if_applied || event_index_ < indices_[cut]) {
        ++counts_[cut];
        event_index_ = indices_[cut];
      }
      return counts_[cut];
    }

    void ResetEvent() {
      event_index_ = -1;
    }

    const auto& Cuts()   const { return cuts_; }
    const auto& Counts() const { return counts_; }

    void Print(std::ostream& out = std::cout) const {
      if(cuts_.empty()) return;
      const auto norm = counts_.at(cuts_.at(0));
      for(size_t i = 0; i < cuts_.size(); ++i) {
        const auto& name = cuts_.at(i);
        const auto count = counts_.at(name);
        out << "  Cut " << std::setw(2) << i << ": "
            << std::setw(20) << name
            << " " << std::setw(10) << count
            << " (" << std::defaultfloat << std::setprecision(4) << std::setw(10)
            << (100.*count)/norm << "%)\n";
      }
    }

  private:
    std::vector<std::string> cuts_;
    std::unordered_map<std::string, int> indices_;
    std::unordered_map<std::string, long> counts_;
    int event_index_;
  };
}
#endif
