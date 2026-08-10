//
// Help track selection cuts
// Michael MacKenzie (2026)

#ifndef MU2EEVTANA_CUTID_HH
#define MU2EEVTANA_CUTID_HH

#include <map>
#include <vector>
#include <string>
#include <iostream>

namespace Mu2eEvtAna {
  class CutID {
  public:

    CutID(int id = 0) : id_(id) {
    }

    void SetBit(const int bit) {
      if(bit >= maxBit_) throw std::runtime_error("Exceeded maximum selection bit!");
      id_ = id_ | (1 << bit);
    }

    // No bits are set
    bool Passes(int mask = 0xffff) const { return (id_ & mask) == 0; }

    // Check if a specific bit is set
    bool CheckBit(const int bit) const {
      if(bit >= maxBit_) throw std::runtime_error("Exceeded maximum selection bit!");
      return ((id_ & (1 << bit)) != 0);
    }

    // Check if only this bit is set
    bool CheckBitExclusive(const int bit) const {
      if(bit >= maxBit_) throw std::runtime_error("Exceeded maximum selection bit!");
      return CheckBit(bit) && ((id_ & ~(1 << bit)) == 0);
    }

    int ID(const int mask = 0xffff) const { return id_ & mask; }

    void Reset() { id_ = 0; }

    void Print(std::ostream& out = std::cout) const {
      out << "ID = " << std::hex << id_ << std::dec;
    }

    friend bool operator==(const CutID& obj, const int val) {
      return obj.ID() == val;
    }

    friend bool operator==(const int val, const CutID& obj) {
      return val == obj.ID();
    }

    friend bool operator==(const CutID& a, CutID& b) {
      return a.ID() == b.ID();
    }

    operator int() const {
        return id_;
    }

  private:
    static constexpr int maxBit_ = 32;
    int id_;
  };
}
#endif
