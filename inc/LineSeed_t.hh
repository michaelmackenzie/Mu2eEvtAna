//
// Line seed information (cosmic/straight-line track seed)
// Michael MacKenzie (2026)

#ifndef MU2EEVTANA_LINESEED_T_HH
#define MU2EEVTANA_LINESEED_T_HH

// ROOT includes
#include "Rtypes.h"
#include "TString.h"

// EventNtuple includes
#include "EventNtuple/inc/LineSeedInfo.hh"

namespace Mu2eEvtAna {
  struct LineSeed_t {
    const mu2e::LineSeedInfo* seed_;

    //-------------------------------------------------
    // Accessors

    int   Status()     const { return (seed_) ? seed_->status     :  0    ; }
    int   NHits()      const { return (seed_) ? seed_->nhits      : -1    ; }
    int   NStrawHits() const { return (seed_) ? seed_->nStrawHits : -1    ; }
    float T0()         const { return (seed_) ? seed_->t0         :  0.f  ; }
    float D0()         const { return (seed_) ? seed_->d0         : -1.e6 ; }
    float Phi0()       const { return (seed_) ? seed_->phi0       : -1.e6 ; }
    float Z0()         const { return (seed_) ? seed_->z0         : -1.e6 ; }
    float Cos()        const { return (seed_) ? seed_->cos        : -1.e6 ; }
    float A0()         const { return (seed_) ? seed_->A0         : -1.e6 ; }
    float B0()         const { return (seed_) ? seed_->B0         : -1.e6 ; }
    float A1()         const { return (seed_) ? seed_->A1         : -1.e6 ; }
    float B1()         const { return (seed_) ? seed_->B1         : -1.e6 ; }
    float ECalo()      const { return (seed_) ? seed_->ecalo      : -1.f  ; }
    float TCalo()      const { return (seed_) ? seed_->tcalo      :  0.f  ; }
    bool  HasCalo()    const { return ECalo() >= 0.f; }

    //-------------------------------------------------
    // Additional functions

    bool IsGood() const { return Status() > 0; }

    void Reset() {
      seed_ = nullptr;
    }

    LineSeed_t() { Reset(); }

    //----------------------------------------------
    // Print the line seed
    void Print(TString opt = "") const {
      opt.ToLower();
      if(opt.Contains("banner")) {
        std::string filler(110, '-');
        printf("%s\n", filler.c_str());
        printf("Idx: %6s %5s %5s %10s %10s %10s %10s %10s\n", "Status", "NHit", "NSt", "T0", "d0", "phi0", "z0", "cos");
        printf("%s\n", filler.c_str());
      }
      if(!seed_) return;
      printf("Idx: %6i %5i %5i %10.1f %10.2f %10.3f %10.2f %10.3f\n",
             Status(), NHits(), NStrawHits(), T0(), D0(), Phi0(), Z0(), Cos());
    }
  };
}
#endif
