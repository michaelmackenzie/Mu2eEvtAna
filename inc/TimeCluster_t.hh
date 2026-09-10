//
// Time cluster information
// Michael MacKenzie (2026)

#ifndef MU2EEVTANA_TIMECLUSTER_T_HH
#define MU2EEVTANA_TIMECLUSTER_T_HH

// ROOT includes
#include "Rtypes.h"
#include "TString.h"

// EventNtuple includes
#include "EventNtuple/inc/TimeClusterInfo.hh"

namespace Mu2eEvtAna {
  struct TimeCluster_t {
    const mu2e::EventNtupleTimeClusterInfo* cluster_;

    //-------------------------------------------------
    // Accessors

    int   NHits()      const { return (cluster_) ? cluster_->nhits      : -1   ; }
    int   NStrawHits() const { return (cluster_) ? cluster_->nStrawHits : -1   ; }
    float T0()         const { return (cluster_) ? cluster_->t0         :  0.f ; }
    float X()          const { return (cluster_) ? cluster_->pos.x()    :  0.f ; }
    float Y()          const { return (cluster_) ? cluster_->pos.y()    :  0.f ; }
    float Z()          const { return (cluster_) ? cluster_->pos.z()    :  0.f ; }
    float ECalo()      const { return (cluster_) ? cluster_->ecalo      : -1.f ; }
    float TCalo()      const { return (cluster_) ? cluster_->tcalo      :  0.f ; }
    bool  HasCalo()    const { return ECalo() >= 0.f; }

    //-------------------------------------------------
    // Additional functions

    float R() const { return std::sqrt(X()*X() + Y()*Y()); }

    void Reset() {
      cluster_ = nullptr;
    }

    TimeCluster_t() { Reset(); }

    //----------------------------------------------
    // Print the time cluster
    void Print(TString opt = "") const {
      opt.ToLower();
      if(opt.Contains("banner")) {
        std::string filler(90, '-');
        printf("%s\n", filler.c_str());
        printf("Idx: %5s %5s %10s %10s %10s %10s %10s %10s\n", "NHit", "NSt", "T0", "x", "y", "z", "E(calo)", "T(calo)");
        printf("%s\n", filler.c_str());
      }
      if(!cluster_) return;
      printf("Idx: %5i %5i %10.1f %10.2f %10.2f %10.2f %10.2f %10.1f\n",
             NHits(), NStrawHits(), T0(), X(), Y(), Z(), ECalo(), TCalo());
    }
  };
}
#endif
