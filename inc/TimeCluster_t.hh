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
#include "EventNtuple/inc/ComboHitInfo.hh"

// local includes
#include "Mu2eEvtAna/inc/GlobalConstants.h"

namespace Mu2eEvtAna {
  struct TimeCluster_t {
    const mu2e::EventNtupleTimeClusterInfo* cluster_;
    // The combo hits of this time cluster, only available if the job stored this collection's hit
    // list (EventNtupleMaker's timeclusters.fillHitsFor); null otherwise
    const std::vector<mu2e::EventNtupleComboHitInfo>* hits_;

    //-------------------------------------------------
    // Accessors

    int   NHits()      const { return (cluster_) ? cluster_->nhits      : -1   ; }
    int   NStrawHits() const { return (cluster_) ? cluster_->nStrawHits : -1   ; }
    float T0()         const { return (cluster_) ? cluster_->t0         :  0.f ; }
    float X()          const { return (cluster_) ? cluster_->pos.x()    :  0.f ; }
    float Y()          const { return (cluster_) ? cluster_->pos.y()    :  0.f ; }
    float Z()          const { return (cluster_) ? cluster_->pos.z()    :  0.f ; }
    float AvgEDep()    const { return (cluster_) ? cluster_->edep       : -1.f ; }
    float ECalo()      const { return (cluster_) ? cluster_->ecalo      : -1.f ; }
    float TCalo()      const { return (cluster_) ? cluster_->tcalo      :  0.f ; }
    bool  HasCalo()    const { return ECalo() >= 0.f; }

    //-------------------------------------------------
    // Combo hit accessors

    bool   HasHits()    const { return hits_ != nullptr; }
    int    NComboHits() const { return (hits_) ? int(hits_->size()) : -1; }

    //-------------------------------------------------
    // Additional functions

    float R() const { return std::sqrt(X()*X() + Y()*Y()); }

    // N(combo hits) of this cluster downstream of ZMin (tracker coordinates, mm). The downstream
    // end of the tracker is the side facing the calorimeter, so this counts the hits of a cluster
    // that reached it -- a handle on charged activity pointing at a calo cluster. Returns -1 if
    // this collection's hit list was not stored, so "no hits stored" is distinguishable from
    // "no hits above ZMin"; check HasHits() before cutting on it.
    int NHitsAboveZ(const float ZMin = 1300.) const {
      if(!hits_) return -1;
      int nhits(0);
      for(const auto& hit : *hits_) if(hit.pos.z() > ZMin) ++nhits;
      return nhits;
    }

    void Reset() {
      cluster_ = nullptr;
      hits_    = nullptr;
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
