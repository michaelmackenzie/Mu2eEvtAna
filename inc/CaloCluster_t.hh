//
// Calo cluster information
// Michael MacKenzie (2025)

#ifndef MU2EEVTANA_CALOLUSTER_T_HH
#define MU2EEVTANA_CALOLUSTER_T_HH

// ROOT includes
#include "Rtypes.h"

// EventNtuple includes
#include "EventNtuple/inc/CaloClusterInfo.hh"
#include "EventNtuple/inc/CaloHitInfo.hh"

namespace Mu2eEvtAna {
  // Forward declarations only: CaloCluster_t just holds pointers to these (set by whichever
  // analysis module does the matching, e.g. Run1BAna), so it doesn't need their full definitions.
  struct Track_t;
  struct LineSeed_t;
  struct TimeCluster_t;
  struct CRVCluster_t;

  struct CaloCluster_t {
    const mu2e::CaloClusterInfo* cluster_;
    const mu2e::CaloClusterInfoMC* cluster_mc_;

    // Best-matched objects for this cluster (nullptr if none/not searched for). Populated by an
    // analysis module, not by the base Mu2eEvtAna -- e.g. Run1BAna::MatchCaloClusters().
    Track_t*       line_        ; // best matched track (line fit)
    LineSeed_t*    line_seed_   ; // best matched line seed
    TimeCluster_t* time_cluster_; // best matched time cluster
    CRVCluster_t*  crv_cluster_ ; // best matched CRV coincidence cluster


    //-------------------------------------------------
    // Accessors

    float Energy   () const { return (cluster_) ? cluster_->energyDep_    : -1.f ; }
    float EnergyErr() const { return (cluster_) ? cluster_->energyDepErr_ : -1.f ; }
    float Time     () const { return (cluster_) ? cluster_->time_         : -1.f ; }
    float TimeErr  () const { return (cluster_) ? cluster_->timeErr_      : -1.f ; }
    float X        () const { return (cluster_) ? cluster_->cog_.x()      : 0.f  ; }
    float Y        () const { return (cluster_) ? cluster_->cog_.y()      : 0.f  ; }
    float Z        () const { return (cluster_) ? cluster_->cog_.z()      : 0.f  ; }
    int   DiskID   () const { return (cluster_) ? cluster_->diskID_       : -1   ; }
    int   NCrystals() const { return (cluster_) ? cluster_->size_         : -1   ; }
    int   IsSplit  () const { return (cluster_) ? cluster_->isSplit_      : -1   ; }


    //-------------------------------------------------
    // Additional functions

    float R        () const { return std::sqrt(X()*X() + Y()*Y()); }


    void Reset() {
      cluster_ = nullptr;
      cluster_mc_ = nullptr;
      line_ = nullptr;
      line_seed_ = nullptr;
      time_cluster_ = nullptr;
      crv_cluster_ = nullptr;
    }

    CaloCluster_t() { Reset(); }

    //----------------------------------------------
    // Print the cluster
    void Print(TString opt = "") const {
      opt.ToLower();
      if(opt.Contains("banner")) {
        std::string filler(130, '-');
        printf("%s\n", filler.c_str());
        printf("Idx: %10s %10s %10s %10s %4s %10s %10s %10s %10s\n",
               "Energy", "E err", "Time", "T err", "Disk", "N(crystals)", "x", "y", "z");
        printf("%s\n", filler.c_str());
      }
      if(!cluster_) return;
      printf("Idx: %10.1f %10.3f %10.1f %10.3f %4i %10i %10.2f %10.2f %10.2f\n",
             Energy(), EnergyErr(), Time(), TimeErr(), DiskID(), NCrystals(), X(), Y(), Z());
    }
  };
}
#endif
