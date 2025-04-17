//
// A signle Track
// Michael MacKenzie (2024)

#ifndef MU2EEVTANA_TRACK_T_HH
#define MU2EEVTANA_TRACK_T_HH

// ROOT includes
#include "Rtypes.h"

// Event Ntuple includes
#include "EventNtuple/inc/TrkInfo.hh"
#include "EventNtuple/inc/TrkCaloHitInfo.hh"
#include "EventNtuple/inc/TrkSegInfo.hh"
#include "EventNtuple/utils/rooutil/inc/Track.hh"

// local includes
#include "Mu2eEvtAna/inc/GlobalConstants.h"

namespace Mu2eEvtAna {
  struct Track_t {
    Track*   track_             = nullptr; // pointer to the EventNtuple::Track object
    mu2e::TrkInfo* info_        = nullptr; // info about track fit
    mu2e::TrkCaloHitInfo* tch_  = nullptr; // info about track-calo hit

    //----------------------------------------------
    // Track information accessors

    //----------------------------------------------
    // Basic track fit checks
    bool IsGood() {
      if(!track_) return false;
      if(!track_->trk) return false;
      if(track_->trk->status < 0) return false;
      if(track_->trk->goodfit == 0) return false;
      return true;
    }

    //----------------------------------------------
    // Particle hypothesis used in the fit
    int FitPDG() {
      if(!track_ || !track_->trk) return 0;
      return track_->trk->pdg;
    }

    //----------------------------------------------
    // Tracker front segment info
    mu2e::TrkSegInfo* FrontSeg() {
      if(!track_ || !track_->trksegs) return nullptr;
      for(auto& seg : *(track_->trksegs)) {
        if(seg.sid == mu2e::SurfaceIdDetail::TT_Front) return &seg;
      }
      return nullptr;
    }

    //----------------------------------------------
    // Reset the input info
    void Reset() {
      track_   = nullptr;
      info_    = nullptr;
      tch_     = nullptr;
    }

    Track_t() { Reset(); }
  };
}
#endif
