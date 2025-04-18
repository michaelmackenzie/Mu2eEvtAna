//
// A single Track
// Michael MacKenzie (2025)

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
    // Sgment info
    mu2e::TrkSegInfo* Segment(mu2e::SurfaceIdDetail::enum_type surface) {
      if(!track_ || !track_->trksegs) return nullptr;
      for(auto& seg : *(track_->trksegs)) {
        if(seg.sid == surface) return &seg;
      }
      return nullptr;
    }

    //----------------------------------------------
    // Tracker front segment info
    mu2e::TrkSegInfo* FrontSeg() {
      return Segment(mu2e::SurfaceIdDetail::TT_Front);
    }

    //----------------------------------------------
    // Track kinematics at a given segment
    float PSegment     (mu2e::SurfaceIdDetail::enum_type surface) { auto seg = Segment(surface); return (seg) ? seg->mom.r()               :  0.; }
    float PTSegment    (mu2e::SurfaceIdDetail::enum_type surface) { auto seg = Segment(surface); return (seg) ? seg->mom.rho()             :  0.; }
    float PZSegment    (mu2e::SurfaceIdDetail::enum_type surface) { auto seg = Segment(surface); return (seg) ? seg->mom.z()               :  0.; }
    float DMomSegment  (mu2e::SurfaceIdDetail::enum_type surface) { auto seg = Segment(surface); return (seg) ? seg->dmom                  :  0.; }
    float MomErrSegment(mu2e::SurfaceIdDetail::enum_type surface) { auto seg = Segment(surface); return (seg) ? seg->momerr                : -1.; }
    float TSegment     (mu2e::SurfaceIdDetail::enum_type surface) { auto seg = Segment(surface); return (seg) ? seg->time                  :  0.; }

    //----------------------------------------------
    // Track kinematics at the tracker front
    float PFront     () { return PSegment     (mu2e::SurfaceIdDetail::TT_Front); }
    float PTFront    () { return PTSegment    (mu2e::SurfaceIdDetail::TT_Front); }
    float PZFront    () { return PZSegment    (mu2e::SurfaceIdDetail::TT_Front); }
    float DMomFront  () { return DMomSegment  (mu2e::SurfaceIdDetail::TT_Front); }
    float MomErrFront() { return MomErrSegment(mu2e::SurfaceIdDetail::TT_Front); }
    float TFront     () { return TSegment     (mu2e::SurfaceIdDetail::TT_Front); }

    //----------------------------------------------
    // Reset the input info
    void Reset() {
      track_   = nullptr;
      info_    = nullptr;
      tch_     = nullptr;
    }

    //----------------------------------------------
    // Print the track
    void Print(TString opt = "") {
      opt.ToLower();
      if(opt.Contains("banner")) {
        std::string filler(100, '-');
        printf("%s\n", filler.c_str());
        printf("Idx: %5s %10s %10s %10s %10s %5s\n", "Hyp", "p", "pT", "pz", "t", "good");
        printf("%s\n", filler.c_str());
      }
      printf("Idx: %5i %10.2f %10.2f %10.2f %10.1f %5i\n", FitPDG(), PFront(), PTFront(), PZFront(), TFront(), IsGood());
    }

    Track_t() { Reset(); }
  };
}
#endif
