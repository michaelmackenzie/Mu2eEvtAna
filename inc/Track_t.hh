//
// A single Track
// Michael MacKenzie (2025)

#ifndef MU2EEVTANA_TRACK_T_HH
#define MU2EEVTANA_TRACK_T_HH

// Offline includes
#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"

// ROOT includes
#include "Rtypes.h"
#include "TString.h"

// Event Ntuple includes
#include "EventNtuple/inc/TrkInfo.hh"
#include "EventNtuple/inc/TrkCaloHitInfo.hh"
#include "EventNtuple/inc/TrkSegInfo.hh"
#include "EventNtuple/inc/LoopHelixInfo.hh"
#include "EventNtuple/inc/SimInfo.hh"

#include "EventNtuple/rooutil/inc/Track.hh"

// local includes
#include "Mu2eEvtAna/inc/GlobalConstants.h"
#include "Mu2eEvtAna/inc/CRVCluster_t.hh"

namespace Mu2eEvtAna {
  struct Track_t {

    rooutil::Track* track_      = nullptr; // pointer to the EventNtuple::Track object
    CRVCluster_t* stub_         = nullptr; // best matched CRV cluster
    Track_t*      upstream_     = nullptr; // pointer to a matched upstream track

    // Local MVA scores
    float trkqual_;
    float pid_;
    float trkpid_;
    float cosmic_id_;

    // Evaluated values
    float tz_slope_;
    float tz_slope_unc_;

    // Track IDs
    int id_[kMaxTrackIDs];

    // Fit observables
    double obs_[kMaxObservables];

    //----------------------------------------------
    // Track information accessors

    //----------------------------------------------
    // Basic track fit checks
    bool IsGood() const {
      if(!track_) return false;
      if(!track_->trk) return false;
      if(track_->trk->status < 0) return false;
      if(track_->trk->goodfit == 0) return false;
      return true;
    }

    //----------------------------------------------
    // Track fit info
    float Chi2Dof   () const { return (track_ && track_->trk && track_->trk->ndof > 0) ? track_->trk->chisq / track_->trk->ndof : -1.f; }
    float FitCon    () const { return (track_ && track_->trk) ? track_->trk->fitcon : -1.f; }
    int   NActive   () const { return (track_ && track_->trk) ? track_->trk->nactive : -1; }
    int   NHits     () const { return (track_ && track_->trk) ? track_->trk->nhits : -1; }
    int   NNull     () const { return (track_ && track_->trk) ? track_->trk->nnullambig : -1; }
    int   NMatActive() const { return (track_ && track_->trk) ? track_->trk->nmatactive : -1; }
    float TrkQual   () const { return (track_ && track_->trkqual && track_->trkqual->valid) ? track_->trkqual->result : -1000.f; }
    float PID       () const { return -1000.f; } // (track_ && track_->trkpid && track_->trkpid->valid) ? track_->trkpid->result : -1000.f; }
    float AltTrkQual() const { return (track_) ? trkqual_      : -1000.f; }
    float AltPID    () const { return (track_) ? pid_          : -1000.f; }
    float TrkPID    () const { return (track_) ? trkpid_       : -1000.f; }
    float CosmicID  () const { return (track_) ? cosmic_id_    : -1000.f; }
    float TZSlope   () const { return (track_) ? tz_slope_     :     0.f; }
    float TZSlopeUnc() const { return (track_) ? tz_slope_unc_ :    -1.f; }
    bool  OPAInter  () const { return (track_ && track_->trk) ? track_->trk->opainter : false; }
    int   NSTInter  () const {
      return (track_ && track_->trk) ? track_->trk->nstup + track_->trk->nstdown : 0;
    }

    // Static functions
    static float Velocity(double p, double m) {
      if(p < 0. || m < 0. || (p <= 0. && m <= 0.))
        return 0.;
      return p / std::sqrt(p * p + m * m);
    }

    float TZSlopeSig() const {
      if(!track_) return 0.;
      const float slope = TZSlope();
      const float unc   = TZSlopeUnc();
      return (unc <= 0.) ? 0. : slope/unc;
    }

    float TZSlopeRatio() const {
      if(!track_) return 0.;
      const int pdg = std::abs(FitPDG());
      double mass = 0.511; // default to electron
      if     (pdg ==   13) mass = 105.66;
      else if(pdg ==  211) mass = 139.57;
      else if(pdg == 2212) mass = 938.27;
      constexpr float v_light = 300.f; // ~300 mm / ns
      float p = PMiddle();
      float cz = CosThetaMiddle();
      if(p <= 0.) {  // use front as first backup
        p = PFront();
        cz = CosThetaFront();
      }
      if(p <= 0.) {  // use back as last backup
        p = PBack();
        cz = CosThetaBack();
      }
      const float vz = v_light*Velocity(p,mass)*cz;
      const float ratio = TZSlope() * vz;
      return ratio;
    }

    // Evaluate the TZ slope from the track hits and store it
    void EvaluateTZSlope() {
      tz_slope_ = 0.;
      tz_slope_unc_ = -1.;
      if(!track_) return;
      if(!track_->trkhits || track_->trkhits->empty()) return;
      ::LsqSums2 fitDtDz;
      const size_t nhits = track_->trkhits->size();
      for(size_t index = 0; index < nhits; ++index) {
        const auto& hit = track_->trkhits->at(index);
        const double t = hit.etime[hit.earlyend] - hit.tottdrift; // use less biased time
        // const double t = hit.ptoca; // from the track fit
        const double z = hit.poca.z();
        //   double timeErrSquared = hit->timeVar();//ns^2
        //   double hitWeight      = 1./timeErrSquared;
        fitDtDz.addPoint(z, t, 1.);
      }
      tz_slope_ = fitDtDz.dydx();
      // Scale up the errors by chisq / dof to get an estimate of what the averged weight should have been
      const float chisq = fitDtDz.chi2Dof();
      tz_slope_unc_ = fitDtDz.dydxErr()*chisq;
      // Chi2ndof = fitDtDz.chi2Dof();
    }

    //----------------------------------------------
    // Particle hypothesis used in the fit
    int FitPDG() const {
      if(!track_ || !track_->trk) return 0;
      return track_->trk->pdg;
    }
    int Charge() const {
      const int pdg = FitPDG();
      const int abs_pdg = std::abs(pdg);
      const bool negative = std::signbit(pdg); // true if negative
      switch(abs_pdg) {
      case 11: case 13: case 15:
        return (negative) ?  1 : -1;
      case 211: case 2212:
        return (negative) ? -1 :  1;
      case 0: return 0;
      default:
        break;
      }
      printf("Track_t::%s: Unknown fit hypothesis %i --> Returning 0 charge\n", __func__, pdg);
      return 0;
    }

    //----------------------------------------------
    // MC Particle info
    const mu2e::SimInfo* SimInfo() const {
      if(!track_) return nullptr;
      if(!track_->trkmcsim) return nullptr;

      // Find the sim particle with the most active hits
      int max_hits(-1);
      const mu2e::SimInfo* sim_info(nullptr);
      for(const auto& info : *(track_->trkmcsim)) {
        if(info.nactive > max_hits) {
          max_hits = info.nactive;
          sim_info = &info;
        }
      }
      return sim_info;
    }

    int MCPDG    () const { auto sim_info = SimInfo(); return (sim_info) ? sim_info->pdg            : 0; }
    int MCHits   () const { auto sim_info = SimInfo(); return (sim_info) ? sim_info->nhits          : 0; }
    int MCActive () const { auto sim_info = SimInfo(); return (sim_info) ? sim_info->nactive        : 0; }
    int MCProcess() const { auto sim_info = SimInfo(); return (sim_info) ? sim_info->startCode      : 0; }
    int MCGenP   () const { auto sim_info = SimInfo(); return (sim_info) ? sim_info->mom.r()        : 0; }
    int MCGenE   () const {
      auto sim_info = SimInfo();
      if(!sim_info) return 0.;
      // Retrieve the particle mass and add it to the momentum (if available)
      const double mass((std::abs(sim_info->pdg) > 10000) ? 0. : ParticleMass(sim_info->pdg));
      return std::sqrt(std::pow(sim_info->mom.r(), 2) + mass*mass);
    }

    //----------------------------------------------
    // Track-Calo hit retrieval
    const mu2e::TrkCaloHitInfo* TCH() const { return (track_) ? track_->trkcalohit : nullptr; }

    //----------------------------------------------
    // Segment info
    const mu2e::TrkSegInfo* Segment(mu2e::SurfaceIdDetail::enum_type surface) const {
      if(!track_ || !track_->trksegs) return nullptr;
      for(const auto& seg : *(track_->trksegs)) {
        if(seg.sid == surface) return &seg;
      }
      return nullptr;
    }
    int SegmentIndex(mu2e::SurfaceIdDetail::enum_type surface) const {
      if(!track_ || !track_->trksegs) return -1;
      for(size_t index = 0; index < track_->trksegs->size(); ++index) {
        const auto& seg =  track_->trksegs->at(index);
        if(seg.sid == surface) return int(index);
      }
      return -1;
    }
    const mu2e::LoopHelixInfo* LHSegment(mu2e::SurfaceIdDetail::enum_type surface) const {
      if(!track_ || !track_->trksegpars_lh) return nullptr;
      const int index = SegmentIndex(surface);
      if(index < 0 || index >= int(track_->trksegpars_lh->size())) return nullptr;
      return &(track_->trksegpars_lh->at(index));
    }
    const mu2e::SurfaceStepInfo* MCSegment(mu2e::SurfaceIdDetail::enum_type surface) const {
      if(!track_ || !track_->trksegsmc) return nullptr; // No MC info
      auto reco_seg = Segment(surface);
      if(!reco_seg) return nullptr; // Can't do anything if there's no reco info

      // Search for an MC segment matched to this reco segment
      const mu2e::SurfaceStepInfo* seg(nullptr);
      for(const auto& mc_seg : *(track_->trksegsmc)) {
        if(mc_seg.sid != surface) continue;
        // Found this MC surface
        if(seg) { // check which is segment is a better match
          const float dt_curr = std::fabs(reco_seg->time - seg->time);
          const float dt_new  = std::fabs(reco_seg->time - mc_seg.time);
          if(dt_curr > dt_new) seg = &mc_seg;
        } else seg = &mc_seg;
      }
      return seg;
    }

    //----------------------------------------------
    // Tracker front segment info
    const mu2e::TrkSegInfo* FrontSeg() const { return Segment(mu2e::SurfaceIdDetail::TT_Front); }

    //----------------------------------------------
    // Track kinematics at a given surface
    float PSegment      (mu2e::SurfaceIdDetail::enum_type surface) const { auto seg =   Segment(surface); return (seg) ? seg->mom.r()               :  0.; }
    float PTSegment     (mu2e::SurfaceIdDetail::enum_type surface) const { auto seg =   Segment(surface); return (seg) ? seg->mom.rho()             : -1.; }
    float PZSegment     (mu2e::SurfaceIdDetail::enum_type surface) const { auto seg =   Segment(surface); return (seg) ? seg->mom.z()               :  0.; }
    float DMomSegment   (mu2e::SurfaceIdDetail::enum_type surface) const { auto seg =   Segment(surface); return (seg) ? seg->dmom                  :  0.; }
    float MomErrSegment (mu2e::SurfaceIdDetail::enum_type surface) const { auto seg =   Segment(surface); return (seg) ? seg->momerr                : -1.; }
    float TSegment      (mu2e::SurfaceIdDetail::enum_type surface) const { auto seg =   Segment(surface); return (seg) ? seg->time                  :  0.; }
    float TErrSegment   (mu2e::SurfaceIdDetail::enum_type surface) const { auto seg = LHSegment(surface); return (seg) ? seg->t0err                 : -1.; }
    float RMaxSegment   (mu2e::SurfaceIdDetail::enum_type surface) const { auto seg = LHSegment(surface); return (seg) ? seg->maxr                  :  0.; }
    float RadiusSegment (mu2e::SurfaceIdDetail::enum_type surface) const { auto seg = LHSegment(surface); return (seg) ? std::fabs(seg->rad)        :  0.; }
    float RMinSegment   (mu2e::SurfaceIdDetail::enum_type surface) const { return std::fabs(RMaxSegment(surface) - 2.f*RadiusSegment(surface)); }

    float D0Segment(mu2e::SurfaceIdDetail::enum_type surface) const {
      auto seg = LHSegment(surface);
      // return (seg) ? seg->d0 : -1.e6;
      if(!seg) return -1.e6;
      // FIXME: Evaluating this locally to get the sign
      const double radius = std::fabs(seg->rad);
      const double max_r  = seg->maxr;
      return max_r - 2.*radius;
    }

    float TanDipSegment(mu2e::SurfaceIdDetail::enum_type surface) const {
      auto seg = LHSegment(surface);
      if(seg) return seg->tanDip;
      // Estimate it by hand if the segment is not available
      // FIXME: Check that this definition is correct
      // tan(theta) = x/y = pz / pt
      const float pt = PTSegment(surface);
      const float pz = PZSegment(surface);
      if(pt <= 0. || pz == 0.) return -100.;
      return pz/pt;
    }

    float CosThetaSegment(mu2e::SurfaceIdDetail::enum_type surface) const {
      // cos(theta) = x/r = pz / p
      const float p  = PSegment (surface);
      const float pz = PZSegment(surface);
      if(p <= 0. || pz == 0.) return -100.;
      return pz/p;
    }

    //----------------------------------------------
    // MC Track kinematics at a given surface
    float MCPSegment     (mu2e::SurfaceIdDetail::enum_type surface) const { auto seg = MCSegment(surface); return (seg) ? seg->mom.r()               :  0.; }
    float MCPTSegment    (mu2e::SurfaceIdDetail::enum_type surface) const { auto seg = MCSegment(surface); return (seg) ? seg->mom.rho()             : -1.; }
    float MCPZSegment    (mu2e::SurfaceIdDetail::enum_type surface) const { auto seg = MCSegment(surface); return (seg) ? seg->mom.z()               :  0.; }
    float MCTSegment     (mu2e::SurfaceIdDetail::enum_type surface) const { auto seg = MCSegment(surface); return (seg) ? seg->time                  :  0.; }

    //----------------------------------------------
    // Track kinematics at the tracker front
    float PFront        () const { return PSegment       (mu2e::SurfaceIdDetail::TT_Front); }
    float PTFront       () const { return PTSegment      (mu2e::SurfaceIdDetail::TT_Front); }
    float PZFront       () const { return PZSegment      (mu2e::SurfaceIdDetail::TT_Front); }
    float DMomFront     () const { return DMomSegment    (mu2e::SurfaceIdDetail::TT_Front); }
    float MomErrFront   () const { return MomErrSegment  (mu2e::SurfaceIdDetail::TT_Front); }
    float TFront        () const { return TSegment       (mu2e::SurfaceIdDetail::TT_Front); }
    float TErrFront     () const { return TErrSegment    (mu2e::SurfaceIdDetail::TT_Front); }
    float D0Front       () const { return D0Segment      (mu2e::SurfaceIdDetail::TT_Front); }
    float TanDipFront   () const { return TanDipSegment  (mu2e::SurfaceIdDetail::TT_Front); }
    float CosThetaFront () const { return CosThetaSegment(mu2e::SurfaceIdDetail::TT_Front); }
    float RMaxFront     () const { return RMaxSegment    (mu2e::SurfaceIdDetail::TT_Front); }
    float RadiusFront   () const { return RadiusSegment  (mu2e::SurfaceIdDetail::TT_Front); }
    int   Trajectory    () const {
      const float pz = PZFront();
      if(pz == 0.f) return 0;
      return (pz < 0.f) ? -1 : 1;
    }

    float MCPFront     () const { return MCPSegment   (mu2e::SurfaceIdDetail::TT_Front); }
    float MCPTFront    () const { return MCPTSegment  (mu2e::SurfaceIdDetail::TT_Front); }
    float MCPZFront    () const { return MCPZSegment  (mu2e::SurfaceIdDetail::TT_Front); }
    float MCTFront     () const { return MCTSegment   (mu2e::SurfaceIdDetail::TT_Front); }
    float MCDeltaPFront() const { return PFront() - MCPFront(); }
    int   MCTrajectory () const {
      const float pz = MCPZFront();
      if(pz == 0.f) return 0;
      return (pz < 0.f) ? -1 : 1;
    }

    //----------------------------------------------
    // Track kinematics at the tracker middle
    float PMiddle       () const { return PSegment       (mu2e::SurfaceIdDetail::TT_Mid); }
    float PTMiddle      () const { return PTSegment      (mu2e::SurfaceIdDetail::TT_Mid); }
    float PZMiddle      () const { return PZSegment      (mu2e::SurfaceIdDetail::TT_Mid); }
    float DMomMiddle    () const { return DMomSegment    (mu2e::SurfaceIdDetail::TT_Mid); }
    float MomErrMiddle  () const { return MomErrSegment  (mu2e::SurfaceIdDetail::TT_Mid); }
    float TMiddle       () const { return TSegment       (mu2e::SurfaceIdDetail::TT_Mid); }
    float TErrMiddle    () const { return TErrSegment    (mu2e::SurfaceIdDetail::TT_Mid); }
    float D0Middle      () const { return D0Segment      (mu2e::SurfaceIdDetail::TT_Mid); }
    float TanDipMiddle  () const { return TanDipSegment  (mu2e::SurfaceIdDetail::TT_Mid); }
    float CosThetaMiddle() const { return CosThetaSegment(mu2e::SurfaceIdDetail::TT_Mid); }
    float RMaxMiddle    () const { return RMaxSegment    (mu2e::SurfaceIdDetail::TT_Mid); }

    //----------------------------------------------
    // Track kinematics at the tracker back
    float PBack       () const { return PSegment       (mu2e::SurfaceIdDetail::TT_Back); }
    float PTBack      () const { return PTSegment      (mu2e::SurfaceIdDetail::TT_Back); }
    float PZBack      () const { return PZSegment      (mu2e::SurfaceIdDetail::TT_Back); }
    float DMomBack    () const { return DMomSegment    (mu2e::SurfaceIdDetail::TT_Back); }
    float MomErrBack  () const { return MomErrSegment  (mu2e::SurfaceIdDetail::TT_Back); }
    float TBack       () const { return TSegment       (mu2e::SurfaceIdDetail::TT_Back); }
    float TErrBack    () const { return TErrSegment    (mu2e::SurfaceIdDetail::TT_Back); }
    float D0Back      () const { return D0Segment      (mu2e::SurfaceIdDetail::TT_Back); }
    float TanDipBack  () const { return TanDipSegment  (mu2e::SurfaceIdDetail::TT_Back); }
    float CosThetaBack() const { return CosThetaSegment(mu2e::SurfaceIdDetail::TT_Back); }
    float RMaxBack    () const { return RMaxSegment    (mu2e::SurfaceIdDetail::TT_Back); }

    //----------------------------------------------
    // Track kinematics at the stopping target exit
    float PSTBack     () const { return PSegment     (mu2e::SurfaceIdDetail::ST_Back); }
    float PTSTBack    () const { return PTSegment    (mu2e::SurfaceIdDetail::ST_Back); }
    float PZSTBack    () const { return PZSegment    (mu2e::SurfaceIdDetail::ST_Back); }
    float DMomSTBack  () const { return DMomSegment  (mu2e::SurfaceIdDetail::ST_Back); }
    float MomErrSTBack() const { return MomErrSegment(mu2e::SurfaceIdDetail::ST_Back); }
    float TSTBack     () const { return TSegment     (mu2e::SurfaceIdDetail::ST_Back); }
    float TErrSTBack  () const { return TErrSegment  (mu2e::SurfaceIdDetail::ST_Back); }
    float D0STBack    () const { return D0Segment    (mu2e::SurfaceIdDetail::ST_Back); }
    float TanDipSTBack() const { return TanDipSegment(mu2e::SurfaceIdDetail::ST_Back); }
    float RMaxSTBack  () const { return RMaxSegment  (mu2e::SurfaceIdDetail::ST_Back); }

    float MCPSTBack     () const { return MCPSegment   (mu2e::SurfaceIdDetail::ST_Back); }
    float MCPTSTBack    () const { return MCPTSegment  (mu2e::SurfaceIdDetail::ST_Back); }
    float MCPZSTBack    () const { return MCPZSegment  (mu2e::SurfaceIdDetail::ST_Back); }
    float MCTSTBack     () const { return MCTSegment   (mu2e::SurfaceIdDetail::ST_Back); }
    float MCDeltaPSTBack() const { return PSTBack() - MCPSTBack(); }

    //----------------------------------------------
    // Track CaloHit info
    float ECluster() const { auto tch = TCH(); return (!tch) ?    0. : tch->edep  ; }
    float Dt      () const { auto tch = TCH(); return (!tch) ? -1.e6 : tch->dt    ; }
    float DOCA    () const { auto tch = TCH(); return (!tch) ? -1.e6 : tch->doca  ; }
    float PTOCA   () const { auto tch = TCH(); return (!tch) ? -1.e6 : tch->ptoca ; }
    float CDepth  () const { auto tch = TCH(); return (!tch) ? -1.e6 : tch->cdepth; }
    float EPFront () const {
      const float ecl(ECluster()), p(PFront());
      return (p > 0.) ? ecl/p : 0.;
    }

    //----------------------------------------------
    // Accessing/setting the track IDs
    void SetID(const int ID, const int index = 0) {
      if(index < 0 || index >= kMaxTrackIDs) throw std::runtime_error(Form("Accessing a track ID index (%i) out of bounds!", index));
      id_[index] = ID;
    }
    int ID(const int index = 0) const {
      if(index < 0 || index >= kMaxTrackIDs) throw std::runtime_error(Form("Accessing track ID index (%i) out of bounds!", index));
      return id_[index];
    }

    //----------------------------------------------
    // Accessing/setting the observable
    void SetObs(const double val, const int index = 0) {
      if(index < 0 || index >= kMaxObservables) throw std::runtime_error(Form("Accessing an observable index (%i) out of bounds!", index));
      obs_[index] = val;
    }
    double Obs(const int index = 0) const {
      if(index < 0 || index >= kMaxObservables) throw std::runtime_error(Form("Accessing an observable index (%i) out of bounds!", index));
      return obs_[index];
    }

    //----------------------------------------------
    // Reset the input info
    void Reset() {
      track_        = nullptr;
      stub_         = nullptr;
      upstream_     = nullptr;
      trkqual_      = -1000.f;
      pid_          = -1000.f;
      trkpid_       = -1000.f;
      cosmic_id_    = -1000.f;
      tz_slope_     = 0.f;
      tz_slope_unc_ = -1.f;
      for(int iid = 0; iid < kMaxTrackIDs; ++iid) id_[iid] = 0;
      for(int iobs = 0; iobs < kMaxObservables; ++iobs) obs_[iobs] = 0.;
    }

    //----------------------------------------------
    // Print the track
    void Print(TString opt = "") const {
      opt.ToLower();
      if(opt.Contains("banner")) {
        std::string filler(130, '-');
        printf("%s\n", filler.c_str());
        printf("Idx: %5s %10s %10s %10s %10s %7s %6s %10s %5s %5s %8s %8s\n", "Hyp", "p", "pT", "pz", "t", "Ecl", "tandip", "p(MC)", "PDG", "good", "fitcon", "trkqual");
        printf("%s\n", filler.c_str());
      }
      if(!track_) return;
      printf("Idx: %5i %10.2f %10.2f %10.2f %10.1f %7.1f %6.2f %10.1f %5i %5i %.2e %8.5f\n", FitPDG(), PFront(), PTFront(), PZFront(), TFront(), ECluster(), TanDipFront(),
             MCPFront(), MCPDG(), IsGood(), FitCon(), TrkQual());
    }

    Track_t() { Reset(); }
  };
}
#endif
