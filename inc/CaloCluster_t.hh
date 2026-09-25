//
// Calo cluster information
// Michael MacKenzie (2025)

#ifndef MU2EEVTANA_CALOLUSTER_T_HH
#define MU2EEVTANA_CALOLUSTER_T_HH

// ROOT includes
#include "Rtypes.h"

// EventNtuple includes
#include "EventNtuple/rooutil/inc/CaloCluster.hh"
#include "EventNtuple/rooutil/inc/CaloHit.hh"
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
    const rooutil::CaloCluster* cc_;

    // Best-matched objects for this cluster (nullptr if none/not searched for). Populated by an
    // analysis module, not by the base Mu2eEvtAna -- e.g. Run1BAna::MatchCaloClusters().
    Track_t*       line_        ; // best matched track (line fit)
    LineSeed_t*    line_seed_   ; // best matched line seed
    TimeCluster_t* time_cluster_; // best matched time cluster, from any collection
    // Best matched time cluster from the nominal (target-origin electron) collection specifically.
    // time_cluster_ above searches the collections in whatever order they were discovered, so it
    // is not necessarily from the nominal one; a selection that depends on which finder produced
    // the cluster wants this instead.
    TimeCluster_t* nom_time_cluster_;
    CRVCluster_t*  crv_cluster_ ; // best matched CRV coincidence cluster

    // Derived hit-based quantities, cached by Init() so repeated accessor calls don't
    // re-walk the hit collection. Indexed by usewt for the (un)weighted variants.
    bool  is_init_        ;
    float tmean_[2]       ; // [0] = unweighted, [1] = energy-weighted
    float tvar_ [2]       ; // [0] = unweighted, [1] = energy-weighted
    float hit_max_r_      ;
    float hit_max_extent_ ;
    float e_out_ring_     ; // energy beyond 600 mm
    float e9_             ; // energy of 3x3 around the main hit
    float e25_            ; // energy of 5x5 around the main hit
    float second_moment_  ; // hit position second moment, weighted by hit energy

    // crystal size
    constexpr static float crystal_dx_ = 34.; // crystal width


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
    // Hit information

    const rooutil::CaloHits& Hits() const { return cc_->GetHits(); }
    const rooutil::CaloHit& Hit(size_t index) const { return Hits().at(index); }
    bool ValidHit(size_t index) const {
      if(!cluster_ || index >= Hits().size()) return false;
      const auto& hit = Hit(index);
      if(!hit.reco) return false;
      return true;
    }

    float HitE(size_t index) const {
      if(!ValidHit(index)) return 0.f;
      const auto& hit = Hit(index);
      return hit.reco->eDep_;
    }
    float HitT(size_t index) const {
      if(!ValidHit(index)) return 0.f;
      const auto& hit = Hit(index);
      return hit.reco->time_;
    }
    XYZVectorF HitPos(size_t index) const {
      if(!ValidHit(index)) return  XYZVectorF();
      const auto& hit = Hit(index);
      return hit.reco->crystalPos_;
    }
    float HitR(size_t index) const {
      if(!ValidHit(index)) return 0.f;
      const auto& pos = HitPos(index);
      const float x = pos.x();
      const float y = pos.y();
      return std::sqrt(x*x+y*y);
    }
    float HitMCEDep(size_t index) const {
      if(!ValidHit(index)) return 0.f;
      const auto& hit = Hit(index);
      if(!hit.mc) return 0.f;
      return hit.mc->eDep;
    }
    const std::vector<int>& HitMCSimIDs(size_t index) const {
      static const std::vector<int> vec;
      if(!ValidHit(index)) return vec;
      const auto& hit = Hit(index);
      if(!hit.mc) {
        return vec;
      }
      return hit.mc->simParticleIds;
    }
    float HitSimEDep(size_t index, int sim_id) const {
      if(!ValidHit(index)) return 0.f;
      const auto& hit = Hit(index);
      if(!hit.mc) return 0.f;
      const auto& ids = hit.mc->simParticleIds;
      const size_t nsims = ids.size();
      float edep = 0.f;
      for(size_t sim_index = 0; sim_index < nsims; ++sim_index) {
        if(sim_id == ids[sim_index]) edep += hit.mc->eDeps[sim_index];
      }
      return edep;
    }

    //-------------------------------------------------
    // Additional functions

    float R        () const { return std::sqrt(X()*X() + Y()*Y()); }
    float E1       () const { return HitE(0); }
    float E2       () const { return E1() + HitE(1); }
    float TMean    (bool usewt = false) const { // mean hit time, with/without energy weights
      if(is_init_) return tmean_[usewt ? 1 : 0];
      return CalcTMean(usewt);
    }
    float TVar     (bool usewt = false) const { // hit time variance, with/without energy weights
      if(is_init_) return tvar_[usewt ? 1 : 0];
      return CalcTVar(usewt);
    }
    float MaxHitR() const { // radial position of highest radius hit
      if(is_init_) return hit_max_r_;
      return CalcMaxHitR();
    }
    float MaxHitExtent() const { // distance from the main hit and the farthest hit
      if(is_init_) return hit_max_extent_;
      return CalcMaxHitExtent();
    }
    float EOutRing() const { // energy at high radius (uses default value)
      if(is_init_) return e_out_ring_;
      return CalcEOutRing();
    }
    float E9() const {
      if(is_init_) return e9_;
      return CalcENeighbors(1.5*crystal_dx_);
    }
    float E25() const {
      if(is_init_) return e25_;
      return CalcENeighbors(2.5*crystal_dx_);
    }
    float ERing() const { return E9() - E1(); }
    float SecondMoment() const {
      if(is_init_) return second_moment_;
      return CalcSecondMoment();
    }

    //-------------------------------------------------
    // Uncached implementations of the derived hit-based quantities above, used both as the
    // fallback when the cluster hasn't been Init()-ed and to fill the cache inside Init().

    float CalcTMean(bool usewt) const {
      if(!cc_) return 0.f;
      const size_t nhits = Hits().size();
      if(nhits == 0) return 0.f;
      float tmean = 0.f;
      float sumwt = 0.f;
      for(size_t ihit = 0; ihit < nhits; ++ihit) {
        const float wt = (usewt) ? HitE(ihit) : 1.f;
        tmean += wt*HitT(ihit);
        sumwt += wt;
      }
      tmean /= sumwt;
      return tmean;
    }
    float CalcTVar(bool usewt) const {
      if(!cc_) return -1.f;
      const size_t nhits = Hits().size();
      if(nhits == 0) return 0.f;
      const float tmean = CalcTMean(false);
      float var = 0.f;
      float sumwt = 0.f;
      for(size_t ihit = 0; ihit < nhits; ++ihit) {
        const float dt = HitT(ihit) - tmean;
        const float wt = (usewt) ? HitE(ihit) : 1.f;
        var += wt*dt*dt;
        sumwt += wt;
      }
      var /= sumwt;
      return var;
    }
    float CalcMaxHitR() const {
      if(!cc_) return 0.f;
      const size_t nhits = Hits().size();
      if(nhits == 0) return 0.f;
      float max_r = 0.f;
      for(size_t ihit = 0; ihit < nhits; ++ihit) {
        const float r = HitR(ihit);
        max_r = std::max(r, max_r);
      }
      return max_r;
    }
    float CalcMaxHitExtent() const {
      if(!cc_) return 0.f;
      const size_t nhits = Hits().size();
      if(nhits <= 1) return 0.f;
      const auto main_pos = HitPos(0);
      float max_r = 0.f;
      for(size_t ihit = 1; ihit < nhits; ++ihit) {
        const auto hit_pos = HitPos(ihit);
        const float r = (main_pos - hit_pos).r();
        max_r = std::max(r, max_r);
      }
      return max_r;
    }
    float CalcEOutRing(const float r_cut = 600.) const {
      if(!cc_) return 0.f;
      const size_t nhits = Hits().size();
      if(nhits == 0) return 0.f;
      float energy = 0.f;
      for(size_t ihit = 1; ihit < nhits; ++ihit) {
        const float r = HitR(ihit);
        if(r > r_cut) energy += HitE(ihit);
      }
      return energy;
    }
    float CalcENeighbors(const float dx) const { // calculate energy within a (dx, dx) matrix around the main hit
      if(!cc_) return 0.f;
      const size_t nhits = Hits().size();
      if(nhits == 0) return 0.f;
      const auto& main_pos = HitPos(0);
      float energy = HitE(0);
      for(size_t ihit = 1; ihit < nhits; ++ihit) {
        const auto hit_pos = HitPos(ihit);
        const auto dp = hit_pos - main_pos;
        // within the (dx, dx) matrix
        if(std::fabs(dp.x()) < dx && std::fabs(dp.y()) < dx)
          energy += HitE(ihit);
      }
      return energy;
    }
    float CalcSecondMoment() const { // calculate energy second moment
      if(!cc_) return 0.f;
      const size_t nhits = Hits().size();
      if(nhits == 0) return 0.f;
      double sx(0.),sy(0.),sx2(0.),sy2(0.),sw(0.);
      for(size_t ihit = 0; ihit < nhits; ++ihit) {
        const double energy(HitE(ihit));
        const auto pos = HitPos(ihit);
        double xCrystal = pos.x();
        double yCrystal = pos.y();

        double weight = energy;

        sw  += weight;
        sx  += xCrystal*weight;
        sy  += yCrystal*weight;
        sx2 += xCrystal*xCrystal*weight;
        sy2 += yCrystal*yCrystal*weight;
      }
      const float moment = (sx2-sx*sx/sw + sy2-sy*sy/sw)/sw;
      return moment;
    }

    const mu2e::SimInfo* MCSim() const {
      if(!cluster_mc_) return nullptr;
      if(!cc_) return nullptr;
      const auto& sims = cc_->GetMCParticles();
      if(sims.empty()) return 0;
      const size_t nsims = sims.size();
      if(cluster_mc_->simRelRels.size() != nsims) return nullptr;
      for(size_t index = 0; index < nsims; ++index) {
        if(!sims.at(index).mcsim) continue;
        if(cluster_mc_->simRelRels.at(index) == 0) {
          return sims.at(index).mcsim;
        }
      }
      return nullptr;
    }

    int MCPDG() const {
      const auto sim = MCSim();
      if(!sim) return 0;
      return sim->pdg;
    }
    int MCSimID() const {
      const auto sim = MCSim();
      if(!sim) return -1;
      return sim->id;
    }
    float MCEDep() const {
      if(!cluster_mc_) return 0.f;
      return cluster_mc_->etot;
    }
    float MCTime() const {
      if(!cluster_mc_) return 0.f;
      return cluster_mc_->tavg;
    }

    float MCSimEDep() const {
      if(!cluster_mc_) return 0.f;
      return cluster_mc_->eprimary;
      // const size_t nhits = Hits().size();
      // if(nhits == 0) return 0.f;
      // const auto sim = MCSim();
      // if(!sim) return 0.f;
      // const int sim_id = sim->id;
      // float edep = 0.f;
      // for(size_t ihit = 0; ihit < nhits; ++ihit) {
      //   edep += HitSimEDep(ihit, sim_id);
      // }
      // return edep;
    }

    void Reset() {
      cluster_ = nullptr;
      cluster_mc_ = nullptr;
      cc_ = nullptr;
      line_ = nullptr;
      line_seed_ = nullptr;
      time_cluster_ = nullptr;
      nom_time_cluster_ = nullptr;
      crv_cluster_ = nullptr;

      is_init_ = false;
      tmean_[0] = tmean_[1] = 0.f;
      tvar_ [0] = tvar_ [1] = 0.f;
      hit_max_r_ = 0.f;
      hit_max_extent_ = 0.f;
      e_out_ring_ = 0.f;
      e9_ = 0.f;
      e25_ = 0.f;
      second_moment_ = 0.f;
    }

    // Set the underlying cluster pointers and cache the derived hit-based quantities above
    // (TMean, TVar, MaxHitR, MaxHitExtent) so later accessor calls are O(1).
    void Init(const rooutil::CaloCluster* cluster) {
      Reset();
      if(!cluster) return;
      cluster_ = cluster->calocluster;
      cluster_mc_ = cluster->caloclustermc;
      cc_ = cluster;

      tmean_[0] = CalcTMean(false);
      tmean_[1] = CalcTMean(true);
      tvar_ [0] = CalcTVar(false);
      tvar_ [1] = CalcTVar(true);
      hit_max_r_ = CalcMaxHitR();
      hit_max_extent_ = CalcMaxHitExtent();
      e_out_ring_ = CalcEOutRing();
      e9_ = CalcENeighbors(1.5*crystal_dx_);
      e25_ = CalcENeighbors(2.5*crystal_dx_);
      second_moment_ = CalcSecondMoment();
      is_init_ = true;
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
