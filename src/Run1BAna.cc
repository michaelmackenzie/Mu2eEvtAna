#include "Mu2eEvtAna/inc/Run1BAna.hh"
#include <algorithm>

using namespace mu2e;
namespace Mu2eEvtAna {

  //------------------------------------------------------------------------------------
  // Constructor
  Run1BAna::Run1BAna(int verbose) : Mu2eEvtAna(verbose) {
    // Run1B tracks are KinematicLine fits: keep trksegpars_kl, which the base class disables by
    // default (Run-1A tracks are LoopHelix and don't need it).
    disabled_branches_.erase(
      std::remove(disabled_branches_.begin(), disabled_branches_.end(), TString("trksegpars_kl")),
      disabled_branches_.end());
    // Calo hits are used by default
    disabled_branches_.erase(
      std::remove(disabled_branches_.begin(), disabled_branches_.end(), TString("calohits")),
      disabled_branches_.end());
  }

  //------------------------------------------------------------------------------------
  // Initialize the input ntuple information
  int Run1BAna::InitializeInput() {
    if(Mu2eEvtAna::InitializeInput()) return -1;

    // Discover whatever time cluster / line seed collections this ntuple contains -- Run1B jobs
    // can fill several of each (e.g. timeclusters/protontimeclusters/tztimeclusters/
    // tphitimeclusters and lineseeds/protonlineseeds/cosmiclineseeds), under fhicl-configurable
    // names, so nothing is hardcoded here.
    tc_names_.clear();
    ls_names_.clear();
    for(const auto& name : event_->TimeClusterCollectionNames()) tc_names_.emplace_back(name);
    for(const auto& name : event_->LineSeedCollectionNames())    ls_names_.emplace_back(name);

    if(int(tc_names_.size()) > kMaxTCCollections)
      throw std::runtime_error(Form("Run1BAna::%s: Exceeded the maximum number of time cluster collections (%i)!", __func__, kMaxTCCollections));
    if(int(ls_names_.size()) > kMaxLSCollections)
      throw std::runtime_error(Form("Run1BAna::%s: Exceeded the maximum number of line seed collections (%i)!", __func__, kMaxLSCollections));

    time_clusters_.assign(tc_names_.size(), std::vector<TimeCluster_t>());
    line_seeds_.assign(ls_names_.size(), std::vector<LineSeed_t>());

    if(verbose_ > -1) {
      printf("Run1BAna::%s: Found %zu time cluster collection(s):", __func__, tc_names_.size());
      for(const auto& name : tc_names_) printf(" %s", name.Data());
      printf("\n");
      printf("Run1BAna::%s: Found %zu line seed collection(s):", __func__, ls_names_.size());
      for(const auto& name : ls_names_) printf(" %s", name.Data());
      printf("\n");
    }
    return 0;
  }

  //------------------------------------------------------------------------------------
  // Initialize the output ntuple information
  int Run1BAna::InitializeOutput() {
    Mu2eEvtAna::InitializeOutput();
    return 0;
  }

  //------------------------------------------------------------------------------------
  // Define the histogram selections
  void Run1BAna::InitHistSelections() {
    //-----------------------------------------------------------------------------
    // book histogram selections
    //-----------------------------------------------------------------------------
    struct hist_info_t {
      TString _dsc; // description of the selection
      bool    _cls; // calo cluster histograms
      bool    _trk; // track histograms
      bool    _tcs; // time cluster histograms
      bool    _lns; // line seed histograms
      bool    _smp; // sim particle histograms
      bool    _gnp; // gen particle histograms
      bool    _crv; // CRV histograms
      bool    _crs; // control regions included
      hist_info_t(TString dsc = "", bool cls = false, bool trk = false, bool tcs = false, bool lns = false,
                  bool smp = false, bool gnp = false,
                  bool crv = false, bool crs = false)
        : _dsc(dsc), _cls(cls), _trk(trk), _tcs(tcs), _lns(lns), _smp(smp), _gnp(gnp), _crv(crv),
          _crs(crs) {}
    };

    hist_info_t* hist_sets[kMaxHists];
    for (int i=0; i<kMaxHists; i++) {
      hist_sets[i] = nullptr;
    }

    //                                 description                         cls    trk    tcs    lns    simp   genp   crv    crs
    hist_sets[  0] = new hist_info_t("All events"                      ,  true,  true,  true,  true,  true,  true, false, false);
    hist_sets[  1] = new hist_info_t("E > 50 MeV"                      ,  true,  true,  true,  true,  true,  true, false, false);
    hist_sets[  2] = new hist_info_t("E > 70 MeV"                      ,  true,  true,  true,  true,  true,  true, false, false);

    for (int i=0; i<kMaxHists; i++) {
      const int index = i % 1000; // base index, ignoring control region offset
      if(!hist_sets[index]) continue;
      const bool is_cr = i >= 1000 && !hist_sets[index]->_crs;
      if(is_cr && !hist_sets[index]->_crs) continue; // control region
      if(i >= 4000) break; //Control regions above 4000 not yet implemented
      evt_hists_[i] = new EventHist_t;
      if(hist_sets[index]->_cls) cls_hists_[i] = new CaloClusterHist_t;
      if(hist_sets[index]->_trk) trk_hists_[i] = new TrackHist_t;
      if(hist_sets[index]->_crv && ! is_cr) crv_hists_[i] = new CRVHist_t;
      if(hist_sets[index]->_tcs) tcs_hists_[i] = new TimeClusterHist_t;
      if(hist_sets[index]->_lns) lns_hists_[i] = new LineSeedHist_t;
      // FIXME: Add missing histogram types
    }
  }

  //------------------------------------------------------------------------------------
  // Book the histogram sets
  void Run1BAna::BookHistograms(TDirectory* dir) {
    Mu2eEvtAna::BookHistograms(dir);
  }

  //------------------------------------------------------------------------------------
  // Initialize event information: base class tracks/calo/CRV, plus every discovered time
  // cluster / line seed collection
  void Run1BAna::InitializeEvent() {
    Mu2eEvtAna::InitializeEvent();

    for(size_t icoll = 0; icoll < tc_names_.size(); ++icoll) {
      auto& coll = time_clusters_[icoll];
      coll.clear();
      const auto& clusters = event_->GetTimeClusters(tc_names_[icoll].Data());
      if(int(clusters.size()) > kMaxTimeClusters) {
        throw std::runtime_error(Form("Run1BAna::%s: Exceeded the maximum number of time clusters (%i) in collection %s!",
                                      __func__, kMaxTimeClusters, tc_names_[icoll].Data()));
      }
      for(const auto& cluster : clusters) {
        TimeCluster_t tc;
        tc.cluster_ = cluster.timecluster;
        coll.emplace_back(tc);
      }
    }

    for(size_t icoll = 0; icoll < ls_names_.size(); ++icoll) {
      auto& coll = line_seeds_[icoll];
      coll.clear();
      const auto& seeds = event_->GetLineSeeds(ls_names_[icoll].Data());
      if(int(seeds.size()) > kMaxLineSeeds) {
        throw std::runtime_error(Form("Run1BAna::%s: Exceeded the maximum number of line seeds (%i) in collection %s!",
                                      __func__, kMaxLineSeeds, ls_names_[icoll].Data()));
      }
      for(const auto& seed : seeds) {
        LineSeed_t ls;
        ls.seed_ = seed.lineseed;
        coll.emplace_back(ls);
      }
    }

    // Match every calo cluster to its best line/line seed/time cluster/CRV cluster. Must run
    // last: it needs tracks_, crv_clusters_ (both filled by Mu2eEvtAna::InitializeEvent() above)
    // and time_clusters_/line_seeds_ (just filled above) all populated for this event.
    MatchCaloClusters();
  }

  //------------------------------------------------------------------------------------
  // Populate each calo cluster's best-matched line/line seed/time cluster/CRV cluster.
  //
  // A line, line seed, or time cluster that has an associated calo cluster stores a copy of that
  // cluster's own energy/time (LineSeedInfo::ecalo/tcalo, EventNtupleTimeClusterInfo::ecalo/tcalo,
  // TrkCaloHitInfo::edep/ctime) -- literal copies made when Offline filled that object, not an
  // independent measurement -- so matching on (near-)equality of those copies against a
  // calo_clusters_[] entry's own Energy()/Time() recovers the *same* association Offline already
  // made, rather than reconstructing an approximate one. CRV clusters carry no such calo
  // association at all, so that match instead reuses the extrapolated-time-window approach
  // Mu2eEvtAna::InitTrack() already uses for track-to-CRV matching.
  void Run1BAna::MatchCaloClusters() {
    const float kCaloMatchEps(1.e-3f); // energy [MeV] / time [ns] tolerance for the copied-value match
    const float kCRVMatchWindow(250.f); // ns, matching Mu2eEvtAna::InitTrack()'s CRV time window

    for(int icls = 0; icls < evt_.ncalo_clusters_; ++icls) {
      auto& cls = calo_clusters_[icls];
      const float E(cls.Energy());
      const float T(cls.Time());

      // Best matched line (track)
      cls.line_ = nullptr;
      for(int itrk = 0; itrk < evt_.ntracks_ && !cls.line_; ++itrk) {
        Track_t* trk = &tracks_[itrk];
        if(!trk->IsGood()) continue;
        const auto* tch = trk->TCH();
        if(!tch) continue;
        if(std::fabs(tch->edep - E) < kCaloMatchEps && std::fabs(tch->ctime - T) < kCaloMatchEps) cls.line_ = trk;
      }

      // Best matched line seed (any collection)
      cls.line_seed_ = nullptr;
      for(size_t icoll = 0; icoll < ls_names_.size() && !cls.line_seed_; ++icoll) {
        for(auto& seed : line_seeds_[icoll]) {
          if(!seed.HasCalo()) continue;
          if(std::fabs(seed.ECalo() - E) < kCaloMatchEps && std::fabs(seed.TCalo() - T) < kCaloMatchEps) { cls.line_seed_ = &seed; break; }
        }
      }

      // Best matched time cluster (any collection)
      cls.time_cluster_ = nullptr;
      for(size_t icoll = 0; icoll < tc_names_.size() && !cls.time_cluster_; ++icoll) {
        for(auto& tc : time_clusters_[icoll]) {
          if(!tc.HasCalo()) continue;
          if(std::fabs(tc.ECalo() - E) < kCaloMatchEps && std::fabs(tc.TCalo() - T) < kCaloMatchEps) { cls.time_cluster_ = &tc; break; }
        }
      }

      // Best matched CRV cluster: no direct calo association recorded on either side, so search
      // by extrapolated time instead (same style as the track-to-CRV match in InitTrack()).
      CRVCluster_t* crv_match(nullptr);
      float match_min_dt(1.e10);
      for(int icrv = 0; icrv < evt_.ncrv_clusters_; ++icrv) {
        CRVCluster_t* stub = &crv_clusters_[icrv];
        const float dt = std::fabs(stub->TimeViaCaloFront() - T);
        if(dt < kCRVMatchWindow && (!crv_match || dt < match_min_dt)) { crv_match = stub; match_min_dt = dt; }
      }
      cls.crv_cluster_ = crv_match;
    }
  }

  //------------------------------------------------------------------------------------
  // Find the index of the best-matched line seed (by |T0(seed) - T0(track)|) in collection icoll
  int Run1BAna::BestLineSeed(const Track_t* track, int icoll) const {
    if(!track || icoll < 0 || icoll >= int(line_seeds_.size())) return -1;
    const float trk_time(track->TFront());
    int best(-1);
    float best_dt(1.e10);
    const auto& seeds = line_seeds_[icoll];
    for(size_t iseed = 0; iseed < seeds.size(); ++iseed) {
      const float dt = std::fabs(seeds[iseed].T0() - trk_time);
      if(dt < best_dt) { best_dt = dt; best = iseed; }
    }
    return best;
  }

  //------------------------------------------------------------------------------------
  // Find the index of the best-matched time cluster (by |T0(cluster) - T0(seed)|) in collection icoll
  int Run1BAna::BestTimeCluster(const LineSeed_t* seed, int icoll) const {
    if(!seed || icoll < 0 || icoll >= int(time_clusters_.size())) return -1;
    const float seed_time(seed->T0());
    int best(-1);
    float best_dt(1.e10);
    const auto& clusters = time_clusters_[icoll];
    for(size_t icluster = 0; icluster < clusters.size(); ++icluster) {
      const float dt = std::fabs(clusters[icluster].T0() - seed_time);
      if(dt < best_dt) { best_dt = dt; best = icluster; }
    }
    return best;
  }

  //------------------------------------------------------------------------------------
  // Straight-line track selection
  CutID Run1BAna::TrackID(Track_t* track) {
    if(!track || !track->track_) return 0;
    CutID ID;

    // FIXME: initial, wide windows -- tune once real Run1B efficiency/background studies are
    // available. Replaces the helix RMax window with the fit-agnostic RTrackerMax()/RTrackerMin(),
    // and drops the tan(dip) window entirely (undefined for a field-off/straight-line track).
    if(track->PFront() < 50.f || track->PFront() > 150.f)          ID.SetBit(kP);
    if(track->RTrackerMax() < 300. || track->RTrackerMax() > 700.) ID.SetBit(kRMax);
    if(track->TFront() < 400. || track->TFront() > 1700.)          ID.SetBit(kT0);
    if(track->FitCon() < 1.e-5)                                    ID.SetBit(kFitCon);
    if(track->ECluster() <= 0.)                                    ID.SetBit(kClusterE);

    if(track->OPAInter())                                          ID.SetBit(kRMax);
    if(track->TSDAInter())                                         ID.SetBit(kRMax);

    return ID;
  }

  //------------------------------------------------------------------------------------
  // Main event-by-event processing
  bool Run1BAna::ProcessEvent() {
    cut_flow_.ResetEvent();
    FillEventHist(evt_hists_[0]); // all events with well defined inputs

    for(int icls = 0; icls < evt_.ncalo_clusters_; ++icls) {
      const auto cluster = &calo_clusters_[icls];
      cut_flow_.Increment("has_cluster");
      FillCaloClusterHist(cls_hists_[0], cluster);
      if(cluster->Energy() > 50.) FillCaloClusterHist(cls_hists_[1], cluster);
      if(cluster->Energy() > 70.) FillCaloClusterHist(cls_hists_[2], cluster);
      if(cluster->Energy() > 50.) {
        cut_flow_.Increment("e_50");
        if(cluster->Energy() > 70.) {
          cut_flow_.Increment("e_70");
        }
      }
    }

    return false; // default to not writing output trees, matching the base class
  }

  //------------------------------------------------------------------------------------
  // After the processing loop
  void Run1BAna::EndJob() {
    printf("Run1BAna::%s\n", __func__);
    cut_flow_.Print();
  }
}
