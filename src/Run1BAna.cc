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
    Mu2eEvtAna::InitHistSelections();
    tc_hists_.assign(tc_names_.size(), nullptr);
    for(size_t icoll = 0; icoll < tc_names_.size(); ++icoll) tc_hists_[icoll] = new TimeClusterHist_t;
    ls_hists_.assign(ls_names_.size(), nullptr);
    for(size_t icoll = 0; icoll < ls_names_.size(); ++icoll) ls_hists_[icoll] = new LineSeedHist_t;
  }

  //------------------------------------------------------------------------------------
  // Book the time cluster histograms
  void Run1BAna::BookTimeClusterHist(TimeClusterHist_t* Hist, const char* Folder) {
    if(!Hist) throw std::runtime_error("Attempting to book histograms in a null TimeClusterHist_t\n");
    Hist->fNHits      = new TH1F("nhits"     , Form("%s: N(combo hits)"    , Folder), 100,    0.,  200.);
    Hist->fNStrawHits = new TH1F("nstrawhits", Form("%s: N(straw hits)"    , Folder), 100,    0.,  200.);
    Hist->fT0         = new TH1F("t0"        , Form("%s: T0 (ns)"          , Folder), 200,    0., 2000.);
    Hist->fX          = new TH1F("x"         , Form("%s: x (mm)"           , Folder), 100,-1000., 1000.);
    Hist->fY          = new TH1F("y"         , Form("%s: y (mm)"           , Folder), 100,-1000., 1000.);
    Hist->fZ          = new TH1F("z"         , Form("%s: z (mm)"           , Folder), 100,-5000., 5000.);
    Hist->fR          = new TH1F("r"         , Form("%s: r (mm)"           , Folder), 100,    0., 1000.);
    Hist->fECalo      = new TH1F("ecalo"     , Form("%s: E(calo) (MeV)"    , Folder), 100,    0.,  200.);
    Hist->fTCalo      = new TH1F("tcalo"     , Form("%s: T(calo) (ns)"     , Folder), 200,    0., 2000.);
    Hist->fHasCalo    = new TH1F("hascalo"   , Form("%s: Has calo cluster?", Folder),   2,    0.,    2.);
  }

  //------------------------------------------------------------------------------------
  // Book the line seed histograms
  void Run1BAna::BookLineSeedHist(LineSeedHist_t* Hist, const char* Folder) {
    if(!Hist) throw std::runtime_error("Attempting to book histograms in a null LineSeedHist_t\n");
    Hist->fStatus         = new TH1F("status"        , Form("%s: TrkFitFlag status"  , Folder),  10,   -2.,    8.);
    Hist->fNHits          = new TH1F("nhits"         , Form("%s: N(combo hits)"      , Folder), 100,    0.,  200.);
    Hist->fNStrawHits     = new TH1F("nstrawhits"    , Form("%s: N(straw hits)"      , Folder), 100,    0.,  200.);
    Hist->fT0             = new TH1F("t0"            , Form("%s: T0 (ns)"            , Folder), 200,    0., 2000.);
    Hist->fD0             = new TH1F("d0"            , Form("%s: d_{0} (mm)"         , Folder), 100,-1000., 1000.);
    Hist->fPhi0           = new TH1F("phi0"          , Form("%s: #phi_{0}"           , Folder), 100,   -4.,    4.);
    Hist->fZ0             = new TH1F("z0"            , Form("%s: z_{0} (mm)"         , Folder), 100,-5000., 5000.);
    Hist->fCos            = new TH1F("cos"           , Form("%s: cos(#theta)"        , Folder), 100,   -1.,    1.);
    Hist->fA0             = new TH1F("A0"            , Form("%s: Fit parameter A0"   , Folder), 100,-2000., 2000.);
    Hist->fB0             = new TH1F("B0"            , Form("%s: Fit parameter B0"   , Folder), 100,-8000., 8000.);
    Hist->fA1             = new TH1F("A1"            , Form("%s: Fit parameter A1"   , Folder), 100,  -10.,   10.);
    Hist->fB1             = new TH1F("B1"            , Form("%s: Fit parameter B1"   , Folder), 100,  -60.,   60.);
    Hist->fECalo          = new TH1F("ecalo"         , Form("%s: E(calo) (MeV)"      , Folder), 100,    0.,  200.);
    Hist->fTCalo          = new TH1F("tcalo"         , Form("%s: T(calo) (ns)"       , Folder), 200,    0., 2000.);
    Hist->fHasCalo        = new TH1F("hascalo"       , Form("%s: Has calo cluster?"  , Folder),   2,    0.,    2.);
    Hist->fMatchedTrackDt = new TH1F("matchedtrackdt", Form("%s: T0(seed) - T0(track) (ns)", Folder), 200, -200., 200.);
    Hist->fMatchedTrackDD0= new TH1F("matchedtrackdd0",Form("%s: d0(seed) - d0(track) (mm)", Folder), 200, -200., 200.);
    Hist->fMatchedTCDt    = new TH1F("matchedtcdt"   , Form("%s: T0(seed) - T0(time cluster) (ns)", Folder), 200, -200., 200.);
  }

  //------------------------------------------------------------------------------------
  // Book the histogram sets, extending the base class' event/track/calo/CRV histograms
  void Run1BAna::BookHistograms(TDirectory* dir) {
    Mu2eEvtAna::BookHistograms(dir);

    tc_dirs_.assign(tc_hists_.size(), nullptr);
    for(size_t icoll = 0; icoll < tc_hists_.size(); ++icoll) {
      const TString folder = "tc_" + tc_names_[icoll];
      auto subdir = dir->mkdir(folder);
      subdir->cd();
      BookTimeClusterHist(tc_hists_[icoll], folder.Data());
      dir->cd();
      tc_dirs_[icoll] = subdir;
    }

    ls_dirs_.assign(ls_hists_.size(), nullptr);
    for(size_t icoll = 0; icoll < ls_hists_.size(); ++icoll) {
      const TString folder = "ls_" + ls_names_[icoll];
      auto subdir = dir->mkdir(folder);
      subdir->cd();
      BookLineSeedHist(ls_hists_[icoll], folder.Data());
      dir->cd();
      ls_dirs_[icoll] = subdir;
    }
  }

  //------------------------------------------------------------------------------------
  // Fill the time cluster histograms
  void Run1BAna::FillTimeClusterHist(TimeClusterHist_t* Hist, const TimeCluster_t* Cluster) {
    if(!Hist) {
      if(verbose_ > 0) printf("Run1BAna::%s: Filling time cluster histogram set with null hist par\n", __func__);
      return;
    }
    if(!Cluster || !Cluster->cluster_) {
      if(verbose_ > 0) printf("Run1BAna::%s: Filling time cluster histogram set with null time cluster\n", __func__);
      return;
    }
    Hist->fNHits     ->Fill(Cluster->NHits());
    Hist->fNStrawHits->Fill(Cluster->NStrawHits());
    Hist->fT0        ->Fill(Cluster->T0());
    Hist->fX         ->Fill(Cluster->X());
    Hist->fY         ->Fill(Cluster->Y());
    Hist->fZ         ->Fill(Cluster->Z());
    Hist->fR         ->Fill(Cluster->R());
    Hist->fECalo     ->Fill(Cluster->ECalo());
    if(Cluster->HasCalo()) Hist->fTCalo->Fill(Cluster->TCalo());
    Hist->fHasCalo   ->Fill(Cluster->HasCalo());
  }

  //------------------------------------------------------------------------------------
  // Fill the line seed histograms
  void Run1BAna::FillLineSeedHist(LineSeedHist_t* Hist, const LineSeed_t* Seed) {
    if(!Hist) {
      if(verbose_ > 0) printf("Run1BAna::%s: Filling line seed histogram set with null hist par\n", __func__);
      return;
    }
    if(!Seed || !Seed->seed_) {
      if(verbose_ > 0) printf("Run1BAna::%s: Filling line seed histogram set with null line seed\n", __func__);
      return;
    }
    Hist->fStatus    ->Fill(Seed->Status());
    Hist->fNHits     ->Fill(Seed->NHits());
    Hist->fNStrawHits->Fill(Seed->NStrawHits());
    Hist->fT0        ->Fill(Seed->T0());
    Hist->fD0        ->Fill(Seed->D0());
    Hist->fPhi0      ->Fill(Seed->Phi0());
    Hist->fZ0        ->Fill(Seed->Z0());
    Hist->fCos       ->Fill(Seed->Cos());
    Hist->fA0        ->Fill(Seed->A0());
    Hist->fB0        ->Fill(Seed->B0());
    Hist->fA1        ->Fill(Seed->A1());
    Hist->fB1        ->Fill(Seed->B1());
    Hist->fECalo     ->Fill(Seed->ECalo());
    if(Seed->HasCalo()) Hist->fTCalo->Fill(Seed->TCalo());
    Hist->fHasCalo   ->Fill(Seed->HasCalo());
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

    // Track selection, using the straight-line-appropriate TrackID()
    for(int itrk = 0; itrk < evt_.ntracks_; ++itrk) {
      FillTrackHist(trk_hists_[0], &tracks_[itrk]);
      if(!tracks_[itrk].IsGood()) continue;
      cut_flow_.Increment("a_track");
      const auto ID = TrackID(&tracks_[itrk]);
      if(ID.Passes()) {
        cut_flow_.Increment("track_id");
        FillTrackHist(trk_hists_[1], &tracks_[itrk]);
      }
    }

    // Calo cluster histograms (base-class hist set, never filled by the base since
    // Mu2eEvtAna::ProcessEvent() isn't called here; matched line/line seed/time cluster/CRV
    // cluster info is already populated on each cluster, from InitializeEvent()'s MatchCaloClusters())
    for(int icls = 0; icls < evt_.ncalo_clusters_; ++icls) {
      FillCaloClusterHist(cls_hists_[0], &calo_clusters_[icls]);
    }

    // Time cluster / line seed histograms, and matching for efficiency studies
    bool has_good_seed(false);
    for(size_t icoll = 0; icoll < tc_names_.size(); ++icoll) {
      for(const auto& cluster : time_clusters_[icoll]) FillTimeClusterHist(tc_hists_[icoll], &cluster);
    }
    for(size_t icoll = 0; icoll < ls_names_.size(); ++icoll) {
      for(const auto& seed : line_seeds_[icoll]) {
        FillLineSeedHist(ls_hists_[icoll], &seed);
        if(seed.IsGood()) has_good_seed = true;

        // Match to the best time cluster in the same collection index, if one exists
        const int itc = BestTimeCluster(&seed, icoll);
        if(itc >= 0) ls_hists_[icoll]->fMatchedTCDt->Fill(seed.T0() - time_clusters_[icoll][itc].T0());

        // Match to the best reconstructed track
        int best_track(-1);
        float best_dt(1.e10);
        for(int itrk = 0; itrk < evt_.ntracks_; ++itrk) {
          if(!tracks_[itrk].IsGood()) continue;
          const float dt = std::fabs(seed.T0() - tracks_[itrk].TFront());
          if(dt < best_dt) { best_dt = dt; best_track = itrk; }
        }
        if(best_track >= 0) {
          ls_hists_[icoll]->fMatchedTrackDt ->Fill(seed.T0() - tracks_[best_track].TFront());
          ls_hists_[icoll]->fMatchedTrackDD0->Fill(seed.D0() - tracks_[best_track].D0Front());
        }
      }
    }

    if(has_good_seed) cut_flow_.Increment("has_good_line_seed");

    return false; // default to not writing output trees, matching the base class
  }

  //------------------------------------------------------------------------------------
  // After the processing loop
  void Run1BAna::EndJob() {
    printf("Run1BAna::%s\n", __func__);
    cut_flow_.Print();
  }
}
