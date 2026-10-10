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
    disabled_branches_.erase(
      std::remove(disabled_branches_.begin(), disabled_branches_.end(), TString("calohitsmc")),
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

    // Locate each species' time cluster / line seed collection, and check whether the electron
    // time cluster combo hit lists were stored -- the hit-based cuts below are only evaluable if
    // they were.
    for(int isp = 0; isp < kNSpecies; ++isp) {
      species_tc_[isp] = -1;
      species_ls_[isp] = -1;
      for(size_t icoll = 0; icoll < tc_names_.size(); ++icoll) {
        if(tc_names_[icoll] == kSpeciesTimeClusters[isp]) { species_tc_[isp] = int(icoll); break; }
      }
      for(size_t icoll = 0; icoll < ls_names_.size(); ++icoll) {
        if(ls_names_[icoll] == kSpeciesLineSeeds[isp]) { species_ls_[isp] = int(icoll); break; }
      }
    }
    electron_tc_hits_ = species_tc_[kElectron] >= 0 && event_->HasTimeClusterHits(kSpeciesTimeClusters[kElectron]);

    if(verbose_ > -1) {
      printf("Run1BAna::%s: Found %zu time cluster collection(s):", __func__, tc_names_.size());
      for(const auto& name : tc_names_) printf(" %s%s", name.Data(),
                                               (event_->HasTimeClusterHits(name.Data())) ? " (+hits)" : "");
      printf("\n");
      printf("Run1BAna::%s: Found %zu line seed collection(s):", __func__, ls_names_.size());
      for(const auto& name : ls_names_) printf(" %s%s", name.Data(),
                                               (event_->HasLineSeedHits(name.Data())) ? " (+hits)" : "");
      printf("\n");
    }

    // Warn once, regardless of verbosity: a missing collection leaves that species' matches
    // nullptr, and without the electron hit lists the hit-based part of the set 74 veto is
    // silently skipped, which changes what that set counts
    for(int isp = 0; isp < kNSpecies; ++isp) {
      if(species_tc_[isp] < 0)
        printf("Run1BAna::%s: WARNING: no \"%s\" collection in this ntuple\n", __func__, kSpeciesTimeClusters[isp]);
      if(species_ls_[isp] < 0)
        printf("Run1BAna::%s: WARNING: no \"%s\" collection in this ntuple\n", __func__, kSpeciesLineSeeds[isp]);
    }
    if(species_tc_[kElectron] >= 0 && !electron_tc_hits_) {
      printf("Run1BAna::%s: WARNING: the \"%s\" hit lists were not stored (timeclusters.fillHitsFor)\n",
             __func__, kSpeciesTimeClusters[kElectron]);
    }

    // Physics weights for flat-spectrum samples
    delete spectrum_weight_;
    spectrum_weight_ = nullptr;
    n_no_primary_ = 0;
    if(name_.Contains("fele")) { // flat e- from Al stops, generated flat in momentum over 50-110 MeV/c
      spectrum_weight_ = new SpectrumWeight("Offline/EventGenerator/data/czarnecki_szafron_Al_2016.tbl", 0.51099895);
      spectrum_pdg_  = 11;
      spectrum_pmin_ = 50.;
      spectrum_pmax_ = 110.;
      // The stored events are filtered (calo energy, reco), so their generated momenta are not flat: each event needs
      // its own primary momentum, from the primary branch (EventNtuple versions without it, e.g. Run1Baw, can't be weighted)
      if(!ntuple_->FindBranch("primary"))
        throw std::runtime_error(Form("Run1BAna::%s: %s needs DIO weights, but the ntuple has no primary branch", __func__, name_.Data()));
      printf("Run1BAna::%s: Weighting %s to the DIO spectrum (%s), generated flat over [%.0f, %.0f] MeV/c\n",
             __func__, name_.Data(), spectrum_weight_->File().Data(), spectrum_pmin_, spectrum_pmax_);
    } else if(name_.Contains("fgam") || name_.Contains("pgam")) {
      printf("Run1BAna::%s: WARNING: %s is generated with a flat photon spectrum and is not weighted: "
             "the RMC phase-space model weights are not available\n", __func__, name_.Data());
    }

    // Pion survival weights for RPC samples generated with pion decay off
    rpc_weight_ = name_.Contains("rpce") || name_.Contains("rpci");
    wt_tree_  = nullptr;
    wt_index_ = -1;
    if(rpc_weight_) {
      if(!ntuple_->FindBranch(kEvtWtBranch))
        throw std::runtime_error(Form("Run1BAna::%s: %s needs the pion survival weight, but the ntuple has no %s branch",
                                      __func__, name_.Data(), kEvtWtBranch));
      ntuple_->SetBranchAddress(kEvtWtBranch, &evtwt_);
      printf("Run1BAna::%s: Weighting %s by the pion survival probability (%s%s)\n", __func__, name_.Data(), kEvtWtBranch, kRPCWeightLeaf);
    }
    return 0;
  }

  //------------------------------------------------------------------------------------
  // Physics weight of the event: spectrum weights for flat-spectrum samples, survival weights for RPC
  void Run1BAna::SetEventWeight() {
    if(spectrum_weight_) {
      const float p = PrimaryMomentum(spectrum_pdg_);
      if(p > 0.f) {
        evt_.weight_ *= spectrum_weight_->Weight(p);
      } else { // no primary with this PDG ID: not expected for these samples, counted and reported in EndJob()
        evt_.weight_ = 0.;
        ++n_no_primary_;
      }
    }

    if(rpc_weight_) {
      // The evtwt leaves follow the input's EventWeight module labels, so find the generator's leaf in each input tree
      TTree* tree = ntuple_->GetTree();
      if(tree != wt_tree_) {
        wt_tree_  = tree;
        wt_index_ = -1;
        TBranch* branch = (tree) ? tree->GetBranch(kEvtWtBranch) : nullptr;
        TObjArray* leaves = (branch) ? branch->GetListOfLeaves() : nullptr;
        for(int ileaf = 1; leaves && ileaf < leaves->GetEntries(); ++ileaf) { // leaf 0 is the weight count
          if(TString(leaves->At(ileaf)->GetName()) == kRPCWeightLeaf) { wt_index_ = ileaf - 1; break; }
        }
        if(wt_index_ < 0 || wt_index_ >= EvtWt_t::kMaxWeights)
          throw std::runtime_error(Form("Run1BAna::%s: No %s%s weight in %s", __func__, kEvtWtBranch, kRPCWeightLeaf,
                                        (ntuple_->GetCurrentFile()) ? ntuple_->GetCurrentFile()->GetName() : "the input"));
      }
      // The weight is exp(-pion proper time), in [0, 1]. Anything else (-1 for an event without weights, inf, or stale
      // values from EventNtuple versions that filled evtwt from the previous event's handles) can't be used.
      const float w = evtwt_.weights[wt_index_];
      if(!std::isfinite(w) || w < 0.f || w > 1.f)
        throw std::runtime_error(Form("Run1BAna::%s: Invalid pion survival weight %g in event %i:%i:%i -- this EventNtuple "
                                      "version may fill evtwt from the previous event (fixed in EventNtupleMaker)",
                                      __func__, w, evt_.run_, evt_.subrun_, evt_.event_));
      evt_.weight_ *= w;
    }
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
    hist_sets[ 40] = new hist_info_t("base_proton"                     ,  true,  true,  true,  true,  true,  true, false, false);
    hist_sets[ 43] = new hist_info_t("proton_edep"                     ,  true,  true,  true,  true,  true,  true, false, false);
    hist_sets[ 44] = new hist_info_t("proton_id"                       ,  true,  true,  true,  true,  true,  true, false, false);
    hist_sets[ 70] = new hist_info_t("60 < E < 120 MeV"                ,  true,  true,  true,  true,  true,  true, false, false);
    hist_sets[ 71] = new hist_info_t("rmc_base"                        ,  true,  true,  true,  true,  true,  true, false, false);
    hist_sets[ 72] = new hist_info_t("rmc_r_cut"                       ,  true,  true,  true,  true,  true,  true, false, false);
    hist_sets[ 73] = new hist_info_t("rmc_line_cut"                    ,  true,  true,  true,  true,  true,  true, false, false);
    hist_sets[ 74] = new hist_info_t("rmc_seed_cut"                    ,  true,  true,  true,  true,  true,  true, false, false);
    hist_sets[ 80] = new hist_info_t("ce_id"                           ,  true,  true,  true,  true,  true,  true, false, false);
    hist_sets[ 90] = new hist_info_t("rpc_base"                        ,  true,  true,  true,  true,  true,  true, false, false);
    hist_sets[ 93] = new hist_info_t("rpc_line_cut"                    ,  true,  true,  true,  true,  true,  true, false, false);
    hist_sets[ 94] = new hist_info_t("rpc_seed_cut"                    ,  true,  true,  true,  true,  true,  true, false, false);

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
        tc.hits_    = cluster.hits; // null unless this collection's hit list was stored
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
        ls.hits_ = seed.hits; // null unless this collection's hit list was stored
        coll.emplace_back(ls);
      }
    }

    // Match every calo cluster to its best line/line seed/time cluster/CRV cluster. Must run
    // last: it needs tracks_, crv_clusters_ (both filled by Mu2eEvtAna::InitializeEvent() above)
    // and time_clusters_/line_seeds_ (just filled above) all populated for this event.
    MatchCaloClusters();

    SetEventWeight();
  }

  //------------------------------------------------------------------------------------
  // Populate each calo cluster's best-matched electron/proton/cosmic line, line seed, and time
  // cluster, and its CRV cluster. Each species' line seed / time cluster is searched for only in
  // that species' collection (species_ls_/species_tc_), and its line only among the tracks fit
  // with that species' hypothesis (kSpeciesFitPDG); a missing collection leaves the match nullptr.
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

      Track_t*       lines        [kNSpecies] = {nullptr, nullptr, nullptr};
      LineSeed_t*    line_seeds   [kNSpecies] = {nullptr, nullptr, nullptr};
      TimeCluster_t* time_clusters[kNSpecies] = {nullptr, nullptr, nullptr};
      for(int isp = 0; isp < kNSpecies; ++isp) {
        // Best matched line (track) with this species' fit hypothesis
        for(int itrk = 0; itrk < evt_.ntracks_ && !lines[isp]; ++itrk) {
          Track_t* trk = &tracks_[itrk];
          if(!trk->IsGood()) continue;
          if(std::abs(trk->FitPDG()) != kSpeciesFitPDG[isp]) continue;
          const auto* tch = trk->TCH();
          if(!tch) continue;
          if(std::fabs(tch->edep - E) < kCaloMatchEps && std::fabs(tch->ctime - T) < kCaloMatchEps) lines[isp] = trk;
        }

        // Best matched line seed from this species' collection
        if(species_ls_[isp] >= 0) {
          for(auto& seed : line_seeds_[species_ls_[isp]]) {
            if(!seed.HasCalo()) continue;
            if(std::fabs(seed.ECalo() - E) < kCaloMatchEps && std::fabs(seed.TCalo() - T) < kCaloMatchEps) { line_seeds[isp] = &seed; break; }
          }
        }

        // Best matched time cluster from this species' collection
        if(species_tc_[isp] >= 0) {
          for(auto& tc : time_clusters_[species_tc_[isp]]) {
            if(!tc.HasCalo()) continue;
            if(std::fabs(tc.ECalo() - E) < kCaloMatchEps && std::fabs(tc.TCalo() - T) < kCaloMatchEps) { time_clusters[isp] = &tc; break; }
          }
        }
      }
      cls.electron_line_         = lines        [kElectron];
      cls.proton_line_           = lines        [kProton  ];
      cls.cosmic_line_           = lines        [kCosmic  ];
      cls.electron_line_seed_    = line_seeds   [kElectron];
      cls.proton_line_seed_      = line_seeds   [kProton  ];
      cls.cosmic_line_seed_      = line_seeds   [kCosmic  ];
      cls.electron_time_cluster_ = time_clusters[kElectron];
      cls.proton_time_cluster_   = time_clusters[kProton  ];
      cls.cosmic_time_cluster_   = time_clusters[kCosmic  ];

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

    if(track->TSDAInter())                                         ID.SetBit(kRMax);

    return ID;
  }

  //------------------------------------------------------------------------------------
  // Main event-by-event processing
  bool Run1BAna::ProcessEvent() {
    cut_flow_ce_     .ResetEvent();
    cut_flow_rmc_    .ResetEvent();
    cut_flow_rpc_    .ResetEvent();
    cut_flow_proton_ .ResetEvent();
    cut_flow_neutron_.ResetEvent();

    // Identify the sample type
    const bool is_pu = name_.Contains("mnbs");


    FillEventHist(evt_hists_[0]); // all events with well defined inputs

    // Every time cluster / line seed of the event, from every collection. Only set 0 ("All
    // events") is event-level; the selections above it are per-calo-cluster, so filling a
    // per-event object into them would not mean the same thing.
    // FIXME: the time cluster / line seed histogram sets booked for the other selections are
    // still never filled -- they need a per-calo-cluster definition first.
    for(size_t icoll = 0; icoll < tc_names_.size(); ++icoll) {
      for(const auto& tc : time_clusters_[icoll]) FillTimeClusterHist(tcs_hists_[0], &tc);
    }
    for(size_t icoll = 0; icoll < ls_names_.size(); ++icoll) {
      for(const auto& seed : line_seeds_[icoll]) FillLineSeedHist(lns_hists_[0], &seed);
    }

    for(int icls = 0; icls < evt_.ncalo_clusters_; ++icls) {
      const auto cluster = &calo_clusters_[icls];
      cut_flow_ce_     .Increment("has_cluster");
      cut_flow_rmc_    .Increment("has_cluster");
      cut_flow_rpc_    .Increment("has_cluster");
      cut_flow_proton_ .Increment("has_cluster");
      cut_flow_neutron_.Increment("has_cluster");

      //------------------------------------------------
      // Common variables
      //------------------------------------------------

      const float energy        = cluster->Energy();
      const float time          = cluster->Time();
      const int   ncr           = cluster->NCrystals();
      const float e1_r          = cluster->E1() / energy;
      const float e2_r          = cluster->E2() / energy;
      const float e24_r         = (cluster->E25() - cluster->E1()) / energy;
      const float e2p_r         = (cluster->E2()-cluster->E1()) / energy;
      const float tvar          = cluster->TVar();
      const float second_moment = cluster->SecondMoment();
      const int   disk          = cluster->DiskID();
      const float r             = cluster->R();
      const auto e_tc           = cluster->electron_time_cluster_;
      const auto e_line_seed    = cluster->electron_line_seed_;
      const auto e_line         = cluster->electron_line_;
      const auto p_tc           = cluster->proton_time_cluster_;
      const auto p_line_seed    = cluster->proton_line_seed_;
      const auto p_line         = cluster->proton_line_;
      const auto c_line         = cluster->cosmic_line_;
      const int   n_e_hits_high = (e_tc) ? e_tc->NHitsAboveZ(1300.) : -1;

      const int sim_pdg = cluster->MCPDG();
      const float sim_edep = cluster->MCSimEDep();
      bool mc_veto = false;
      if(is_pu) {
        // skip high energy tail RMC/protons/neutrons
        mc_veto |= (sim_pdg == 22 || sim_pdg == 2212 || sim_pdg == 2112) && sim_edep > 60.f;
      }

      if(!mc_veto) {
        cut_flow_ce_     .Increment("mc_veto");
        cut_flow_rmc_    .Increment("mc_veto");
        cut_flow_rpc_    .Increment("mc_veto");
        cut_flow_proton_ .Increment("mc_veto");
        cut_flow_neutron_.Increment("mc_veto");
      }

      // Base selection cuts
      const float min_energy         =   60.f;
      const float max_energy         =  120.f;
      const float min_time           =  500.f;
      const float max_time           = 1650.f;
      const int   min_cr             =      2;
      const int   max_cr             =      5;
      const float min_e1_r           =   0.6f;
      const float min_e2_r           =   0.8f;
      const float max_tvar           =    1.f;
      const float max_sec_mom        =  1.e3f;
      const float min_r              =  500.f;
      const float max_r              =  580.f;
      const int   max_e_tc_hits_high = 2;
      const bool  base_id            = (!mc_veto &&
                                        energy > min_energy && energy < max_energy
                                        && time > min_time && time < max_time);
      const bool  pu_veto            = (ncr >= min_cr && ncr <= max_cr &&
                                        e1_r > min_e1_r && e2_r > min_e2_r &&
                                        tvar < max_tvar &&
                                        second_moment < max_sec_mom &&
                                        disk == 0);
      const bool  pu_r_veto          = r > min_r && r < max_r;

      //------------------------------------------------
      // Basic selections
      //------------------------------------------------

      FillCaloClusterHist(cls_hists_[0], cluster);
      if(energy > 50.) FillCaloClusterHist(cls_hists_[1], cluster);
      if(energy > 70.) FillCaloClusterHist(cls_hists_[2], cluster);
      if(base_id) {
        FillCaloClusterHist(cls_hists_[70], cluster);
        FillTimeClusterHist(tcs_hists_[70], e_tc);
        if(pu_veto) {
          FillCaloClusterHist(cls_hists_[71], cluster);
          FillTimeClusterHist(tcs_hists_[71], e_tc);
          if(pu_r_veto) {
            FillCaloClusterHist(cls_hists_[72], cluster);
            FillTimeClusterHist(tcs_hists_[72], e_tc);
          }
        }
      }

      //------------------------------------------------
      // RMC selection
      //------------------------------------------------

      if(base_id && pu_veto && pu_r_veto) {
        bool trk_veto = false;
        trk_veto |= e_line != nullptr;
        if(!trk_veto) {
          FillCaloClusterHist(cls_hists_[73], cluster);
          FillTimeClusterHist(tcs_hists_[73], e_tc);
          trk_veto |= e_line_seed != nullptr;
          if(electron_tc_hits_ && e_tc && e_tc->HasHits()) {
            trk_veto |= e_tc->NHitsAboveZ(1300.) >= 3;
          }
          if(!trk_veto) {
            if(is_pu && energy > 70.f) PrintClusterInfo("[Accepted RMC: PU]", cluster);
            FillCaloClusterHist(cls_hists_[74], cluster);
            FillTimeClusterHist(tcs_hists_[74], e_tc);
          }
        }
      }

      //------------------------------------------------
      // RPC selection
      //------------------------------------------------

      const float rpc_energy_min =  60.;
      const float rpc_energy_max = 140.;
      const float rpc_time_min   = 300.;
      const float rpc_time_max   = 550.;

      // cut-flow
      // Cut Group  Events Passing  Absolute [%]  Relative [%]                                Description
      //                             No cuts   N/A         2351533    100.000000    100.000000                       No selection applied
      //                  cluster_energy_cut   NaN          791591     33.662764     33.662764                           60 < E < 150 MeV
      //                    cluster_time_cut   NaN           67572      2.873530      8.536226                           300 < T < 550 ns
      //                    num_crystals_cut   NaN           34187      1.453818     50.593441                        1 < N(crystals) < 6
      //                         disk_id_cut   NaN           17381      0.739135     50.840963                                  Disk ID 0
      //                  cluster_radius_cut   NaN            3596      0.152922     20.689258           500 mm < Cluster radius < 580 mm
      //           cluster_second_moment_cut   NaN             883      0.037550     24.555061    Cluster energy second moment < 1000 MeV
      //     max_crystal_energy_fraction_cut   NaN             135      0.005741     15.288788 Ei / E > 0.6 (max crystal energy fraction)
      // top_two_crystal_energy_fraction_cut   NaN             130      0.005528     96.296296                            (E1+E2)/E > 0.8
      //        high_z_time_cluster_hits_cut   NaN              64      0.002722     49.230769            N(high z time cluster hits) < 3
      //            e5x5_energy_fraction_cut   NaN              58      0.002466     90.625000                         (E5x5-E1)/E > 7.5%
      //  second_crystal_energy_fraction_cut   NaN              31      0.001318     53.448276                                 E2/E < 30%
      if(energy > rpc_energy_min && energy < rpc_energy_max) {
        cut_flow_rpc_.Increment("cluster_energy");
        if(time > rpc_time_min && time < rpc_time_max) {
          cut_flow_rpc_.Increment("cluster_time");
          if(ncr >= min_cr && ncr <= max_cr) {
            cut_flow_rpc_.Increment("num_crystals");
            if(disk == 0) {
              cut_flow_rpc_.Increment("disk_id");
              if(second_moment < max_sec_mom) {
                cut_flow_rpc_.Increment("second_moment");
                if(e1_r > min_e1_r) {
                  cut_flow_rpc_.Increment("frac_e1");
                  if(e2_r > min_e2_r) {
                    cut_flow_rpc_.Increment("frac_e1+2");
                    if(n_e_hits_high <= max_e_tc_hits_high) {
                      cut_flow_rpc_.Increment("high_z_hits");
                      if(e24_r > 0.075) {
                        cut_flow_rpc_.Increment("frac_e24");
                        if(e2p_r < 0.3) {
                          cut_flow_rpc_.Increment("frac_e2p");
                          if(e_line == nullptr) {
                            cut_flow_rpc_.Increment("no_ele_line");
                            if(c_line == nullptr) {
                              cut_flow_rpc_.Increment("no_cosmic_line");
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }

      // fill histograms
      if(!mc_veto &&
         energy > rpc_energy_min && energy < rpc_energy_max &&
         time > rpc_time_min && time < rpc_time_max &&
         pu_veto && pu_r_veto) {
        FillCaloClusterHist(cls_hists_[90], cluster);
        FillTimeClusterHist(tcs_hists_[90], e_tc);
        bool trk_veto = false;
        trk_veto |= e_line != nullptr;
        if(!trk_veto) {
          FillCaloClusterHist(cls_hists_[93], cluster);
          FillTimeClusterHist(tcs_hists_[93], e_tc);
          trk_veto |= e_line_seed != nullptr;
          if(electron_tc_hits_ && e_tc && e_tc->HasHits()) {
            trk_veto |= e_tc->NHitsAboveZ(1300.) >= 3;
          }
          if(!trk_veto) {
            if(is_pu && energy > 70.f) PrintClusterInfo("[Accepted RPC: PU]", cluster);
            FillCaloClusterHist(cls_hists_[94], cluster);
            FillTimeClusterHist(tcs_hists_[94], e_tc);
          }
        }
      }

      //------------------------------------------------
      // Electron selection
      //------------------------------------------------

      if(base_id && pu_veto && pu_r_veto && e_tc && e_line) {
        const int tc_nhits = e_tc->NHits();
        bool ce_id = tc_nhits > 10 && tc_nhits < 40;
        ce_id &= e_line->NActive() > 0;
        ce_id &= std::fabs(e_line->CosThetaFront()) > 0.985;
        if(ce_id) {
          FillCaloClusterHist(cls_hists_[80], cluster);
          FillTimeClusterHist(tcs_hists_[80], e_tc);
          FillLineSeedHist   (lns_hists_[80], e_line_seed);
          FillTrackHist      (trk_hists_[80], e_line);
        }
      }

      //------------------------------------------------
      // Proton selection
      //------------------------------------------------

      if(base_id && p_tc && p_line) {
        FillCaloClusterHist(cls_hists_[40], cluster);
        FillTimeClusterHist(tcs_hists_[40], p_tc);
        FillLineSeedHist   (lns_hists_[40], p_line_seed);
        FillTrackHist      (trk_hists_[40], p_line);
        if(p_tc->AvgEDep() > 0.0028) { // FIXME: Use line hits
          FillCaloClusterHist(cls_hists_[43], cluster);
          FillTimeClusterHist(tcs_hists_[43], p_tc);
          FillLineSeedHist   (lns_hists_[43], p_line_seed);
          FillTrackHist      (trk_hists_[43], p_line);
          if(disk == 0
             && ncr < 4
             && std::fabs(p_line->CosThetaFront()) > 0.985) {
            FillCaloClusterHist(cls_hists_[44], cluster);
            FillTimeClusterHist(tcs_hists_[44], p_tc);
            FillLineSeedHist   (lns_hists_[44], p_line_seed);
            FillTrackHist      (trk_hists_[44], p_line);
          }
        }
      }

      //------------------------------------------------
      // Neutron selection
      //------------------------------------------------

      //------------------------------------------------
      // RPC selection
      //------------------------------------------------

    } // end of cluster loop

    return false; // default to not writing output trees, matching the base class
  }

  //------------------------------------------------------------------------------------
  // After the processing loop
  void Run1BAna::EndJob() {
    printf("Run1BAna::%s\n", __func__);
    if(spectrum_weight_ && n_no_primary_ > 0)
      printf("  WARNING: %lld events had no primary with PDG ID %i and were given zero weight\n", n_no_primary_, spectrum_pdg_);
    printf("  CE cut-flow:\n");
    cut_flow_ce_.Print();
    printf("  RMC cut-flow:\n");
    cut_flow_rmc_.Print();
    printf("  RPC cut-flow:\n");
    cut_flow_rpc_.Print();
    printf("  Proton cut-flow:\n");
    cut_flow_proton_.Print();
    printf("  Neutron cut-flow:\n");
    cut_flow_neutron_.Print();
  }
}
