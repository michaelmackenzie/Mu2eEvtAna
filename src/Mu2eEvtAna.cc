#include "Mu2eEvtAna/inc/Mu2eEvtAna.hh"

using namespace mu2e;
namespace Mu2eEvtAna {

  //------------------------------------------------------------------------------------
  // Constructor
  Mu2eEvtAna::Mu2eEvtAna(int verbose) : ntuple_(nullptr), evt_info_(nullptr), evt_info_mc_(nullptr), trk_(nullptr),
                                        trk_calo_hit_(nullptr), trk_qual_(nullptr), hit_cnt_(nullptr), trk_cnt_(nullptr),
                                        fout_(nullptr), tout_(nullptr), tnorm_(nullptr), name_("test"),
                                        report_rate_(1000), verbose_(verbose) {
    //initialize the arrays to null
    for(int ihist = 0; ihist < kMaxEventHists; ++ihist) evt_hists_[ihist] = nullptr;

    //Default to one histogram set/selection
    evt_hists_[0] = new EventHist_t;

    // Track collections to keep
    tracks_  .keep_ =  true;
    electron_.keep_ =  true;
    muon_    .keep_ =  true;
    proton_  .keep_ = false;
  }

  //------------------------------------------------------------------------------------
  // Initialize the input ntuple information
  int Mu2eEvtAna::InitializeInput() {
    if(verbose_ > 1) printf("Mu2eEvtAna::%s: Adding input information\n", __func__);
    if(!ntuple_) {
      if(verbose_ > -2) printf("Mu2eEvtAna::%s: No ntuple is defined to configure!\n", __func__);
      return -1;
    }

    //-----------------------------
    // Initialize event info branches

    if(ntuple_->GetBranch("evtinfo")) ntuple_->SetBranchAddress("evtinfo", &evt_info_);
    else if(verbose_ > -1) {
      printf("Mu2eEvtAna::%s: Warning! No event info branch is found in the ntuple!\n", __func__);
    }
    // (likely) only defined for simulations
    if(ntuple_->GetBranch("evtinfomc")) ntuple_->SetBranchAddress("evtinfomc", &evt_info_mc_);

    //-----------------------------
    // Initialize count branches

    ntuple_->SetBranchAddress("hitcount.", &hit_cnt_);
    if(ntuple_->GetBranch("trkcount.")) ntuple_->SetBranchAddress("trkcount.", &trk_cnt_);

    //-----------------------------
    // Initialize track branches

    ntuple_->SetBranchAddress("trk", &trk_);
    ntuple_->SetBranchAddress("trkcalohit", &trk_calo_hit_);
    ntuple_->SetBranchAddress("trkqual", &trk_qual_);


    if(verbose_ > 1) printf("Mu2eEvtAna::%s: Initialized input data\n", __func__);
    return 0;
  }

  //------------------------------------------------------------------------------------
  // Initialize the histograms for an event selection
  void Mu2eEvtAna::BookEventHist(EventHist_t* Hist) {
    if(!Hist) {
      throw std::runtime_error("Attempting to book histograms in a null EventHist_t\n");
    }
    Hist->hNTrks     = new TH1D("ntrks"    , "N(tracks)"  ,  10, 0,  10);
    Hist->hNDigis    = new TH1D("ndigis"   , "N(digis)"   , 200, 0, 200);
    Hist->hNClusters = new TH1D("nclusters", "N(clusters)",  10, 0,  10);
  }

  //------------------------------------------------------------------------------------
  // Initialize the histogram sets
  void Mu2eEvtAna::BookHistograms(TDirectory* dir) {

    for(int ihist = 0; ihist < kMaxEventHists; ++ihist) {
      if(evt_hists_[ihist]) {
        auto subdir = dir->mkdir(Form("hist_%i", ihist));
        subdir->cd();
        BookEventHist(evt_hists_[ihist]);
        dir->cd();
        evt_dirs_[ihist] = subdir;
      }
    }

  }

  //------------------------------------------------------------------------------------
  // Fill the histograms for an event selection
  void Mu2eEvtAna::FillEventHist(EventHist_t* Hist) {
    Hist->hNTrks     ->Fill(evt_.ntrks_, evt_.weight_);
    Hist->hNDigis    ->Fill(evt_.ndigis_, evt_.weight_);
    Hist->hNClusters ->Fill(evt_.nclusters_, evt_.weight_);
  }

  //------------------------------------------------------------------------------------
  // Add an output track collection
  void Mu2eEvtAna::AddTrackCollection(TTree* t, Tracks_t& trk, const char* name) {
    t->Branch(Form("n%s", name), &trk.size_);
    t->Branch(Form("%s_status"         , name), trk.status_         , Form("%s_status[n%s]/I"         , name, name));
    t->Branch(Form("%s_pdg"            , name), trk.pdg_            , Form("%s_pdg[n%s]/I"            , name, name));
    t->Branch(Form("%s_nhits"          , name), trk.nhits_          , Form("%s_nhits[n%s]/I"          , name, name));
    t->Branch(Form("%s_nactive"        , name), trk.nactive_        , Form("%s_nactive[n%s]/I"        , name, name));
    t->Branch(Form("%s_ndouble"        , name), trk.ndouble_        , Form("%s_ndouble[n%s]/I"        , name, name));
    t->Branch(Form("%s_nplanes"        , name), trk.nplanes_        , Form("%s_nplanes[n%s]/I"        , name, name));
    t->Branch(Form("%s_nnull"          , name), trk.nnull_          , Form("%s_nnull[n%s]/I"          , name, name));
    t->Branch(Form("%s_first_hit"      , name), trk.first_hit_      , Form("%s_first_hit[n%s]/I"      , name, name));
    t->Branch(Form("%s_last_hit"       , name), trk.last_hit_       , Form("%s_last_hit[n%s]/I"       , name, name));
    t->Branch(Form("%s_ndof"           , name), trk.ndof_           , Form("%s_ndof[n%s]/I"           , name, name));
    t->Branch(Form("%s_chisq"          , name), trk.chisq_          , Form("%s_chisq[n%s]/F"          , name, name));
    t->Branch(Form("%s_fitcon"         , name), trk.fitcon_         , Form("%s_fitcon[n%s]/F"         , name, name));
    t->Branch(Form("%s_reco_id"        , name), trk.reco_id_        , Form("%s_reco_id[n%s]/I"        , name, name));
    t->Branch(Form("%s_trk_qual"       , name), trk.trk_qual_       , Form("%s_trk_qual[n%s]/F"       , name, name));
    t->Branch(Form("%s_p_corr"         , name), trk.p_corr_         , Form("%s_p_corr[n%s]/F"         , name, name));
    t->Branch(Form("%s_id_weight"      , name), trk.id_weight_      , Form("%s_id_weight[n%s]/F"      , name, name));
    t->Branch(Form("%s_gen_trk_index"  , name), trk.gen_trk_index_  , Form("%s_gen_trk_index[n%s]/I"  , name, name));
    t->Branch(Form("%s_seg_first_index", name), trk.seg_first_index_, Form("%s_seg_first_index[n%s]/I", name, name));
    t->Branch(Form("%s_trk_calo_index" , name), trk.trk_calo_index_ , Form("%s_trk_calo_index[n%s]/I" , name, name));

    //Track+Calo hit info for the tracks
    t->Branch(Form("%s_ncalohit", name), &trk.calo_hits_.size_);
    t->Branch(Form("%s_calohit_index"    , name), trk.calo_hits_.index_     , Form("%s_calohit_index[%s_ncalohit]/I"    , name, name));
    t->Branch(Form("%s_calohit_disk_id"  , name), trk.calo_hits_.disk_id_   , Form("%s_calohit_disk_id[%s_ncalohit]/I"  , name, name));
    t->Branch(Form("%s_calohit_cl_t"     , name), trk.calo_hits_.cl_t_      , Form("%s_calohit_cl_t[%s_ncalohit]/F"     , name, name));
    t->Branch(Form("%s_calohit_cl_t_err" , name), trk.calo_hits_.cl_t_err_  , Form("%s_calohit_cl_t_err[%s_ncalohit]/F" , name, name));
    t->Branch(Form("%s_calohit_cl_e"     , name), trk.calo_hits_.cl_e_      , Form("%s_calohit_cl_e[%s_ncalohit]/F"     , name, name));
    t->Branch(Form("%s_calohit_cl_e_err" , name), trk.calo_hits_.cl_e_err_  , Form("%s_calohit_cl_e_err[%s_ncalohit]/F" , name, name));
    t->Branch(Form("%s_calohit_cl_size"  , name), trk.calo_hits_.cl_size_   , Form("%s_calohit_cl_size[%s_ncalohit]/F"  , name, name));
    t->Branch(Form("%s_calohit_cl_ncr"   , name), trk.calo_hits_.cl_ncr_    , Form("%s_calohit_cl_ncr[%s_ncalohit]/I"   , name, name));
    t->Branch(Form("%s_calohit_x"        , name), trk.calo_hits_.x_         , Form("%s_calohit_x[%s_ncalohit]/F"        , name, name));
    t->Branch(Form("%s_calohit_y"        , name), trk.calo_hits_.y_         , Form("%s_calohit_y[%s_ncalohit]/F"        , name, name));
    t->Branch(Form("%s_calohit_z"        , name), trk.calo_hits_.z_         , Form("%s_calohit_z[%s_ncalohit]/F"        , name, name));
    t->Branch(Form("%s_calohit_px"       , name), trk.calo_hits_.px_        , Form("%s_calohit_px[%s_ncalohit]/F"       , name, name));
    t->Branch(Form("%s_calohit_py"       , name), trk.calo_hits_.py_        , Form("%s_calohit_py[%s_ncalohit]/F"       , name, name));
    t->Branch(Form("%s_calohit_pz"       , name), trk.calo_hits_.pz_        , Form("%s_calohit_pz[%s_ncalohit]/F"       , name, name));
    t->Branch(Form("%s_calohit_doca"     , name), trk.calo_hits_.doca_      , Form("%s_calohit_doca[%s_ncalohit]/F"     , name, name));
    t->Branch(Form("%s_calohit_dt"       , name), trk.calo_hits_.dt_        , Form("%s_calohit_dt[%s_ncalohit]/F"       , name, name));
  }

  //------------------------------------------------------------------------------------
  // Setup the output ntuple structure
  void Mu2eEvtAna::AddOutputBranches(TTree* t) {

    // Event branches
    t->Branch("ntrks", &evt_.ntrks_);
    t->Branch("ndigis", &evt_.ndigis_);
    t->Branch("nclusters", &evt_.nclusters_);
    t->Branch("nelectrons", &evt_.nelectrons_);
    t->Branch("npositrons", &evt_.npositrons_);
    t->Branch("nmuminus", &evt_.nmu_minus_);
    t->Branch("nmuplus", &evt_.nmu_plus_);

    // Track branches
    if(tracks_  .keep_) AddTrackCollection(t, tracks_, "track");
    if(electron_.keep_) AddTrackCollection(t, electron_, "electron");
    if(muon_    .keep_) AddTrackCollection(t, muon_, "muon");
    if(proton_  .keep_) AddTrackCollection(t, proton_, "proton");

    // Cluster branches

    // Trk+calo hit branches

    // CRV branches

  }

  //------------------------------------------------------------------------------------
  // Initialize the input ntuple information
  int Mu2eEvtAna::InitializeOutput() {
    fout_ = new TFile(OutputFileName(), "RECREATE");
    top_dir_ = fout_->mkdir("Data");
    top_dir_->cd();

    // Initialize the output ntuple
    tout_ = new TTree("evtana", "Slim Mu2e event analysis tree");
    AddOutputBranches(tout_);

    // Initialize the normalization tree
    tnorm_ = new TTree("Norm", "Normalization information");
    tnorm_->Branch("ngen"   , &norm_.ngen_   );
    tnorm_->Branch("nntuple", &norm_.nntuple_);
    tnorm_->Branch("nseen"  , &norm_.nseen_  );
    tnorm_->Branch("naccept", &norm_.naccept_);
    tnorm_->Branch("nneg"   , &norm_.nneg_   );

    // Initialize histograms
    BookHistograms(top_dir_);

    if(verbose_ > 1) printf("Mu2eEvtAna::%s: Created output file and trees/histograms\n", __func__);
    return 0;
  }

  //------------------------------------------------------------------------------------
  // Process output for an accepted event
  void Mu2eEvtAna::FillOutput() {
    tout_->Fill();
  }

  //------------------------------------------------------------------------------------
  // Add a track to a track collection
  void Mu2eEvtAna::AddTrack(Tracks_t& trks, const Int_t index) {
    if(trks.size_ >= kMaxTracks) {
      throw std::runtime_error(Form("Exceeded maximum allowed tracks! N(tracks) = %i", trks.size_));
    }
    if(trks.size_ < 0) {
      printf("Mu2eEvtAna::%s: Entry %lld: Warning! Negative track count (%i), resetting to 0\n",
             __func__, entry_, trks.size_);
      trks.size_ = 0;
    }
    if(index < 0) return;

    const auto& track = (*trk_)[index];

    // initialize the track information
    trks.index_          [trks.size_] = index;
    trks.status_         [trks.size_] = track.status;
    trks.pdg_            [trks.size_] = track.pdg;
    trks.nhits_          [trks.size_] = track.nhits;
    trks.nactive_        [trks.size_] = track.nactive;
    trks.ndouble_        [trks.size_] = track.ndouble;
    trks.nplanes_        [trks.size_] = track.nplanes;
    trks.nnull_          [trks.size_] = track.nnullambig;
    trks.first_hit_      [trks.size_] = track.firsthit;
    trks.last_hit_       [trks.size_] = track.lasthit;
    trks.ndof_           [trks.size_] = track.ndof;
    trks.chisq_          [trks.size_] = track.chisq;
    trks.fitcon_         [trks.size_] = track.fitcon;
    trks.reco_id_        [trks.size_] = 0;
    trks.trk_qual_       [trks.size_] = (trk_qual_ && (*trk_qual_)[index].valid) ? (*trk_qual_)[index].result : -2.f;
    trks.p_corr_         [trks.size_] = 0.f;
    trks.id_weight_      [trks.size_] = 1.f;
    trks.gen_trk_index_  [trks.size_] = -1;
    trks.seg_first_index_[trks.size_] = -1;

    // add track calo hit information if available
    if(trk_calo_hit_ && index < int(trk_calo_hit_->size()) && (*trk_calo_hit_)[index].active) {
      trks.trk_calo_index_[trks.size_] = trks.calo_hits_.size_;

      const auto& calo_hit = (*trk_calo_hit_)[index];
      trks.calo_hits_.index_    [trks.calo_hits_.size_] = index;
      trks.calo_hits_.disk_id_  [trks.calo_hits_.size_] = calo_hit.did;
      trks.calo_hits_.cl_t_     [trks.calo_hits_.size_] = calo_hit.ctime;
      trks.calo_hits_.cl_t_err_ [trks.calo_hits_.size_] = calo_hit.ctime;
      trks.calo_hits_.cl_e_     [trks.calo_hits_.size_] = calo_hit.edep;
      trks.calo_hits_.cl_e_err_ [trks.calo_hits_.size_] = calo_hit.edeperr;
      trks.calo_hits_.cl_size_  [trks.calo_hits_.size_] = calo_hit.csize;
      trks.calo_hits_.cl_ncr_   [trks.calo_hits_.size_] = 0; //FIXME: Get N(crystals)
      trks.calo_hits_.x_        [trks.calo_hits_.size_] = calo_hit.poca.x();
      trks.calo_hits_.y_        [trks.calo_hits_.size_] = calo_hit.poca.y();
      trks.calo_hits_.z_        [trks.calo_hits_.size_] = calo_hit.poca.z();
      trks.calo_hits_.px_       [trks.calo_hits_.size_] = calo_hit.mom.x();
      trks.calo_hits_.py_       [trks.calo_hits_.size_] = calo_hit.mom.y();
      trks.calo_hits_.pz_       [trks.calo_hits_.size_] = calo_hit.mom.z();
      trks.calo_hits_.doca_     [trks.calo_hits_.size_] = calo_hit.doca;
      trks.calo_hits_.dt_       [trks.calo_hits_.size_] = calo_hit.dt;

      // increment the collection size
      trks.calo_hits_.size_ = trks.calo_hits_.size_+1;
    } else {
      trks.trk_calo_index_[trks.size_] = -1;
    }

    // increment the collection size
    ++trks.size_;
  }

  //------------------------------------------------------------------------------------
  // Initialize event information
  void Mu2eEvtAna::InitializeEvent() {
    evt_.ntrks_ = (trk_) ? trk_->size() : 0;
    evt_.ndigis_ = (hit_cnt_) ? hit_cnt_->nsd : 0;
    evt_.nelectrons_ = 0;
    evt_.npositrons_ = 0;
    evt_.nmu_minus_  = 0;
    evt_.nmu_plus_   = 0;

    //---------------------------------------------------
    // Add tracks to the output collections

    // Reset the collections
    tracks_  .reset();
    electron_.reset();
    muon_    .reset();
    proton_  .reset();

    // Loop through the tracks
    if(trk_) {
      for(Int_t index = 0; index < evt_.ntrks_; ++index) {
        const auto& trk = (*trk_)[index];
        // reject bad tracks
        if(trk.status < 0 || trk.goodfit == 0) continue;

        // all tracks with accepted statuses
        AddTrack(tracks_, index);
        // electron tracks
        if(std::abs(trk.pdg) == 11) AddTrack(electron_, index);
        // muon tracks
        if(std::abs(trk.pdg) == 13) AddTrack(muon_, index);
        // proton tracks
        if(std::abs(trk.pdg) == 2212) AddTrack(proton_, index);
        if(trk.pdg ==  11) ++evt_.nelectrons_;
        if(trk.pdg == -11) ++evt_.npositrons_;
        if(trk.pdg ==  13) ++evt_.nmu_minus_;
        if(trk.pdg == -13) ++evt_.nmu_plus_;
      }
    }
  }

  //------------------------------------------------------------------------------------
  // Main event-by-event processing
  bool Mu2eEvtAna::ProcessEvent() {
    if(verbose_ > 4) {
      printf("Mu2eEvtAna::%s: Printing event information:\n", __func__);
      if(trk_) {
        for(auto trk : *trk_) {
          printf("    Track: status = %i, pdg = %i\n", trk.status, trk.pdg);
        }
      } else { //no track collection is defined
        printf("  --> Not defined!\n");
      }
    }
    if(!trk_) return false; //require the collection to be defined

    FillEventHist(evt_hists_[0]); //all events with well defined inputs

    return true;
  }

  //------------------------------------------------------------------------------------
  // Main processing loop
  int Mu2eEvtAna::Process(Long64_t nentries, Long64_t first) {

    //---------------------------------------------------
    // Validate the input

    if(!ntuple_) {
      printf("Mu2eEvtAna::%s: No ntuple is initialized to process!\n", __func__);
      return -1;
    }

    const auto entries = ntuple_->GetEntries();
    if(entries == 0) {
      if(verbose_ > -1) printf("Mu2eEvtAna::%s: No entries to process\n", __func__);
      return 1;
    }
    if(first >= entries) {
      if(verbose_ > -2) printf("Mu2eEvtAna::%s: Start entry (%llu) is greater than or equal to the number of entries (%llu)\n",
                               __func__, first, entries);
      return 2;
    }

    //---------------------------------------------------
    // Initialize the input
    if(InitializeInput()) return 3;

    //---------------------------------------------------
    // Initialize output file/ntuple structures
    if(InitializeOutput()) return 4;

    //---------------------------------------------------
    // Process each requested event, storing accepted events

    Long64_t nseen     = 0; //track processing
    Long64_t naccepted = 0; //track the acceptance rate
    Long64_t nnegative = 0; //track the number of negative weight events

    const Long64_t max_entry = (nentries < 0) ? entries : std::min(entries, first+nentries);

    for(Long64_t entry = first; entry < max_entry; ++entry) {
      entry_ = entry;
      ntuple_->GetEntry(entry);
      if((verbose_ > -1 && nseen % report_rate_ == 0) || verbose_ > 5) {
        printf("Mu2eEvtAna::%s: Processing event %7lld (entry %8lld, event %6i/%6i/%6i): N(accept) = %7lld (%6.2f%%)\n", __func__, nseen, entry,
               evt_info_->run,evt_info_->subrun,evt_info_->event, naccepted, naccepted*100./((nseen <= 0) ? 1 : nseen));
      }
      ++nseen;

      // Initialize event information
      InitializeEvent();

      // Decide whether or not to accept the event
      if(ProcessEvent()) {
        FillOutput();
        ++naccepted;
      }
    }

    //---------------------------------------------------
    // Store the normalization information for the ntuple

    // FIXME: Figure out how to access gen count objects
    norm_.ngen_    = entries;
    norm_.nntuple_ = entries;
    norm_.nseen_   = nseen;
    norm_.naccept_ = naccepted;
    norm_.nneg_    = nnegative;
    tnorm_->Fill();

    // Close output structures and report results
    tout_->Write();
    tnorm_->Write();
    fout_->Write();
    fout_->Close();

    if(verbose_ > -2) printf("Mu2eEvtAna::%s: Processed %8lld events, accepted %8lld (%6.2f%%)\n", __func__, nseen, naccepted,
                             naccepted*100./((nseen <= 0) ? 1 : nseen));

    return 0;
  }

}
