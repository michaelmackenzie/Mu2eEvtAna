#include "Mu2eEvtAna/inc/Mu2eEvtAna.hh"

using namespace mu2e;
namespace Mu2eEvtAna {

  //------------------------------------------------------------------------------------
  // Constructor
  Mu2eEvtAna::Mu2eEvtAna(int verbose) : ntuple_(nullptr),
                                        fout_(nullptr), tout_(nullptr), tnorm_(nullptr), name_("test"),
                                        report_rate_(1000), verbose_(verbose) {
    //initialize the arrays to null
    for(int ihist = 0; ihist < kMaxEventHists; ++ihist) {
      evt_hists_[ihist] = nullptr;
      trk_hists_[ihist] = nullptr;
    }

    //Default histogram selections
    evt_hists_[0] = new EventHist_t;
    trk_hists_[0] = new TrackHist_t;

    trk_hists_[1] = new TrackHist_t; // good electron tracks
  }

  //------------------------------------------------------------------------------------
  // Retrieve the input ntuple from the given file/file list
  int Mu2eEvtAna::AddFile(TString file_name, Long64_t max_entries) {
    if(!ntuple_) ntuple_ = new TChain("EventNtuple/ntuple");
    // Check if the given filename contains .root at the end
    if (file_name.EndsWith(".root")) { // assume it's a single file FIXME: Allow for wildcards
      ntuple_->Add(file_name.Data());
    } else { // assume it's a file list
      std::ifstream filelist(file_name.Data());
      if (filelist.is_open()) {
        // Count the number of input files in the text file
        const int nfiles = std::count(std::istreambuf_iterator<char>(filelist),
                                      std::istreambuf_iterator<char>(), '\n');
        // Clear any error flags and reset to the beginning
        filelist.clear();
        filelist.seekg(0, std::ios::beg);

        // Add each file to the TChain
        if(verbose_ > -1) printf("%s: Loading file list %s with %i files\n", __func__, file_name.Data(), nfiles);
        std::string line;
        int ifile = 0;
        while (std::getline(filelist, line)) {
          ++ifile;
          if(verbose_ > -1 && (ifile-1) % 10 == 0) {printf("\r%s: Loading file %3i (%.1f%%)", __func__, ifile, ifile*100./nfiles); fflush(stdout);}
          ntuple_->Add(line.c_str());
          if(max_entries > 0 && ntuple_->GetEntries() > max_entries) {
            if(verbose_ > -1) printf("%s: Loaded %i files of %i with %llu entries\n", __func__, ifile, nfiles, ntuple_->GetEntries());
            break;
          }
        }
        filelist.close();
      } else {
        printf("%s: Error! Unable to read input file list %s\n", __func__, file_name.Data());
        return 1;
      }
    }
    return 0;
  }

  //------------------------------------------------------------------------------------
  // Initialize the input ntuple information
  int Mu2eEvtAna::InitializeInput() {
    if(verbose_ > 1) printf("Mu2eEvtAna::%s: Adding input information\n", __func__);
    if(!ntuple_) {
      if(verbose_ > -2) printf("Mu2eEvtAna::%s: No ntuple is defined to configure!\n", __func__);
      return -1;
    }
    event_ = new Event(ntuple_);
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
    Hist->hNGoodTrks = new TH1D("ngoodtrks", "N(good tracks)"  ,  10, 0,  10);
    Hist->hNIDTrks   = new TH1D("nidtrks", "N(accepted tracks)"  ,  10, 0,  10);
    Hist->hNDigis    = new TH1D("ndigis"   , "N(digis)"   , 200, 0, 200);
    Hist->hNClusters = new TH1D("nclusters", "N(clusters)",  10, 0,  10);
  }

  //------------------------------------------------------------------------------------
  // Initialize the histograms for a track selection
  void Mu2eEvtAna::BookTrackHist(TrackHist_t* Hist) {
    if(!Hist) {
      throw std::runtime_error("Attempting to book histograms in a null TrackHist_t\n");
    }
    Hist->fP[0]     = new TH1F("p"        , "Track momentum", 600, -150,  150);
    Hist->fP[1]     = new TH1F("p_2"      , "Track momentum", 400,   80,  120);
  }

  //------------------------------------------------------------------------------------
  // Initialize the histogram sets
  void Mu2eEvtAna::BookHistograms(TDirectory* dir) {

    for(int ihist = 0; ihist < kMaxEventHists; ++ihist) {
      if(evt_hists_[ihist]) {
        auto subdir = dir->mkdir(Form("evt_%i", ihist));
        subdir->cd();
        BookEventHist(evt_hists_[ihist]);
        dir->cd();
        evt_dirs_[ihist] = subdir;
      }
      if(trk_hists_[ihist]) {
        auto subdir = dir->mkdir(Form("trk_%i", ihist));
        subdir->cd();
        BookTrackHist(trk_hists_[ihist]);
        dir->cd();
        trk_dirs_[ihist] = subdir;
      }
    }

  }

  //------------------------------------------------------------------------------------
  // Fill the histograms for an event selection
  void Mu2eEvtAna::FillEventHist(EventHist_t* Hist) {
    if(!Hist) {
      throw std::runtime_error("Attempting to fill histograms in a null EventHist_t\n");
    }
    Hist->hNTrks     ->Fill(evt_.ntracks_  , evt_.weight_);
    Hist->hNGoodTrks ->Fill(evt_.ngoodtrks_  , evt_.weight_);
    Hist->hNIDTrks   ->Fill(evt_.ntrks_id_   , evt_.weight_);
    Hist->hNDigis    ->Fill(evt_.ndigis_   , evt_.weight_);
    Hist->hNClusters ->Fill(evt_.nclusters_, evt_.weight_);
  }

  //------------------------------------------------------------------------------------
  // Fill the histograms for a track selection
  void Mu2eEvtAna::FillTrackHist(TrackHist_t* Hist, Track_t* Track) {
    if(!Hist) {
      throw std::runtime_error(Form("Mu2eEvtAna::%s: Attempting to fill histograms in a null TrackHist_t\n", __func__));
    }
    if(!Track) {
      if(verbose_ > 0) printf("Mu2eEvtAna::%s: Filling track histogram set with null track par\n", __func__);
      return;
    }
    auto trk = Track->track_;
    if(!trk) {
      if(verbose_ > 0) printf("Mu2eEvtAna::%s: Filling track histogram set with null track\n", __func__);
      return;
    }
    Hist->fP[0] ->Fill(Track->PFront(), evt_.weight_);
    Hist->fP[1] ->Fill(Track->PFront(), evt_.weight_);
  }

  //------------------------------------------------------------------------------------
  // Setup the output ntuple structure
  void Mu2eEvtAna::AddOutputBranches(TTree* t) {

    // Event branches

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
  // Initialize event information
  void Mu2eEvtAna::InitializeEvent() {
    //---------------------------------------------------
    // Add tracks to the output collections

    // Reset the event info
    evt_.ntracks_    = 0;
    evt_.ngoodtrks_  = 0;
    evt_.ntrks_id_   = 0;
    evt_.nelectrons_ = 0;
    evt_.nmuons_     = 0;
    evt_.nprotons_   = 0;

    // Loop through the tracks
    auto tracks = event_->GetTracks();
    evt_.ntracks_ = tracks.size();
    for(int itrk = 0; itrk < evt_.ntracks_; ++itrk) {
      InitTrack(&tracks[itrk], tracks_[itrk]);
      if(tracks_[itrk].IsGood()) ++evt_.ngoodtrks_;
      if(verbose_ > 1) tracks_[itrk].Print((itrk == 0) ? "banner" : "");
    }

    if(verbose_ > 4) {
      printf("Mu2eEvtAna::%s: Printing event information:\n", __func__);
      printf(" N(tracks) = %2i\n", evt_.ntracks_);
    }
  }

  //------------------------------------------------------------------------------------
  // Initialize track information
  void Mu2eEvtAna::InitTrack(Track* track, Track_t& trk_par) {
    trk_par.Reset();
    trk_par.track_ = track;
    if(!track) return;
  }

  //------------------------------------------------------------------------------------
  // Main event-by-event processing
  bool Mu2eEvtAna::ProcessEvent() {
    if(verbose_ > 4) {
      printf("Mu2eEvtAna::%s: Printing event information:\n", __func__);
    }
    FillEventHist(evt_hists_[0]); //all events with well defined inputs
    for(int itrk = 0; itrk < evt_.ntracks_; ++itrk) { // all tracks
      FillTrackHist(trk_hists_[0], &tracks_[itrk]);
      if(tracks_[itrk].IsGood() && tracks_[itrk].FitPDG() == -11) FillTrackHist(trk_hists_[1], &tracks_[itrk]);
    }

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
      event_->Update(verbose_ > 3);
      if((verbose_ > -1 && nseen % report_rate_ == 0) || verbose_ > 1) {
        printf("Mu2eEvtAna::%s: Processing event %7lld (entry %8lld, event %6i/%7i/%7i): N(accept) = %7lld (%6.2f%%)\n", __func__, nseen, entry,
               event_->evtinfo->run,event_->evtinfo->subrun,event_->evtinfo->event, naccepted, naccepted*100./((nseen <= 0) ? 1 : nseen));
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
