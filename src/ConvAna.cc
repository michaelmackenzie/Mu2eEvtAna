#include "Mu2eEvtAna/inc/ConvAna.hh"

using namespace mu2e;
namespace Mu2eEvtAna {

  //------------------------------------------------------------------------------------
  // Constructor
  ConvAna::ConvAna(int verbose) : Mu2eEvtAna(verbose) {
    //initialize the arrays to null
    for(int ihist = 0; ihist < kMaxHists; ++ihist) {
      sys_hists_[ihist] = nullptr;
    }

    // Initialize the MVA models

    trkqual_ = new TMVA::Reader("!Color:!Silent");
    MVATools::InitializeVariables(*trkqual_, "TrkQual", tree_, trkqual_version_);
    trkqual_->BookMVA("TrkQual", "Mu2eEvtAna/data/trkqual_MLP.weights.xml");

    pid_ = new TMVA::Reader("!Color:!Silent");
    MVATools::InitializeVariables(*pid_, "PID", tree_, pid_version_);
    pid_->BookMVA("PID", "Mu2eEvtAna/data/pid_MLP.weights.xml");

    trkpid_ = new TMVA::Reader("!Color:!Silent");
    MVATools::InitializeVariables(*trkpid_, "TrkPID", tree_, trkpid_version_);
    trkpid_->BookMVA("TrkPID", "Mu2eEvtAna/data/trkpid_MLP.weights.xml");

    cosmic_id_ = new TMVA::Reader("!Color:!Silent");
    MVATools::InitializeVariables(*cosmic_id_, "CosmicID", tree_, cosmic_id_version_);
    cosmic_id_->BookMVA("CosmicID", "Mu2eEvtAna/data/cosmicid_MLP.weights.xml");
  }


  //------------------------------------------------------------------------------------
  // Define the histogram selections
  void ConvAna::InitHistSelections() {

    //-----------------------------------------------------------------------------
    // book histogram selections
    //-----------------------------------------------------------------------------
    struct hist_info_t {
      TString _dsc; // description of the selection
      bool    _trk; // track histograms
      bool    _hlx; // helix histograms
      bool    _smp; // sim particle histograms
      bool    _gnp; // gen particle histograms
      bool    _crv; // CRV histograms
      bool    _sys; // systematic histograms
      bool    _crs; // control regions included
      bool    _trs; // output trees
      hist_info_t(TString dsc = "", bool trk = false, bool hlx = false, bool smp = false, bool gnp = false,
                  bool crv = false, bool sys = false, bool crs = false, bool trs = false)
        : _dsc(dsc), _trk(trk), _hlx(hlx), _smp(smp), _gnp(gnp), _crv(crv), _sys(sys),
          _crs(crs), _trs(trs) {}
    };

    hist_info_t* hist_sets[kMaxHists];
    for (int i=0; i<kMaxHists; i++) {
      hist_sets[i] = nullptr;
    }

    //                                 description                         trk    hlx    simp   genp   crv    sys    crs    trs
    hist_sets[  0] = new hist_info_t("All events, Offline track"        ,  true,  true,  true,  true,  true, false, false, false);
    hist_sets[  1] = new hist_info_t("All events, APR track"            ,  true,  true,  true,  true, false, false, false, false);
    hist_sets[  2] = new hist_info_t("All events, CPR track"            ,  true,  true,  true,  true, false, false, false, false);
    hist_sets[  3] = new hist_info_t("No weights, Offline track"        ,  true,  true,  true,  true, false, false, false, false);
    hist_sets[  4] = new hist_info_t("e- and p > 95"                    ,  true, false,  true, false, false, false, false, false);
    hist_sets[  5] = new hist_info_t("e+ and p > 75"                    ,  true, false,  true, false, false, false, false, false);
    hist_sets[  6] = new hist_info_t("Gen(E) > 95"                      ,  true, false,  true, false, false, false, false, false);
    hist_sets[  7] = new hist_info_t("e+-, TrkID, no MC cut"            ,  true, false,  true, false, false, false, false, false);
    hist_sets[  8] = new hist_info_t("Events with a track"              ,  true, false,  true, false, false, false, false, false);
    hist_sets[ 10] = new hist_info_t("Offline, event selection"         ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 11] = new hist_info_t("e+/-: trigger"                    ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 12] = new hist_info_t("e+/-: failed trigger"             ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 13] = new hist_info_t("e+/-: negative"                   ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 14] = new hist_info_t("e+/-: positive"                   ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 15] = new hist_info_t("e+/-: no weights"                 ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 16] = new hist_info_t("e-: p > 95"                       ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 20] = new hist_info_t("e-: full window"                  ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 21] = new hist_info_t("e-: narrow window"                ,  true,  true,  true,  true, false, false,  true, false);
    hist_sets[ 22] = new hist_info_t("e-: full window, no weights"      ,  true,  true,  true,  true, false, false,  true, false);
    hist_sets[ 23] = new hist_info_t("e-: high error"                   ,  true,  true,  true,  true, false, false,  true, false);
    hist_sets[ 25] = new hist_info_t("e-: mom up"                       ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 26] = new hist_info_t("e-: mom down"                     ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 30] = new hist_info_t("e-: no CRV veto"                  ,  true, false, false, false,  true, false,  true, false);
    hist_sets[ 31] = new hist_info_t("e+-: no CRV veto or p cut"        ,  true, false, false, false,  true, false,  true, false);
    hist_sets[ 34] = new hist_info_t("e-: low dP(ST)"                   ,  true, false, false, false, false,  true,  true, false);
    hist_sets[ 35] = new hist_info_t("e-: high dP(ST)"                  ,  true, false, false, false, false,  true,  true, false);
    hist_sets[ 40] = new hist_info_t("e+: full window"                  ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 41] = new hist_info_t("e+: narrow window"                ,  true,  true,  true,  true, false, false,  true, false);
    hist_sets[ 42] = new hist_info_t("e+: full window, no weights"      ,  true,  true,  true,  true, false, false,  true, false);
    hist_sets[ 50] = new hist_info_t("e+: no CRV veto"                  ,  true, false, false, false,  true, false,  true, false);

    // CRV studies histograms
    hist_sets[ 80] = new hist_info_t("CRV: 1"                           ,  true, false, false, false,  true, false, false, false);
    hist_sets[ 81] = new hist_info_t("CRV: Ue/Umu tracks"               ,  true, false, false, false,  true, false, false, false);
    hist_sets[ 82] = new hist_info_t("CRV: calo cluster with track"     ,  true, false, false, false,  true, false, false, false);
    hist_sets[ 83] = new hist_info_t("CRV: dt < -60 ns"                 ,  true, false, false, false,  true, false, false, false);
    hist_sets[ 84] = new hist_info_t("CRV: MC electrons"                ,  true, false, false, false,  true, false, false,  true);
    hist_sets[ 85] = new hist_info_t("CRV: MC muons"                    ,  true, false, false, false,  true, false, false,  true);
    hist_sets[ 86] = new hist_info_t("CRV: Correct cluster"             ,  true, false, false, false,  true, false, false, false);
    hist_sets[ 87] = new hist_info_t("CRV: Upstream veto"               ,  true, false, false, false,  true, false, false, false);

    for (int i=0; i<kMaxHists; i++) {
      const int index = i % 1000;
      if(!hist_sets[index]) continue;
      const bool is_cr = i >= 1000 && !hist_sets[index]->_crs;
      if(is_cr) continue; // control region
      if(i >= 4000) break; //Control regions above 4000 not yet implemented
      evt_hists_[index] = new EventHist_t;
      if(hist_sets[index]->_trk) trk_hists_[index] = new TrackHist_t;
      if(hist_sets[index]->_crv && ! is_cr) crv_hists_[index] = new CRVHist_t;
      if(hist_sets[index]->_sys) sys_hists_[index] = new SysHist_t;
      // FIXME: Add missing histogram types
    }
  }

  //------------------------------------------------------------------------------------
  // Initialize systematic histograms
  void ConvAna::BookSystematicHist(SysHist_t* Hist, const char* Folder) {
    if(!Hist) {
      throw std::runtime_error("Attempting to book histograms in a null SysHist_t\n");
    }
    for(int isys = 0; isys < kMaxSystematics; ++isys) {
      // check if the systematic is defined
      TString name = systematics_.GetName(isys);
      if(name == "") continue;
      Hist->fObs[isys] = new TH1F(Form("obs_%i", isys),Form("%s: Systematic %s",Folder, name.Data()), 300, 80., 110.); //FIXME: This should inherit from the nominal observable binning
      // For debug investigations
      if(fill_verbose_sys_) {
        Hist->fDeltaObs   [isys] = new TH1F(Form("delta_obs_%i"   , isys),Form("%s: Systematic %s: #DeltaObs"   ,Folder, name.Data()), 100, -2., 2.);
        Hist->fWeight     [isys] = new TH1F(Form("weight_%i"      , isys),Form("%s: Systematic %s: Weight"      ,Folder, name.Data()), 100,  0., 2.);
        Hist->fDeltaWeight[isys] = new TH1F(Form("delta_weight_%i", isys),Form("%s: Systematic %s: #DeltaWt/Wt" ,Folder, name.Data()), 100, -1., 1.);
      }
    }
  }

  //------------------------------------------------------------------------------------
  // Initialize the histogram sets
  void ConvAna::BookHistograms(TDirectory* dir) {
    Mu2eEvtAna::BookHistograms(dir);
    for(int ihist = 0; ihist < kMaxHists; ++ihist) {
      if(sys_hists_[ihist]) {
        const char* folder = Form("sys_%i", ihist);
        auto subdir = dir->mkdir(folder);
        subdir->cd();
        BookSystematicHist(sys_hists_[ihist], folder);
        dir->cd();
        sys_dirs_[ihist] = subdir;
      }
    }
  }

  //------------------------------------------------------------------------------------
  // Initialize the input ntuple information
  int ConvAna::InitializeInput() {
    Mu2eEvtAna::InitializeInput();
    return 0;
  }

  //------------------------------------------------------------------------------------
  // Initialize the output ntuple information
  int ConvAna::InitializeOutput() {
    Mu2eEvtAna::InitializeOutput();
    return 0;
  }

  //------------------------------------------------------------------------------------
  // Initialize event information
  void ConvAna::InitializeEvent() {
    Mu2eEvtAna::InitializeEvent();
    track_ = nullptr;
    calo_cluster_ = nullptr;
    crv_cluster_ = nullptr;
  }

  //------------------------------------------------------------------------------------
  // Initialize track information
  void ConvAna::InitTrack(rooutil::Track* track, Track_t& trk_par) {
    Mu2eEvtAna::InitTrack(track, trk_par);

    // Initialize MVA scores
    if(evaluate_mvas_) {
      track_ = &trk_par;
      InitTreeData();
      watch_->SetTime("MVAs");
      trk_par.trkqual_ = (trkqual_) ? trkqual_->EvaluateMVA("TrkQual") : -999.f;
      trk_par.pid_ = (pid_) ? pid_->EvaluateMVA("PID") : -999.f;
      trk_par.trkpid_ = (trkpid_) ? trkpid_->EvaluateMVA("TrkPID") : -999.f;
      trk_par.cosmic_id_ = (cosmic_id_) ? cosmic_id_->EvaluateMVA("CosmicID") : -999.f;
      watch_->StopTime("MVAs");
      // Reset the main ID with the new MVA scores
      trk_par.SetID(TrackID(&trk_par), 0);
    }
  }

  //------------------------------------------------------------------------------------
  // Initialize the tree structure with selected object info
  void ConvAna::InitTreeData() {
    tree_.Reset();

    if(track_) {
      tree_.fTrkP = track_->PFront();
      tree_.fTrkT0 = track_->TFront();
      tree_.fTrkD0 = track_->D0Front();
      tree_.fTrkTanDip = track_->TanDipFront();
      tree_.fTrkCosTheta = track_->CosThetaFront();
      tree_.fTrkFitCon = track_->FitCon();
      tree_.fTrkLogFitCon = (track_->FitCon() > 0.) ? log10(track_->FitCon()) : -100.f;
      tree_.fTrkRMax = track_->RMaxFront();
      tree_.fTrkCluster = track_->ECluster();
      tree_.fTrkEP = track_->EPFront();
      tree_.fTrkDt = track_->Dt();
      tree_.fTrkActiveRatio = track_->NActive() * 1.f / track_->NHits();
      tree_.fTrkNullRatio = track_->NNull() * 1.f / track_->NHits();
      tree_.fTrkTZSlope = 0.f; //fTrkPar.fTZSlope;
      tree_.fTrkTZSlopeSig = 0.f; //fTrkPar.TZSlopeSig();
      tree_.fTrkTZSlopeRatio = 1.f; //fTrkPar.TZSlopeRatio();
      tree_.fTrkPExitDiff = track_->PFront() - track_->PBack();
      tree_.fTrkQual = track_->TrkQual();
      tree_.fTrkPID = track_->PID();
      tree_.fTrkOnlyPID = 1.f; //fTrkPar.fTrkPID;
      tree_.fTrkCosmicID = 1.f; //fTrkPar.fCosmicID;
      tree_.fTrkCharge = track_->Charge();
      tree_.fTrkMCDp = track_->MCDeltaPFront();
      tree_.fTrkMCPDG = track_->MCPDG();

      // For TrkQual
      tree_.fTrkQual_nactive = track_->NActive();
      tree_.fTrkQual_activehitsfraction = tree_.fTrkActiveRatio;
      tree_.fTrkQual_nullhitsfraction = tree_.fTrkNullRatio;
      tree_.fTrkQual_activematsitesfraction = track_->NMatActive() * 1.f / track_->NActive();
      tree_.fTrkQual_fitcons = tree_.fTrkFitCon;
      tree_.fTrkQual_momerr = track_->MomErrFront();
      tree_.fTrkQual_t0err = track_->TErrFront();
    }

    // if(crv_cluster_) {
    //   tree_.fCRVZ = fTrkPar.fCRVStubPar->fZ;
    //   tree_.fCRVDeltaT = fTrkPar.fCRVStubPar->fCorrTime - fTrack->fT0;
    //   tree_.fCRVNPulses = fTrkPar.fCRVStubPar->fCluster->NPulses();
    //   tree_.fCRVNPe = fTrkPar.fCRVStubPar->fCluster->NPe();
    //   tree_.fCRVNPePP = fTrkPar.fCRVStubPar->fNPePP;
    // } else {
    //   tree_.fCRVZ = 0.f;
    //   tree_.fCRVDeltaT = 0.f;
    //   tree_.fCRVNPulses = 0.f;
    //   tree_.fCRVNPe = 0.f;
    //   tree_.fCRVNPePP = 0.f;
    // }

    tree_.fWeight = evt_.weight_;

    // // FIXME: For now just doing 50% splitting by ID
    // const int event = GetEvent()->fEventNumber;
    // tree_.fTrain = (event % 2 == 0) ? 1.f : -1.f;

  }

  //------------------------------------------------------------------------------------
  // Main event-by-event processing
  bool ConvAna::ProcessEvent() {
    FillEventHist(evt_hists_[0]); //all events with well defined inputs

    // Loop through the track collection
    for(int itrk = 0; itrk < evt_.ntracks_; ++itrk) {
      track_ = &tracks_[itrk];
      FillTrackHist(trk_hists_[0], track_); // all tracks

      // Downstream electron sets
      if(track_->IsGood() && track_->FitPDG() == 11 && track_->PZFront() > 0.f) {
        if(track_->PFront() > 95.f && track_->Charge() < 0) FillTrackHist(trk_hists_[4], track_);
        const int ID = track_->ID();
        const int event_id = evt_.nde_tracks_ != 1; // FIXME: Implement

        // Offset to control regions
        int set_offset(0);
        if(ID & (1 << kCRV))
          set_offset += kCRVVetoOffset;
        if(ID & (1 << kT0))
          set_offset += kTimeCutOffset;
        // const int id_no_crv = ID & (~(1 << kCRV)); // ID without the CRV coincidence cluster considered
        // const int id_no_time = ID & (~(1 << kT0)); // ID with a looser timing cut
        // const int id_no_crv_time = id_no_crv & id_no_time;

        if(event_id == 0 && ID == 0) { //id_no_crv_time == 0) {
          FillTrackHist(trk_hists_[10 + set_offset], track_);
          if(track_->Charge() < 0) { // electron selection
            if(track_->PFront() > 100.f && track_->PFront() < 110.f) {
              FillTrackHist(trk_hists_[20 + set_offset], track_);
              if(track_->PFront() > 103.6f && track_->PFront() < 104.9f) {
                FillTrackHist(trk_hists_[21 + set_offset], track_);
              }
            }
          }
        }
      }
    }

    if(evt_.nde_tracks_ != 1) return false; //exactly one positron or electron
    FillEventHist(evt_hists_[1]);
    return true;
  }
}
