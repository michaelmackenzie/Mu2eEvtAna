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
    hist_sets[ 24] = new hist_info_t("e-:  alt ID"                      ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 25] = new hist_info_t("e-: !alt ID"                      ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 30] = new hist_info_t("e-: no CRV veto"                  ,  true, false, false, false,  true, false,  true, false);
    hist_sets[ 31] = new hist_info_t("e+-: no CRV veto or p cut"        ,  true, false, false, false,  true, false,  true, false);
    hist_sets[ 34] = new hist_info_t("e-: low dP(ST)"                   ,  true, false, false, false, false,  true,  true, false);
    hist_sets[ 35] = new hist_info_t("e-: high dP(ST)"                  ,  true, false, false, false, false,  true,  true, false);
    hist_sets[ 40] = new hist_info_t("e+: full window"                  ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 41] = new hist_info_t("e+: narrow window"                ,  true,  true,  true,  true, false, false,  true, false);
    hist_sets[ 42] = new hist_info_t("e+: full window, no weights"      ,  true,  true,  true,  true, false, false,  true, false);
    hist_sets[ 50] = new hist_info_t("e+: no CRV veto"                  ,  true, false, false, false,  true, false,  true, false);
    hist_sets[ 60] = new hist_info_t("e-: Run 1A ID"                    ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 61] = new hist_info_t("e-: Run 1A ID, no weights"        ,  true,  true,  true,  true,  true, false,  true, false);
    hist_sets[ 62] = new hist_info_t("e-: Run 1A ID + upstream veto"    ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 65] = new hist_info_t("e-: Run 1A ID, loose time"        ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 66] = new hist_info_t("e-: Run 1A ID, cut-and-count"     ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 67] = new hist_info_t("e-: Run 1A ID, no tan dip"        ,  true,  true,  true,  true,  true,  true,  true, false);
    hist_sets[ 70] = new hist_info_t("e-: cut-flow ID"                  ,  true, false, false, false, false, false, false, false);
    hist_sets[ 73] = new hist_info_t("e-: Optimized cut-set"            ,  true,  true,  true,  true,  true,  true,  true,  true);

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
      const int index = i % 1000; // base index, ignoring control region offset
      if(!hist_sets[index]) continue;
      const bool is_cr = i >= 1000 && !hist_sets[index]->_crs;
      if(is_cr && !hist_sets[index]->_crs) continue; // control region
      if(i >= 4000) break; //Control regions above 4000 not yet implemented
      evt_hists_[i] = new EventHist_t;
      if(hist_sets[index]->_trk) trk_hists_[i] = new TrackHist_t;
      if(hist_sets[index]->_crv && ! is_cr) crv_hists_[i] = new CRVHist_t;
      if(hist_sets[index]->_sys) sys_hists_[i] = new SysHist_t;
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
  void ConvAna::InitTrack(const rooutil::Track* track, Track_t& trk_par) {
    Mu2eEvtAna::InitTrack(track, trk_par);
    if(!trk_par.IsGood()) return;

    // Initialize MVA scores
    if(evaluate_mvas_) {
      track_ = &trk_par;
      InitTreeData();
      ValidateVariable(tree_.fTrkQual_nactive, "TrkQual_nactive");
      ValidateVariable(tree_.fTrkQual_activehitsfraction, "TrkQual_activehitsfraction");
      ValidateVariable(tree_.fTrkQual_nullhitsfraction, "TrkQual_nullhitsfraction");
      ValidateVariable(tree_.fTrkQual_activematsitesfraction, "TrkQual_activematsitesfraction");
      ValidateVariable(tree_.fTrkQual_fitcons, "TrkQual_fitcons");
      ValidateVariable(tree_.fTrkQual_momerr, "TrkQual_momerr");
      ValidateVariable(tree_.fTrkQual_t0err, "TrkQual_t0err");
      ValidateVariable(tree_.fTrkEP, "TrkEP");
      ValidateVariable(tree_.fTrkDt, "TrkDt");
      ValidateVariable(tree_.fTrkFitCon, "TrkFitCon");
      ValidateVariable(tree_.fTrkLogFitCon, "TrkLogFitCon");
      ValidateVariable(tree_.fTrkActiveRatio, "TrkActiveRatio");
      ValidateVariable(tree_.fTrkNullRatio, "TrkNullRatio");
      ValidateVariable(tree_.fTrkTZSlope, "TrkTZSlope");
      ValidateVariable(tree_.fTrkTZSlopeSig, "TrkTZSlopeSig");
      ValidateVariable(tree_.fTrkTZSlopeRatio, "TrkTZSlopeRatio");
      ValidateVariable(tree_.fTrkD0, "TrkD0");
      ValidateVariable(tree_.fTrkTanDip, "TrkTanDip");
      ValidateVariable(tree_.fTrkCosTheta, "TrkCosTheta");
      ValidateVariable(tree_.fTrkRMax, "TrkRMax");
      watch_->SetTime("MVAs");
      trk_par.trkqual_ = (trkqual_) ? trkqual_->EvaluateMVA("TrkQual") : -999.f;
      trk_par.pid_ = (pid_) ? pid_->EvaluateMVA("PID") : -999.f;
      trk_par.trkpid_ = (trkpid_) ? trkpid_->EvaluateMVA("TrkPID") : -999.f;
      trk_par.cosmic_id_ = (cosmic_id_) ? cosmic_id_->EvaluateMVA("CosmicID") : -999.f;
      watch_->StopTime("MVAs");
      // Reset the main ID with the new MVA scores
      trk_par.SetID(TrackID(&trk_par), 0);
    }

    // Standard Run 1A selection set
    trk_par.SetID(Run1ATrackID(&trk_par), 1);
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
      tree_.fTrkTZSlope = track_->TZSlope();
      tree_.fTrkTZSlopeSig = track_->TZSlopeSig();
      tree_.fTrkTZSlopeRatio = track_->TZSlopeRatio();
      tree_.fTrkPExitDiff = track_->PFront() - track_->PBack();
      tree_.fTrkQual = track_->TrkQual();
      tree_.fTrkPID = track_->PID();
      tree_.fTrkOnlyPID = track_->TrkPID();
      tree_.fTrkCosmicID = track_->CosmicID();
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
  // Track selection
  int ConvAna::Run1ATrackID(Track_t* track) {
    if(!track || !track->track_) return 0;
    int ID(0);
    const float pmin = (track->Charge() < 0) ?  97.f :  85.f;
    const float pmax = (track->Charge() < 0) ? 110.f : 110.f;
    if(track->PFront() < pmin || track->PFront() > pmax)          ID += 1 << kP;        // Loose momentum window
    if(track->OPAInter())                                         ID += 1 << kRMax;     // Cosmic rejection
    if(std::abs(track->FitPDG()) != 11 || track->PZFront() < 0.f) ID += 1 << kFitHyp;   // Downstream electron fit
    if(track->TrkQual() > -10. && track->TrkQual() < 0.2)         ID += 1 << kTrkQual;  // Track quality
    else if(track->TErrFront() > 0.9)                             ID += 1 << kTrkQual;  // FIXME: Give a label
    else if(track->NActive() < 20)                                ID += 1 << kTrkQual;  // FIXME: Give a label
    if(track->TFront() < 475. || track->TFront() > 1650.)         ID += 1 << kT0;       // FIXME: Check selection
    if(track->ECluster() <= 0.)                                   ID += 1 << kClusterE; // Require a cluster
    if(track->NSTInter() == 0)                                    ID += 1 << kD0;       // Consistent with stopping target
    if(track->TanDipFront() < 0.5 || track->TanDipFront() > 0.88) ID += 1 << kTDip;     // Tan(dip) = pz / pt here
    // if(track->TanDipFront() < 0.5 || track->TanDipFront() > 0.95) ID += 1 << kTDip;     // Tan(dip) = pz / pt here
    if(track->TFront() < 475. || track->TFront() > 1650.)         ID += 1 << kT0Loose;  // For RPC control regions

    // PID rejection
    if(track->ECluster() > 0.f) {
      if(track->PID() < -100.f) { // no score
        if(track->EPFront() < 0.65f)                              ID += 1 << kPID;
      } else if(track->PID() < 0.67f)                             ID += 1 << kPID;
    }

    // CRV rejection
    bool fail_crv = false;
    // Look for associated stub
    if(track->stub_) {
      auto stub = track->stub_;
      const float deltat_crv    = track->TFront() - stub->Time();
      if(deltat_crv > -150.f && deltat_crv < 150.f) fail_crv = true;
    }
    // Look for good CRV clusters anywhere in the event
    if(!fail_crv) {
      for(int icrv = 0; icrv < evt_.ncrv_clusters_; ++icrv) {
        CRVCluster_t* stub = &crv_clusters_[icrv];
        if(stub->PEs() >= 25. && stub->NHits() >= 15 && stub->TimeSpan() < 175.) {
          fail_crv = true;
          break;
        }
      }
    }
    if(fail_crv)                                                  ID += 1 << kCRV;

    // Upstream track rejection
    // FIXME: Add the standard matching logic

    // if(track->PFront() > 95. &&
    //    std::abs(track->FitPDG()) == 11 &&
    //    track->PZFront() > 0. &&
    //    ((ID & ~kT0) & ~kClusterE) != 0) {
    //   std::cout << __func__ << ": ID = " << std::bitset<32>(ID) << std::endl;
    // }
    return ID;
  }

  //------------------------------------------------------------------------------------
  // Evaluate the Run 1A selection cut-flow
  bool ConvAna::Run1ACutFlow() {
    if(!track_) return false;

    if(track_->FitPDG() != 11) return false;
    if(track_->Charge() > 0)   return false;
    run1a_cut_flow_.Increment("is_reco_electron");

    // Downstream electron sets
    if(track_->PZFront() <= 0.f) return false;
    run1a_cut_flow_.Increment("has_downstream");
    run1a_cut_flow_.Increment("upstream_veto");
    run1a_cut_flow_.Increment("trk_front_seg");

    const int Run1AID = track_->ID(1);
    const int event_id = (
                          1*(evt_.nde_tracks_ != 1) +
                          2*(!trigger_.FiredAPR() && !trigger_.FiredCPR())
                          );

    if(track_->PID() > 0.67) {
      run1a_cut_flow_.Increment("good_trkpid");
      if(track_->TrkQual() > 0.2) {
        run1a_cut_flow_.Increment("good_trkqual");
        if(track_->TErrFront() < 0.9) {
          run1a_cut_flow_.Increment("within_t0err");
          if(track_->NActive() >= 20) {
            run1a_cut_flow_.Increment("has_hits");
            if(track_->NSTInter() > 0) {
              run1a_cut_flow_.Increment("has_st");
              if(track_->OPAInter() == 0) {
                run1a_cut_flow_.Increment("no_opa");
                bool fail_crv = false;
                for(int icrv = 0; icrv < evt_.ncrv_clusters_; ++icrv) {
                  CRVCluster_t* stub = &crv_clusters_[icrv];
                  if(stub->PEs() >= 25. && stub->NHits() >= 15 && stub->TimeSpan() < 175.) {
                    fail_crv = true;
                    break;
                  }
                }
                if(!fail_crv) {
                  run1a_cut_flow_.Increment("no_crv_quality");
                  run1a_cut_flow_.Increment("no_crv_timewindow");
                  if(track_->stub_) {
                    auto stub = track_->stub_;
                    const float deltat_crv    = track_->TFront() - stub->Time();
                    if(deltat_crv > -150.f && deltat_crv < 150.f) fail_crv = true;
                  }
                  if(!fail_crv) {
                    run1a_cut_flow_.Increment("no_crv_veto");
                    if(track_->TanDipFront() > 0.5 && track_->TanDipFront() < 0.95) {
                      run1a_cut_flow_.Increment("pz_over_pt");
                      if(trigger_.FiredAPR() || trigger_.FiredCPR()) {
                        run1a_cut_flow_.Increment("good_trigger");
                        if(track_->PFront() > 103.34 && track_->PFront() < 104.74) {
                          run1a_cut_flow_.Increment("final_mom_region");
                          if(track_->TFront() > 640. && track_->TFront() < 1650.) {
                            run1a_cut_flow_.Increment("final_time_region");
                            if(Run1AID != 0 || event_id != 0) {
                              std::cout << "[ConvAna::" << __func__ << "] "
                                        << evt_.run_ << ":" << evt_.subrun_ << ":" << evt_.event_
                                        << " Event passes paper cuts but fails Run1AID = "
                                        << std::hex << Run1AID << std::dec
                                        << " or event_id = " << event_id
                                        << std::endl;
                            }
                            return true;
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
    return false;
  }

  //------------------------------------------------------------------------------------
  // Evaluate the standard selection cut-flow
  bool ConvAna::StandardCutFlow() {
    if(!track_) return false;

    const int ID = track_->ID(0);
    const int event_id = (
                          1*(evt_.nde_tracks_ != 1) +
                          2*(!trigger_.FiredAPR() && !trigger_.FiredCPR())
                          );

    if(track_->FitPDG() != 11) return false;
    if(track_->Charge() > 0)   return false;
    cut_flow_.Increment("is_reco_electron");

    // Downstream electron sets
    if(track_->PZFront() <= 0.f) return false;
    cut_flow_.Increment("has_downstream");

    if((ID & (1 << kUpstream)) != 0) return false;
    cut_flow_.Increment("upstream_veto");
    cut_flow_.Increment("trk_front_seg");
    if((ID & (1 << kPID)) != 0) return false;
    cut_flow_.Increment("good_trkpid");
    if((ID & (1 << kClusterE)) != 0) return false;
    cut_flow_.Increment("good_cluster");
    if((ID & (1 << kTrkQual)) != 0) return false;
    cut_flow_.Increment("good_trkqual");
    if((ID & (1 << kFitCon)) != 0) return false;
    cut_flow_.Increment("fitcon");
    cut_flow_.Increment("has_hits");
    if((ID & (1 << kD0)) != 0) return false;
    cut_flow_.Increment("has_st");
    if((ID & (1 << kRMax)) != 0) return false;
    cut_flow_.Increment("no_opa");
    if((ID & (1 << kCRV)) != 0) return false;
    cut_flow_.Increment("no_crv_quality");
    cut_flow_.Increment("no_crv_timewindow");
    cut_flow_.Increment("no_crv_veto");
    if((ID & (1 << kTDip)) != 0) return false;
    cut_flow_.Increment("pz_over_pt");
    if((ID & (1 << kCosmicID)) != 0) return false;
    cut_flow_.Increment("cosmic_id");
    if((event_id & (2)) != 0) return false;
    cut_flow_.Increment("trigger");
    if((ID & (1 << kP)) != 0) return false;
    cut_flow_.Increment("loose_mom_region");
    if(track_->PFront() < 103.5 || track_->PFront() > 105.) return false;
    cut_flow_.Increment("final_mom_region");
    if((ID & (1 << kT0Loose)) != 0) return false;
    cut_flow_.Increment("loose_time_region");
    if((ID & (1 << kT0)) != 0) return false;
    cut_flow_.Increment("final_time_region");
    if((event_id & 1) != 0) return false;
    cut_flow_.Increment("de_track_count");

    if(ID != 0 || event_id != 0) {
      std::cout << "[ConvAna::" << __func__ << "] "
                << evt_.run_ << ":" << evt_.subrun_ << ":" << evt_.event_
                << " Event passes cuts but fails ID = "
                << std::hex << ID << std::dec
                << " or event_id = " << event_id
                << std::endl;
    }
    return true;
  }

  //------------------------------------------------------------------------------------
  // Main event-by-event processing
  bool ConvAna::ProcessEvent() {
    ValidateTracks();

    // Cut-flow for Run 1A paper:
    // No cuts
    // is_reco_electron
    // has_downstream
    // upstream veto
    // has_trk_front_seg
    // good_trkpid (0.67
    // good_trkqual
    // within_t0err
    // has_hits
    // has_st
    // no_opa
    // no_crv_quality
    // no_crv_timewindow
    // no_crv_veto
    // 0.5 < pz_over_pt < 0.95
    // good_trigger
    // final mom region
    // final time
    // SU2020 SR

    const float nominal_weight = evt_.weight_;

    FillEventHist(evt_hists_[0]); //all events with well defined inputs

    // printf("[ConvAna::%s] Event %5i:%6i:%8i\n",
    //        __func__, evt_.run_, evt_.subrun_, evt_.event_);

    // all events
    cut_flow_.ResetEvent();
    cut_flow_.Increment("All");
    run1a_cut_flow_.ResetEvent();
    run1a_cut_flow_.Increment("All");

    // Loop through the track collection
    for(int itrk = 0; itrk < evt_.ntracks_; ++itrk) {
      track_ = &tracks_[itrk];
      // std::cout << "Trk " << itrk << ": "
      //           << tracks_[itrk].track_ << " "
      //           << tracks_[itrk].track_->trk << std::endl;
      if(verbose_ > 2) {
        std::cout << "Address (begin loop): " << track_
                  << " track = " << track_->track_
                  << " qual = " << track_->track_->trkqual << " and " << track_->TrkQual() << std::endl;
        track_->Print((itrk == 0) ? "banner" : ""); // print a second time
      }
      if(!track_->IsGood()) continue; // if not a properly fit track, skip it
      StandardCutFlow();
      if(Run1ACutFlow()) FillTrackHist(trk_hists_[70], track_);
      if(track_->FitPDG() != 11) continue; // skip muon fits for now due to rooutil bug

      // Downstream electron sets
      if(track_->FitPDG() == 11 && track_->PZFront() > 0.f) {
        if(track_->PFront() > 95.f && track_->Charge() < 0) FillTrackHist(trk_hists_[4], track_);
        const int ID = track_->ID(0);
        const int Run1AID = track_->ID(1);
        const int event_id = (
                              1*(evt_.nde_tracks_ != 1) +
                              2*(!trigger_.FiredAPR() && !trigger_.FiredCPR())
                              );


        // Offset to control regions
        int set_offset(0);
        if(ID & (1 << kCRV))
          set_offset += kCRVVetoOffset;
        if(ID & (1 << kT0))
          set_offset += kTimeCutOffset;
        const int id_no_crv = ID & (~(1 << kCRV)); // ID without the CRV coincidence cluster considered
        const int id_no_time = ID & (~(1 << kT0)); // ID with a looser timing cut
        const int id_no_crv_time = id_no_crv & id_no_time;

        if(event_id == 0 && id_no_crv_time == 0) {
          FillTrackHist(trk_hists_[10 + set_offset], track_);
          if(track_->Charge() < 0) { // electron selection
            if(track_->PFront() > 100.f && track_->PFront() < 110.f) {
              FillTrackHist(trk_hists_[20 + set_offset], track_);
              if(track_->PFront() > 103.6f && track_->PFront() < 104.9f) {
                FillTrackHist(trk_hists_[21 + set_offset], track_);
              }
              // alternate ID set on top of standard ID (~75% signal eff)
              const bool alt_id = (track_->CosThetaFront() > 0.5 && track_->CosThetaFront() < 0.646805 &&
                                   track_->TFront() > 594.049 &&
                                   track_->RMaxFront() > 482.031 &&
                                   track_->EPFront() < 0.960571 &&
                                   track_->AltPID() > 0.140625);
              if(alt_id) FillTrackHist(trk_hists_[24 + set_offset], track_);
              else       FillTrackHist(trk_hists_[25 + set_offset], track_);
            }
          }
        }

        // Standard selection set
        // Offset to control regions
        int run1a_set_offset(0);
        if(Run1AID & (1 << kCRV))
          run1a_set_offset += kCRVVetoOffset;
        if(Run1AID & (1 << kT0))
          run1a_set_offset += kTimeCutOffset;
        const int run1a_id_no_crv = Run1AID & (~(1 << kCRV)); // ID without the CRV coincidence cluster considered
        const int run1a_id_no_time = Run1AID & (~(1 << kT0)); // ID with a looser timing cut
        const int run1a_id_no_crv_time = run1a_id_no_crv & run1a_id_no_time;
        const int run1a_id_no_tdip = Run1AID & (~(1 << kTDip)); // ID without the tan dip cut

        bool upstream_veto = false;
        if(track_->upstream_) {
          auto us_trk = track_->upstream_;
          if(us_trk->FitCon() > 1.e-5 &&
             (us_trk->TrkQual() < -10. || us_trk->TrkQual() > 0.01)) upstream_veto = true;
        }

        if(event_id == 0 && run1a_id_no_crv_time == 0) {
          if(track_->TFront() > 640.) { // nominal timing window for 1D fit
            FillTrackHist(trk_hists_[60 + run1a_set_offset], track_);
            evt_.weight_ = 1.f;
            FillTrackHist(trk_hists_[61 + run1a_set_offset], track_);
            evt_.weight_ = nominal_weight;
            if(!upstream_veto) {
              FillTrackHist(trk_hists_[62 + run1a_set_offset], track_);
            }
          }
          // No timing window applied --> 2D fit selection
          FillTrackHist(trk_hists_[65 + run1a_set_offset], track_);
          if(track_->PFront() > 103.34 && track_->PFront() < 104.74
             && track_->TFront() > 640. && track_->TFront() < 1650.) { // FIXME: Timing selection not given for cut-and-count
            FillTrackHist(trk_hists_[66 + run1a_set_offset], track_);
          }
        }
        if(event_id == 0 && run1a_id_no_tdip == 0) {
          FillTrackHist(trk_hists_[67 + run1a_set_offset], track_);
        }

        // no cosmic ID optimized selection
        bool us_cut = true;
        if(track_->upstream_) {
          const float dt = track_->TFront() - track_->upstream_->TFront();
          us_cut = dt < 60.f || dt > 110.f;
        }
        const bool no_csm_opt_id = ((evt_.nde_tracks_ == 1) && (trigger_.FiredAPR() || trigger_.FiredCPR())
                                    && track_->Charge() < 0
                                    && track_->PFront() > 100. && track_->PFront() < 110.
                                    && track_->STBoundary()
                                    && track_->AltPID() > 0.5f && track_->ECluster() > 0.
                                    && us_cut
                                    && track_->CosThetaFront() > 0.525 && track_->CosThetaFront() < 0.649 // start optimized cuts
                                    && track_->RMaxFront() > 482. && track_->RMaxFront() < 642.6
                                    && track_->TFront() > 540. && track_->TFront() < 1650.
                                    && track_->AltTrkQual() > 0.235352
                                    && track_->TrkPID() > 0.078125);
        if(no_csm_opt_id) {
          int cut_opt_offset = 0;
          if(ID & (1 << kCRV)) cut_opt_offset += kCRVVetoOffset;
          FillTrackHist(trk_hists_[73 + cut_opt_offset], track_);
        }
      } // end De selection
    }

    if(evt_.nde_tracks_ != 1) return false; //exactly one positron or electron
    FillEventHist(evt_hists_[1]);
    return true;
  }


  void ConvAna::EndJob() {
    printf("ConvAna::%s\n", __func__);
    cut_flow_.Print();
    run1a_cut_flow_.Print();
  }
}
