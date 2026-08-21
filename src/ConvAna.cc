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
    hist_sets[ 42] = new hist_info_t("e+: broad window"                 ,  true,  true,  true,  true, false, false,  true, false);
    hist_sets[ 50] = new hist_info_t("e+: no CRV veto"                  ,  true, false, false, false,  true, false,  true, false);
    hist_sets[ 60] = new hist_info_t("e-: Run 1A ID"                    ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 62] = new hist_info_t("e-: Run 1A ID + upstream veto"    ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 65] = new hist_info_t("e-: Run 1A ID, loose time"        ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 66] = new hist_info_t("e-: Run 1A ID, cut-and-count"     ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 70] = new hist_info_t("e-: cut-flow ID"                  ,  true, false, false, false, false, false, false, false);
    hist_sets[ 73] = new hist_info_t("e-: Optimized cut-set"            ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 74] = new hist_info_t("e-: Optimized cut-set"            ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 75] = new hist_info_t("e-: Provided cut-set"             ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 76] = new hist_info_t("e-: Provided cut-set"             ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 77] = new hist_info_t("e-: ID, > 3 ST inters"            ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 78] = new hist_info_t("e-: ID, <= 3 ST inters"           ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 79] = new hist_info_t("e-: Test ID set"                  ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 80] = new hist_info_t("e-: 2D (t,p) selection"           ,  true,  true,  true,  true,  true,  true,  true,  true);

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
  // Fill the standard histograms
  void ConvAna::FillAllHistograms(const int index) {
    if(evt_hists_[index]) FillEventHist(evt_hists_[index]);
    if(trk_hists_[index]) FillTrackHist(trk_hists_[index], track_);
    // if(cls_hists_[index]) FillCaloClusterHist(cls_hists_[index], cluster_);
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
  CutID ConvAna::Run1ATrackID(Track_t* track) {
    if(!track || !track->track_) return 0;
    CutID ID;
    const float pmin = (track->Charge() < 0) ?  97.f :  85.f;
    const float pmax = (track->Charge() < 0) ? 110.f : 110.f;
    if(track->PFront() < pmin || track->PFront() > pmax)          ID.SetBit(kP);        // Loose momentum window
    if(track->OPAInter())                                         ID.SetBit(kRMax);     // Cosmic rejection
    if(std::abs(track->FitPDG()) != 11 || track->PZFront() < 0.f) ID.SetBit(kFitHyp);   // Downstream electron fit
    if(track->TrkQual() > -10. && track->TrkQual() < 0.2)         ID.SetBit(kTrkQual);  // Track quality
    else if(track->TErrFront() > 0.9)                             ID.SetBit(kTrkQual);  // FIXME: Give a label
    else if(track->NActive() < 20)                                ID.SetBit(kTrkQual);  // FIXME: Give a label
    if(track->TFront() < 475. || track->TFront() > 1650.)         ID.SetBit(kT0);       // FIXME: Check selection
    if(track->NSTInter() == 0)                                    ID.SetBit(kD0);       // Consistent with stopping target
    if(track->TanDipFront() < 0.5 || track->TanDipFront() > 0.95) ID.SetBit(kTDip);     // Tan(dip) = pz / pt here
    if(track->TFront() < 475. || track->TFront() > 1650.)         ID.SetBit(kT0Loose);  // For RPC control regions
    if(track->ECluster() <= 0.)                                   ID.SetBit(kClusterE); // Require a cluster
    else if(track->PID() < 0.67f)                                 ID.SetBit(kPID);      // PID requirement

    // CRV rejection
    bool fail_crv = false;
    // Look for an assoicated CRV cluster
    for(int icrv = 0; icrv < evt_.ncrv_clusters_; ++icrv) {
      CRVCluster_t* stub = &crv_clusters_[icrv];
      const float deltat_crv    = track->TFront() - stub->Time();
      if(deltat_crv > 0.f && deltat_crv < 150.f) {
        fail_crv = true;
        break;
      }
    }
    if(fail_crv)                                                  ID.SetBit(kCRV);

    // Upstream track rejection
    // FIXME: Add the standard matching logic

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

    const auto Run1AID = track_->ID(1);
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
                        if(evt_.nde_tracks_ == 1) {
                          run1a_cut_flow_.Increment("nde_tracks");
                          if(track_->PFront() > 103.34 && track_->PFront() < 104.74) {
                            run1a_cut_flow_.Increment("final_mom_region");
                            if(track_->TFront() > 640. && track_->TFront() < 1650.) {
                              run1a_cut_flow_.Increment("final_time_region");
                              if(!Run1AID.Passes() || event_id != 0) {
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
    }
    return false;
  }

  //------------------------------------------------------------------------------------
  // Evaluate the standard selection cut-flow
  bool ConvAna::StandardCutFlow() {
    if(!track_) return false;

    const auto ID = track_->ID(0);
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

    if((ID.CheckBit(kUpstream))) return false;
    cut_flow_.Increment("upstream_veto");
    cut_flow_.Increment("trk_front_seg");
    if((ID.CheckBit(kPID))) return false;
    cut_flow_.Increment("good_trkpid");
    if((ID.CheckBit(kClusterE))) return false;
    cut_flow_.Increment("good_cluster");
    if((ID.CheckBit(kTrkQual))) return false;
    cut_flow_.Increment("good_trkqual");
    if((ID.CheckBit(kFitCon))) return false;
    cut_flow_.Increment("fitcon");
    cut_flow_.Increment("has_hits");
    if((ID.CheckBit(kD0))) return false;
    cut_flow_.Increment("has_st");
    if((ID.CheckBit(kRMax))) return false;
    cut_flow_.Increment("no_opa");
    if((ID.CheckBit(kCRV))) return false;
    cut_flow_.Increment("no_crv_quality");
    cut_flow_.Increment("no_crv_timewindow");
    cut_flow_.Increment("no_crv_veto");
    if((ID.CheckBit(kTDip))) return false;
    cut_flow_.Increment("pz_over_pt");
    if((ID.CheckBit(kCosmicID)) != 0) return false;
    cut_flow_.Increment("cosmic_id");
    if((event_id & (2)) != 0) return false;
    cut_flow_.Increment("trigger");
    if((ID.CheckBit(kP))) return false;
    cut_flow_.Increment("loose_mom_region");
    if(track_->PFront() < 103.5 || track_->PFront() > 105.) return false;
    cut_flow_.Increment("final_mom_region");
    if((ID.CheckBit(kT0Loose))) return false;
    cut_flow_.Increment("loose_time_region");
    if((ID.CheckBit(kT0))) return false;
    cut_flow_.Increment("final_time_region");
    if((event_id & 1) != 0) return false;
    cut_flow_.Increment("de_track_count");

    if(!ID.Passes() || event_id != 0) {
      std::cout << "[ConvAna::" << __func__ << "] "
                << evt_.run_ << ":" << evt_.subrun_ << ":" << evt_.event_
                << " Event passes cuts but fails ID = "
                << std::hex << ID.ID() << std::dec
                << " or event_id = " << event_id
                << std::endl;
    }
    return true;
  }

  //------------------------------------------------------------------------------------
  // Main event-by-event processing
  bool ConvAna::ProcessEvent() {
    ValidateTracks();

    // const float nominal_weight = evt_.weight_;

    FillEventHist(evt_hists_[0]); //all events with well defined inputs

    if(verbose_ > 3) printf("[ConvAna::%s] Event %5i:%6i:%8i\n",
                            __func__, evt_.run_, evt_.subrun_, evt_.event_);

    // all events
    cut_flow_.ResetEvent();
    cut_flow_.Increment("All");
    run1a_cut_flow_.ResetEvent();
    run1a_cut_flow_.Increment("All");
    dev_cut_flow_.ResetEvent();
    dev_cut_flow_.Increment("All");

    // Event-level requirements
    const int event_id = (
                          // 1*(evt_.nde_tracks_ != 1) +
                          2*(!trigger_.FiredAPR() && !trigger_.FiredCPR())
                          );

    // Loop through the track collection
    for(int itrk = 0; itrk < evt_.ntracks_; ++itrk) {
      track_ = &tracks_[itrk];
      dev_cut_flow_.Increment("a_track");
      if(!track_->IsGood()) continue; // if not a properly fit track, skip it
      dev_cut_flow_.Increment("a_converged_track");
      if(track_->PFront() <= 0.) continue;
      dev_cut_flow_.Increment("has_front_seg");
      StandardCutFlow();
      if(Run1ACutFlow()) FillAllHistograms(70);
      if(std::abs(track_->FitPDG()) != 11) continue; // skip muon fits for now due to rooutil bug
      dev_cut_flow_.Increment("is_electron");

      // Downstream electron sets
      if(std::abs(track_->FitPDG()) == 11 && track_->PZFront() > 0.f) {
        dev_cut_flow_.Increment("is_downstream");
        // All high momentum electron tracks
        if(track_->PFront() > 95.f && track_->Charge() < 0) FillAllHistograms(4);

        // Standard IDs
        const auto ID = track_->ID(0);
        const auto Run1AID = track_->ID(1);

        // Cosmic multiplicity/reflection vetos
        bool multi_trk(true), upstream_veto(true);
        for(int i = 0; i < evt_.ntracks_; ++i) {
          if(i == itrk) continue; // skip this track
          const auto alt_trk = &tracks_[i];
          if(!alt_trk->IsGood()) continue; // if not a properly fit track, skip it

          // Check for an upstream partner track
          if(alt_trk->PZFront() < 0.f) {
            const float dt = track_->TFront() - alt_trk->TFront();
            upstream_veto &= dt < 40.f || dt > 110.f; // veto events that are reflection candidates
          }

          // Check for other electrons/positrons in-time with this track
          if(std::abs(alt_trk->FitPDG()) == 11 && alt_trk->PZFront() > 0.f) { // downstream electron/positron track
            const float dt = track_->TFront() - alt_trk->TFront();
            multi_trk &= std::fabs(dt) > 150.; // veto events with tracks coincident with the main track
          }
        }

        // Offset to control regions
        int set_offset(0);
        if(ID.CheckBit(kCRV)) set_offset += kCRVVetoOffset;
        if(ID.CheckBit(kT0))  set_offset += kTimeCutOffset;
        const int id_no_crv = ID.ID(~(1 << kCRV)); // ID without the CRV coincidence cluster considered
        const int id_no_time = ID.ID(~(1 << kT0)); // ID with a looser timing cut
        const int id_no_crv_time = id_no_crv & id_no_time;

        if(event_id == 0 && upstream_veto && multi_trk && id_no_crv_time == 0) {
          FillTrackHist(trk_hists_[10 + set_offset], track_);
          if(track_->Charge() < 0) { // electron selection
            if(track_->PFront() > 100.f && track_->PFront() < 110.f) {
              FillAllHistograms(20 + set_offset);
              if(track_->PFront() > 103.6f && track_->PFront() < 104.9f) {
                FillAllHistograms(21 + set_offset);
              }
              // alternate ID set on top of standard ID (~75% signal eff)
              const bool alt_id = (track_->CosThetaFront() > 0.5 && track_->CosThetaFront() < 0.646805 &&
                                   track_->TFront() > 594.049 &&
                                   track_->RMaxFront() > 482.031 &&
                                   track_->EPFront() < 0.960571 &&
                                   track_->AltPID() > 0.140625);
              if(alt_id) FillAllHistograms(24 + set_offset);
              else       FillAllHistograms(25 + set_offset);
            }
          }
        }

        // Standard selection set
        // Offset to control regions
        int run1a_set_offset(0);
        if(Run1AID.CheckBit(kCRV)) run1a_set_offset += kCRVVetoOffset;
        if(Run1AID.CheckBit(kT0))  run1a_set_offset += kTimeCutOffset;
        const int run1a_id_no_crv = Run1AID.ID(~(1 << kCRV)); // ID without the CRV coincidence cluster considered
        const int run1a_id_no_time = Run1AID.ID(~(1 << kT0)); // ID with a looser timing cut
        const int run1a_id_no_crv_time = run1a_id_no_crv & run1a_id_no_time;

        if(event_id == 0 && run1a_id_no_crv_time == 0 && track_->Charge() < 0) {
          if(track_->TFront() > 640.) { // nominal timing window for 1D fit
            if(evt_.nde_tracks_ == 1) FillAllHistograms(60 + run1a_set_offset);
            if(upstream_veto && multi_trk) {
              FillAllHistograms(62 + run1a_set_offset);

              // No timing window applied --> 2D fit selection
              FillTrackHist(trk_hists_[65 + run1a_set_offset], track_);
              if(track_->PFront() > 103.34 && track_->PFront() < 104.74
                 && track_->TFront() > 640. && track_->TFront() < 1650.) { // FIXME: Timing selection not given for cut-and-count
                FillAllHistograms(66 + run1a_set_offset);
              }
            }
          }
        }

        // no cosmic ID optimized selection
        const bool no_csm_opt_id = ((trigger_.FiredAPR() || trigger_.FiredCPR())
                                    && track_->Charge() < 0
                                    && track_->PFront() > 100. && track_->PFront() < 110.
                                    && track_->STBoundary()
                                    && track_->AltPID() > 0.5f && track_->ECluster() > 0.
                                    && upstream_veto
                                    && multi_trk
                                    && track_->CosThetaFront() > 0.525 && track_->CosThetaFront() < 0.649 // start optimized cuts
                                    && track_->RMaxFront() > 482. && track_->RMaxFront() < 642.6
                                    && track_->TFront() > 540. && track_->TFront() < 1650.
                                    && track_->AltTrkQual() > 0.235352
                                    && track_->TrkPID() > 0.078125);
        if(no_csm_opt_id) {
          int cut_opt_offset = 0;
          if(ID.CheckBit(kCRV)) cut_opt_offset += kCRVVetoOffset;
          FillAllHistograms(73 + cut_opt_offset);
        }

        // Version for other EventNtuple users
        const bool evt_opt_id = ((trigger_.FiredAPR() || trigger_.FiredCPR())
                                 && track_->Charge() < 0
                                 && track_->PFront() > 100. && track_->PFront() < 110.
                                 && track_->STBoundary()
                                 && track_->PID() > 0.5f && track_->ECluster() > 0.
                                 && upstream_veto
                                 && multi_trk
                                 && track_->CosThetaFront() > 0.525 && track_->CosThetaFront() < 0.649 // start optimized cuts
                                 && track_->RMaxFront() > 482. && track_->RMaxFront() < 642.6
                                 && track_->TFront() > 540. && track_->TFront() < 1650.
                                 && track_->TrkQual() > 0.2);
        if(evt_opt_id) {
          int cut_opt_offset = 0;
          if(ID.CheckBit(kCRV)) cut_opt_offset += kCRVVetoOffset;
          FillAllHistograms(74 + cut_opt_offset);
        }

        // Version provided by Natalie from detailed selection optimization
        bool prv_opt_id = true;
        prv_opt_id &= track_->Charge() < 0; if(prv_opt_id) dev_cut_flow_.Increment("charge");
        prv_opt_id &= (trigger_.FiredAPR() || trigger_.FiredCPR()); if(prv_opt_id) dev_cut_flow_.Increment("trigger");
        prv_opt_id &= upstream_veto; if(prv_opt_id) dev_cut_flow_.Increment("upstream_reflection");
        prv_opt_id &= multi_trk; if(prv_opt_id) dev_cut_flow_.Increment("multi_trk");
        prv_opt_id &= track_->PID() > 0.54f && track_->ECluster() > 0.; if(prv_opt_id) dev_cut_flow_.Increment("PID");
        prv_opt_id &= track_->TanDipFront() > 0.575 && track_->TanDipFront() < 0.85; if(prv_opt_id) dev_cut_flow_.Increment("tan_dip");
        prv_opt_id &= track_->STBoundary() > 0; if(prv_opt_id) dev_cut_flow_.Increment("st_boundary");
        prv_opt_id &= track_->NSTInter() > 0; if(prv_opt_id) dev_cut_flow_.Increment("st_inter");
        prv_opt_id &= track_->OPAInter() == 0; if(prv_opt_id) dev_cut_flow_.Increment("opa_inter");
        prv_opt_id &= track_->TrkQual() > 0.155; if(prv_opt_id) dev_cut_flow_.Increment("trkqual");
        prv_opt_id &= track_->NActive() >= 20; if(prv_opt_id) dev_cut_flow_.Increment("nactive");
        prv_opt_id &= track_->TErrMiddle() < 0.85; if(prv_opt_id) dev_cut_flow_.Increment("t0_err");
        if(!Run1AID.CheckBit(kCRV)) {
          if(prv_opt_id) dev_cut_flow_.Increment("crv_veto");
          prv_opt_id &= track_->PFront() > 100. && track_->PFront() < 110.;  if(prv_opt_id) dev_cut_flow_.Increment("momentum");
          prv_opt_id &= track_->TFront() > 475. && track_->TFront() < 1650.; if(prv_opt_id) dev_cut_flow_.Increment("t_475");
          if(prv_opt_id && track_->TFront() > 540.) {
            dev_cut_flow_.Increment("t_540"); // don't actually apply this here
            if(prv_opt_id && track_->TFront() > 640) {
              dev_cut_flow_.Increment("t_640"); // don't actually apply this here
              if(track_->PFront() > 103.34 && track_->PFront() < 104.74) dev_cut_flow_.Increment("sr_momentum");
            }
          }
        }
        if(prv_opt_id) {
          int cut_opt_offset = 0;
          if(Run1AID.CheckBit(kCRV)) cut_opt_offset += kCRVVetoOffset;
          if(track_->TFront() > 540.f) {  // nominal 1D selection
            FillAllHistograms(75 + cut_opt_offset);
            // Test splitting on N(ST foil intersections)
            if(track_->NSTInter() > 3) FillAllHistograms(77 + cut_opt_offset);
            else                       FillAllHistograms(78 + cut_opt_offset);
            if(evt_.nde_tracks_ == 1 && track_->TFront() > 640.) // strict track counting and t0 cut
              FillAllHistograms(76 + cut_opt_offset);

            // Print cosmics that pass all cuts
            if(name_.Contains("cry4a") && !Run1AID.CheckBit(kCRV)) {
              printf("[ConvAna::%s] Event %5i:%6i:%8i passes Run 1A selection\n",
                     __func__, evt_.run_, evt_.subrun_, evt_.event_);
              track_->Print("banner");
              printf("CRV clusters:\n");
              for(int istub = 0; istub < evt_.ncrv_clusters_; ++istub) {
                const auto& stub = crv_clusters_[istub];
                stub.Print((istub == 0) ? "banner" : "");
              }
              printf("Calo clusters:\n");
              for(int icls = 0; icls < evt_.ncalo_clusters_; ++icls) {
                const auto& cls = calo_clusters_[icls];
                cls.Print((icls == 0) ? "banner" : "");
              }
            }
          } // end 540 ns selection
          if(track_->TFront() > 475.f) FillAllHistograms(80 + cut_opt_offset); // nominal 2D selection
        }


        bool test_id = true; // As of 2026-08-21 from Natalie
        test_id &= track_->Charge() < 0;
        test_id &= (trigger_.FiredAPR() || trigger_.FiredCPR());
        test_id &= upstream_veto;
        test_id &= multi_trk;
        test_id &= track_->NSTInter() > 0;
        test_id &= track_->OPAInter() == 0;
        test_id &= track_->D0Front() < 90.;
        test_id &= track_->PID() > 0.525f && track_->ECluster() > 0.;
        test_id &= track_->TanDipFront() > 0.57 && track_->TanDipFront() < 0.9;
        test_id &= track_->TrkQual() > 0.20;
        test_id &= track_->NActive() >= 21;
        test_id &= track_->TErrMiddle() < 0.84;
        test_id &= track_->PFront() > 100. && track_->PFront() < 110.;
        test_id &= track_->TFront() > 540. && track_->TFront() < 1650.;
        if(test_id) {
          int cut_opt_offset = 0;
          if(Run1AID.CheckBit(kCRV)) cut_opt_offset += kCRVVetoOffset;
          FillAllHistograms(79 + cut_opt_offset);
        }

        //------------------------------------
        // Positron selections
        //------------------------------------

        if(track_->PFront() > 75.f && track_->Charge() > 0) FillAllHistograms(5);
        if(track_->Charge() > 0 && id_no_crv_time == 0) {
          if(track_->PFront() > 80.f && track_->PFront() < 100.f) {
            FillAllHistograms(40 + set_offset);
          }
          FillAllHistograms(42 + set_offset); // broad momentum window
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
    dev_cut_flow_.Print();
  }
}
