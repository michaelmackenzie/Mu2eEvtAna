#include "Mu2eEvtAna/inc/Mu2eEvtAna.hh"

using namespace mu2e;
namespace Mu2eEvtAna {

  //------------------------------------------------------------------------------------
  // Constructor
  Mu2eEvtAna::Mu2eEvtAna(int verbose) : ntuple_(nullptr),
                                        fout_(nullptr), tout_(nullptr), tnorm_(nullptr), name_("test"),
                                        report_rate_(1000), verbose_(verbose) {
    //initialize the arrays to null
    for(int ihist = 0; ihist < kMaxHists; ++ihist) {
      evt_hists_[ihist] = nullptr;
      trk_hists_[ihist] = nullptr;
    }

    //Default histogram selections
    evt_hists_[0] = new EventHist_t;
    trk_hists_[0] = new TrackHist_t;

    trk_hists_[1] = new TrackHist_t; // good electron tracks
    trk_hists_[2] = new TrackHist_t; // electron tracks + track ID
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
            if(verbose_ > -1) printf("\r%s: Loaded %i files of %i with %llu entries", __func__, ifile, nfiles, ntuple_->GetEntries());
            break;
          }
        }
        filelist.close();
      } else {
        printf("%s: Error! Unable to read input file list %s\n", __func__, file_name.Data());
        return 1;
      }
    }
    if(verbose_ > -1) printf("\n");
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
  void Mu2eEvtAna::BookEventHist(EventHist_t* Hist, const char* Folder) {
    if(!Hist) {
      throw std::runtime_error("Attempting to book histograms in a null EventHist_t\n");
    }

    Hist->fInstLumi         = new TH1F("inst_lumi"        ,Form("%s: POT",Folder), 300,  0.0, 1.5e8);
    Hist->fInstLumiApr      = new TH1F("inst_lumi_apr"    ,Form("%s: POT",Folder), 300,  0.0, 1.5e8);
    Hist->fInstLumiCpr      = new TH1F("inst_lumi_cpr"    ,Form("%s: POT",Folder), 300,  0.0, 1.5e8);
    Hist->fInstLumiAprCpr   = new TH1F("inst_lumi_apr_cpr",Form("%s: POT",Folder), 300,  0.0, 1.5e8);
    Hist->fEventWeight[0]   = new TH1F("event_weight"     ,Form("%s: Event weight",Folder), 100, 0., 2.);
    Hist->fEventWeight[1]   = new TH1F("event_weight_log" ,Form("%s: log10(Event weight)",Folder), 100, -20., 1.);
    Hist->fNAprTracks       = new TH1F("nTracksApr"       ,Form("%s: nTracksApr",Folder), 50, 0.0, 50.0);
    Hist->fNCprTracks       = new TH1F("nTracksCpr"       ,Form("%s: nTracksCpr",Folder), 50, 0.0, 50.0);
    Hist->fNTracks          = new TH1F("nTracks"          ,Form("%s: nTracks",Folder), 50, 0.0, 50.0);
    Hist->fNUeTracks        = new TH1F("nUeTracks"        ,Form("%s: nUeTracks",Folder), 50, 0.0, 50.0);
    Hist->fNDmuTracks       = new TH1F("nDmuTracks"       ,Form("%s: nUmuTracks",Folder), 50, 0.0, 50.0);
    Hist->fNUmuTracks       = new TH1F("nUmuTracks"       ,Form("%s: nUmuTracks",Folder), 50, 0.0, 50.0);
    Hist->fNGoodTrks        = new TH1D("ngoodtrks"        ,Form("%s: N(good tracks)"    , Folder),  10, 0,  10);
    Hist->fNIDTrks          = new TH1D("nidtrks"          ,Form("%s: N(accepted tracks)", Folder),  10, 0,  10);
    Hist->fNAprHelices      = new TH1F("nHelicesApr"      ,Form("%s: nHelicesApr",Folder), 50, 0.0, 50.0);
    Hist->fNCprHelices      = new TH1F("nHelicesCpr"      ,Form("%s: nHelicesCpr",Folder), 50, 0.0, 50.0);
    Hist->fNHelices         = new TH1F("nHelices"         ,Form("%s: nHelices",Folder), 50, 0.0, 50.0);
    Hist->fNCRVClusters     = new TH1F("nCRVClusters"     ,Form("%s: N(CRV clusters)",Folder), 50, 0.0, 50.0);
    Hist->fNGoodCRVClusters = new TH1F("nGoodCRVClusters" ,Form("%s: N(Good CRV clusters)",Folder), 50, 0.0, 50.0);
    Hist->fNonCRVVetoID     = new TH1F("nonCRVVetoID"     ,Form("%s: Non-CRV Veto ID",Folder), 30, 0., 30.);
    Hist->fNDigis           = new TH1D("ndigis"           ,Form("%s: N(digis)"   ,Folder), 200, 0, 200);
    Hist->fNClusters        = new TH1D("nclusters"        ,Form("%s: N(clusters)",Folder),  10, 0,  10);
    Hist->fPrimaryCode      = new TH1F("primary_code"     ,Form("%s: Primary process code",Folder), 200, -0.5, 199.5);
    Hist->fPrimaryGenE      = new TH1F("primary_gene"     ,Form("%s: Primary gen energy",Folder), 500, 50., 150.);

  }

  //------------------------------------------------------------------------------------
  // Initialize the histograms for a track selection
  void Mu2eEvtAna::BookTrackHist(TrackHist_t* Hist, const char* Folder) {
    if(!Hist) {
      throw std::runtime_error("Attempting to book histograms in a null TrackHist_t\n");
    }
    Hist->fP[0]        = new TH1F("p"           ,Form("%s: Track momentum"                       ,Folder),  300,    0.,  150.);
    Hist->fP[1]        = new TH1F("p_2"         ,Form("%s: Track momentum"                       ,Folder),  400,   80.,  120.);
    Hist->fObs         = new TH1F("obs"         ,Form("%s: Track momentum"                       ,Folder),  150,   80.,  110.); // fit histogram
    Hist->fPt          = new TH1F("pt"          ,Form("%s: track transverse momentum"            ,Folder),  300,    0.,  300.);
    Hist->fPCenter[0]  = new TH1F("pCenter"     ,Form("%s: track momentum at tracker center"     ,Folder),  600, -300.,  300.);
    Hist->fPCenter[1]  = new TH1F("pCenter_2"   ,Form("%s: track momentum at tracker center"     ,Folder),  600,   80.,  110.);
    Hist->fPExit       = new TH1F("pExit"       ,Form("%s: track momentum at tracker exit"       ,Folder),  300,    0.,  300.);
    Hist->fPST[0]      = new TH1F("pST"         ,Form("%s: track momentum at ST exit"            ,Folder),  300,    0.,  300.);
    Hist->fPST[1]      = new TH1F("pST_2"       ,Form("%s: track momentum at ST exit"            ,Folder),  300,   80.,  110.);
    Hist->fPSTDiff     = new TH1F("pST_diff"    ,Form("%s: track p(ST) - p(Front)"               ,Folder),  400,   -1.,    9.);
    Hist->fPExitDiff   = new TH1F("pExit_diff"  ,Form("%s: track p(Front) - p(Exit)"             ,Folder),  400,   -1.,    4.);
    Hist->fT0          = new TH1F("t0"          ,Form("%s: track t_{0}"                          ,Folder),  400,    0., 2000.);
    Hist->fT0Err       = new TH1F("t0err"       ,Form("%s: track t_{0} uncertainty"              ,Folder),  100,    0.,   20.);
    Hist->fD0          = new TH1F("d0"          ,Form("%s: track d0"                             ,Folder),  200, -200.,  200.);
    Hist->fDP          = new TH1F("dP"          ,Form("%s: track p_reco - p_mc"                  ,Folder),  400,  -20.,   20.);
    Hist->fChi2NDof    = new TH1F("chi2NDof"    ,Form("%s: track chi2/ndof"                      ,Folder),  200,    0.,   10.);
    Hist->fFitCons[0]  = new TH1F("fitCons"     ,Form("%s: track p(chi2,ndof)"                   ,Folder),  200,    0.,    1.);
    Hist->fFitCons[1]  = new TH1F("fitCons_log" ,Form("%s: track log10(p(chi2,ndof))"            ,Folder),  200,   -6.,    0.);
    Hist->fFitMomErr   = new TH1F("fitMomErr"   ,Form("%s: track momentum uncertainty"           ,Folder),  200,    0.,    5.);
    Hist->fTanDip      = new TH1F("tanDip"      ,Form("%s: track tanDip"                         ,Folder),  200,    0.,    2.);
    Hist->fRadius      = new TH1F("radius"      ,Form("%s: track radius"                         ,Folder), 1000,    0., 1000.);
    Hist->fRMax        = new TH1F("rMax"        ,Form("%s: track rMax"                           ,Folder), 2000,    0., 2000.);
    Hist->fNActive     = new TH1F("nActive"     ,Form("%s: nHits used in fit"                    ,Folder),  150,    0.,  150.);
    Hist->fTrkQual     = new TH1F("trkQual"     ,Form("%s: track MVA score"                      ,Folder),  200,   -1.,    1.);
    Hist->fClusterE    = new TH1F("clusterE"    ,Form("%s: track's cluster energy"               ,Folder),  600,    0.,  300.);
    Hist->fDt          = new TH1F("dt"          ,Form("%s: track - cluster time"                 ,Folder),  200,  -10.,   10.);
    Hist->fEp          = new TH1F("ep"          ,Form("%s: cluster E / track P"                  ,Folder),  200,    0.,    2.);
    Hist->fBestAlg     = new TH1F("bestAlg"     ,Form("%s: Best fit algorithm"                   ,Folder),   10,    0.,   10.);
    Hist->fAlgMask     = new TH1F("algMask"     ,Form("%s: Algorithm mask"                       ,Folder),   10,    0.,   10.);
    Hist->fTrackID     = new TH1F("track_id"    ,Form("%s: Track ID bits"                        ,Folder),   33,    0.,   33.);
    Hist->fExlTrackID  = new TH1F("track_exl_id",Form("%s: Track ID bits for exclusive rejection",Folder),   32,    0.,   32.);

    // Initialize bin labels for the Track ID histograms
    Hist->fTrackID   ->GetXaxis()->SetBinLabel(1, "Passed");
    Hist->fExlTrackID->GetXaxis()->SetBinLabel(1, "Passed");
    int current_bin = 2;
    for(int bit = 0; bit <= Hist->fTrackID->GetNbinsX(); ++bit) {
      TString name(TrackIDBitName(bit));
      if(name.BeginsWith("Unknown")) continue;
      Hist->fTrackID   ->GetXaxis()->SetBinLabel(current_bin, name.Data());
      Hist->fExlTrackID->GetXaxis()->SetBinLabel(current_bin, name.Data());
      ++current_bin;
    }

    // Matched CRV cluster info
    Hist->fCRVDeltaT       = new TH1F("crv_deltat"       ,Form("%s: Track t_{0} - CRV t_{0}"                     ,Folder), 200, -250., 250.);
    Hist->fCRVDeltaTCRV    = new TH1F("crv_deltat_crv"   ,Form("%s: Track t_{0} - CRV t_{0}"                     ,Folder), 200, -250., 500.);
    Hist->fCRVDeltaTST     = new TH1F("crv_deltat_st"    ,Form("%s: Track t_{0} - CRV through ST t_{0}"          ,Folder), 200, -250., 250.);
    Hist->fCRVDeltaTCalo   = new TH1F("crv_deltat_calo"  ,Form("%s: Track t_{0} - CRV through Calo t_{0}"        ,Folder), 200, -250., 250.);
    Hist->fCRVDeltaTExtrap = new TH1F("crv_deltat_extrap",Form("%s: Track t_{0} - CRV through Extrapolated t_{0}",Folder), 200, -250., 250.);
    Hist->fCRVMinDeltaT    = new TH1F("crv_min_deltat"   ,Form("%s: Min #Deltat_{0} for ST and Calo paths"       ,Folder), 200, -250., 250.);
    Hist->fCRVXZ           = new TH2F("crv_x_vs_z"       ,Form("%s: CRV X vs Z"                                  ,Folder), 250, -5000, 20000, 200, -10000, 10000);
    Hist->fCRVYZ           = new TH2F("crv_y_vs_z"       ,Form("%s: CRV Y vs Z"                                  ,Folder), 250, -5000, 20000, 200,      0,  4000);
    Hist->fCRVdTZ          = new TH2F("crv_dt_vs_z"      ,Form("%s: #Deltat vs CRV Z"                            ,Folder), 250, -5000, 20000, 200,   -200,   200);
    Hist->fCRVdTZCRV       = new TH2F("crv_dtcrv_vs_z"   ,Form("%s: #Deltat(CRV) vs CRV Z"                       ,Folder), 250, -5000, 20000, 200,   -200,   200);

    // Matched upstream track info
    Hist->fUpstreamDt      = new TH1F("us_dt",Form("%s: Upstream track #Deltat_{0}"          ,Folder), 150,  -50.,  250.);
    Hist->fUpstreamDp      = new TH1F("us_dp",Form("%s: Upstream track #Deltap"              ,Folder), 100,  -10.,    5.);

    // MC truth
    Hist->fMCPFront      = new TH1F("MC_PFront",Form("%s: MC track P(tracker front)"         ,Folder), 600, -300.,  300.);
    Hist->fMCPStOut      = new TH1F("MC_PSTOut",Form("%s: MC track P(ST exit)"               ,Folder), 600, -300.,  300.);
    Hist->fMCGenE        = new TH1F("MC_GenE"  ,Form("%s: MC generated energy"               ,Folder), 400,    0.,  200.);
    Hist->fMCPSig        = new TH1F("MC_PSig"  ,Form("%s: MC P(front) error / uncertainty"   ,Folder), 400,  -20.,   20.);
    Hist->fMCPdg[0]      = new TH1F("MC_PDG_0",Form("%s: MC Particle PDG code"               ,Folder),  40,  -20.,   20.);
    Hist->fMCPdg[1]      = new TH1F("MC_PDG_1",Form("%s: MC Particle |PDG code|"             ,Folder), 220,    0., 2200.);
    Hist->fMCStrawHits   = new TH1F("MC_strawhits",Form("%s: MC Particle N(straw hits)"      ,Folder), 100,    0.,  100.);
    Hist->fMCGoodHits    = new TH1F("MC_goodhits",Form("%s: MC Particle N(good hits)"        ,Folder), 100,    0.,  100.);
    Hist->fMCTrajectory  = new TH1F("MC_trajectory",Form("%s: MC track p_{z} trajectory"     ,Folder),   3,  -1.5,   1.5);
    Hist->fMCSimProc     = new TH1F("MC_simProc",Form("%s: MC Sim process code"              ,Folder), 200,  -0.5, 199.5);
  }

  //------------------------------------------------------------------------------------
  // Initialize the histogram sets
  void Mu2eEvtAna::BookHistograms(TDirectory* dir) {

    for(int ihist = 0; ihist < kMaxHists; ++ihist) {
      if(evt_hists_[ihist]) {
        const char* folder = Form("evt_%i", ihist);
        auto subdir = dir->mkdir(folder);
        subdir->cd();
        BookEventHist(evt_hists_[ihist], folder);
        dir->cd();
        evt_dirs_[ihist] = subdir;
      }
      if(trk_hists_[ihist]) {
        const char* folder = Form("trk_%i", ihist);
        auto subdir = dir->mkdir(folder);
        subdir->cd();
        BookTrackHist(trk_hists_[ihist], folder);
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
    const double Weight(evt_.weight_);

    Hist->fNTracks   ->Fill(evt_.ntracks_   , Weight);
    Hist->fNGoodTrks ->Fill(evt_.ngoodtrks_ , Weight);
    Hist->fNIDTrks   ->Fill(evt_.ntrks_id_  , Weight);
    Hist->fNDigis    ->Fill(evt_.ndigis_    , Weight);
    Hist->fNClusters ->Fill(evt_.nclusters_ , Weight);

    Hist->fInstLumi->Fill(evt_.inst_lum_, Weight);
    if (evt_.passed_apr_) { Hist->fInstLumiApr->Fill(evt_.inst_lum_, Weight); }
    if (evt_.passed_cpr_) { Hist->fInstLumiCpr->Fill(evt_.inst_lum_, Weight); }
    if (evt_.passed_apr_ || evt_.passed_cpr_) { Hist->fInstLumiAprCpr->Fill(evt_.inst_lum_, Weight); }
    Hist->fEventWeight[0]->Fill(Weight);
    Hist->fEventWeight[1]->Fill((Weight > 0.f) ? std::log10(Weight) : -1000.);
    Hist->fNAprTracks->Fill(evt_.napr_tracks_, Weight);
    Hist->fNCprTracks->Fill(evt_.ncpr_tracks_, Weight);
    Hist->fNUeTracks->Fill(evt_.nue_tracks_, Weight);
    Hist->fNDmuTracks->Fill(evt_.ndmu_tracks_, Weight);
    Hist->fNUmuTracks->Fill(evt_.numu_tracks_, Weight);
    Hist->fNAprHelices->Fill(evt_.napr_helices_, Weight);
    Hist->fNCprHelices->Fill(evt_.ncpr_helices_, Weight);
    Hist->fNHelices->Fill(evt_.noffline_helices_, Weight);
    Hist->fNCRVClusters->Fill(evt_.ncrv_clusters_, Weight);
    Hist->fNGoodCRVClusters->Fill(evt_.ngood_crvclusters_, Weight);
    // // for(int bit = 0; bit < 30; ++bit) {
    // //   if((evt_.fNonCRVVetoID & (1 << bit)) != 0) Hist->fNonCRVVetoID->Fill(bit, Weight);
    // // }
    // // Hist->fPrimaryCode->Fill((fPrimary) ? fPrimary->CreationCode() : 0 , Weight);
    // // Hist->fPrimaryGenE->Fill((fPrimary) ? fPrimary->fStartMom.E()  : 0., Weight);

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
    const float Weight(evt_.weight_);
    Hist->fP[0] ->Fill(Track->PFront(), Weight);
    Hist->fP[1] ->Fill(Track->PFront(), Weight);
    Hist->fObs->Fill(Track->Obs(), Weight); //FIXME: Should this be Ana module dependent?
    Hist->fPt->Fill(Track->PTFront(), Weight);
    Hist->fPCenter[0]->Fill(Track->PMiddle()*Track->Charge(), Weight);
    Hist->fPCenter[1]->Fill(Track->PMiddle(), Weight);
    Hist->fPExit->Fill(Track->PBack(), Weight);
    Hist->fPST[0]->Fill(Track->PSTBack(), Weight);
    Hist->fPST[1]->Fill(Track->PSTBack(), Weight);
    Hist->fPSTDiff->Fill(Track->PSTBack()-Track->PFront(), Weight);
    Hist->fPExitDiff->Fill(Track->PFront()-Track->PBack(), Weight);
    Hist->fT0->Fill(Track->TFront(), Weight);
    Hist->fT0Err->Fill(Track->TErrFront(), Weight);
    Hist->fD0->Fill(Track->D0Front(), Weight);
    Hist->fDP->Fill(Track->MCDeltaPFront(), Weight);
    Hist->fChi2NDof->Fill(Track->Chi2Dof(), Weight);
    Hist->fFitCons[0]->Fill(Track->FitCon(), Weight);
    Hist->fFitCons[1]->Fill(std::log10(std::max(1.e-10f, Track->FitCon())), Weight);
    Hist->fFitMomErr->Fill(Track->MomErrFront(), Weight);
    Hist->fTanDip->Fill(Track->TanDipFront(), Weight);
    Hist->fRadius->Fill(Track->RadiusFront(), Weight);
    Hist->fRMax->Fill(Track->RMaxFront(), Weight);
    Hist->fNActive->Fill(Track->NActive(), Weight);
    Hist->fTrkQual->Fill(Track->TrkQual(), Weight);
    Hist->fClusterE->Fill(Track->ECluster(), Weight);
    Hist->fDt->Fill(Track->Dt(), Weight);
    Hist->fEp->Fill(Track->EPFront(), Weight);
    // Hist->fBestAlg->Fill(TrkPar->fTrack->BestAlg(), Weight);
    // Hist->fAlgMask->Fill(TrkPar->fTrack->AlgMask(), Weight);
    const int ID = Track->ID(0);
    if(ID == 0) {
      Hist->fTrackID   ->Fill("Passed", Weight);
      Hist->fExlTrackID->Fill("Passed", Weight);
    } else {
      for(int bit = 0; bit < 30; ++bit) {
        if((ID & (1 << bit)) != 0) {
          TString name = TrackIDBitName(bit);
          Hist->fTrackID->Fill(name.Data(), Weight);
          if((ID & ~(1 << bit)) == 0) Hist->fExlTrackID->Fill(name.Data(), Weight);
        }
      }
    }

    // auto crv_stub_par = TrkPar->fCRVStubPar;
    // if(crv_stub_par) {
    //   auto crv_stub = crv_stub_par->fCluster;
    //   const float deltat_st   = TrkPar->fTrack->fT0 - crv_stub_par->fApproxTimeSTToFront;
    //   const float deltat_calo = TrkPar->fTrack->fT0 - crv_stub_par->fApproxTimeCaloToFront;
    //   const float deltat_extrap = TrkPar->fTrack->fT0 - crv_stub_par->fApproxTimeExtrapToFront;
    //   const float min_deltat = (std::fabs(deltat_st) < std::fabs(deltat_calo)) ? deltat_st : deltat_calo;
    //   Hist->fCRVDeltaT->Fill(TrkPar->fTrack->fT0 - crv_stub_par->fCorrTime, Weight);
    //   Hist->fCRVDeltaTCRV->Fill(TrkPar->fTrack->fT0 - crv_stub_par->fTime, Weight);
    //   Hist->fCRVDeltaTST->Fill(deltat_st, Weight);
    //   Hist->fCRVDeltaTCalo->Fill(deltat_calo, Weight);
    //   Hist->fCRVDeltaTExtrap->Fill(deltat_extrap, Weight);
    //   Hist->fCRVMinDeltaT->Fill(min_deltat, Weight);
    //   const float x = crv_stub->Position()->X();
    //   const float y = crv_stub->Position()->Y();
    //   const float z = crv_stub->Position()->Z();
    //   Hist->fCRVXZ->Fill(z,x,Weight);
    //   Hist->fCRVYZ->Fill(z,y,Weight);
    //   const float dt = TrkPar->fTrack->fT0 - crv_stub_par->fCorrTime;
    //   const float dt_crv = TrkPar->fTrack->fT0 - crv_stub_par->fTime;
    //   Hist->fCRVdTZ->Fill(z,dt,Weight);
    //   Hist->fCRVdTZCRV->Fill(z,dt_crv,Weight);
    // }

    // auto us_trk = TrkPar->fUpstreamTrack;
    // if(us_trk) {
    //   Hist->fUpstreamDt->Fill(TrkPar->fTrack->fT0 - us_trk->fT0,Weight);
    //   Hist->fUpstreamDp->Fill(TrkPar->fTrack->fP  - us_trk->fP ,Weight);
    // }

    Hist->fMCPFront->Fill(Track->MCPFront(), Weight);
    Hist->fMCPStOut->Fill(Track->MCPSTBack(), Weight);
    Hist->fMCGenE->Fill(Track->MCGenE(), Weight);
    Hist->fMCPSig->Fill((Track->MomErrFront() > 0.) ? Track->MCDeltaPFront() / Track->MomErrFront() : -999., Weight);
    Hist->fMCPdg[0]->Fill(Track->MCPDG(), Weight);
    Hist->fMCPdg[1]->Fill(std::abs(Track->MCPDG()), Weight);
    Hist->fMCStrawHits->Fill(Track->MCHits(), Weight);
    Hist->fMCGoodHits->Fill(Track->MCActive(), Weight);
    Hist->fMCTrajectory->Fill((Track->MCPZFront() < 0.) ? -1 : 1, Weight);
    Hist->fMCSimProc->Fill(Track->MCProcess(), Weight);
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

    // Retrieve global event info
    InitEvent(evt_);

    // Initialize the sim particle info (if available)
    // FIXME: Currently need to do this within the track collection

    // Add track information

    // Loop through the tracks
    auto tracks = event_->GetTracks();
    evt_.ntracks_ = tracks.size();
    for(int itrk = 0; itrk < evt_.ntracks_; ++itrk) {
      // First initialize any new sim particle info
      if(tracks[itrk].trkmcsim) {
        for(size_t index = 0; index < tracks[itrk].trkmcsim->size(); ++index) {
          auto& simp = tracks[itrk].trkmcsim->at(index);
          // Check if this sim particle has already been seen
          bool found = false;
          for(int isimp = 0; isimp < evt_.nsimps_; ++isimp) {
            if(simp.id == simps_[isimp].id_) {
              found = true;
              break;
            }
          }
          if(!found) {
            if(evt_.nsimps_ + 1 >= kMaxSimps) throw std::runtime_error("Exceeded the maximum number of sim particles!\n");
            if(verbose_ > 2) printf("%s: Adding Sim particle: ID = %2i, PDG = %5i, Start Code = %3i\n", __func__, simp.id, simp.pdg, simp.startCode);
            auto& simp_t = simps_[evt_.nsimps_];
            ++evt_.nsimps_;
            simp_t.Initialize(&simp);
          }
        }
      }

      // Initialize the track info
      InitTrack(&tracks[itrk], tracks_[itrk]);
      if(tracks_[itrk].IsGood()) ++evt_.ngoodtrks_;
      if(verbose_ > 1) tracks_[itrk].Print((itrk == 0) ? "banner" : "");
    }
    // Reset remaining tracks
    for(int itrk = evt_.ntracks_; itrk < kMaxTracks; ++itrk) {
      tracks_[itrk].Reset();
    }

    if(verbose_ > 4) {
      printf("Mu2eEvtAna::%s: Printing event information:\n", __func__);
      printf(" N(tracks) = %2i\n", evt_.ntracks_);
    }
  }

  //------------------------------------------------------------------------------------
  // Initialize event information
  void Mu2eEvtAna::InitEvent(Event_t& evt) {
    evt.Reset();
    if(!event_) return;
    if(event_->evtinfo) {
      evt.run_      = event_->evtinfo->run   ;
      evt.subrun_   = event_->evtinfo->subrun;
      evt.event_    = event_->evtinfo->event ;
    }
    if(event_->evtinfomc) {
      evt.inst_lum_ = event_->evtinfomc->nprotons;
      evt.pb_time_  = event_->evtinfomc->pbtime  ;
    }
    if(event_->hitcount) {
      evt.ndigis_   = event_->hitcount->nsd;
    }
  }

  //------------------------------------------------------------------------------------
  // Initialize track information
  void Mu2eEvtAna::InitTrack(Track* track, Track_t& trk_par) {
    trk_par.Reset();
    trk_par.track_ = track;
    if(!track) return;
    trk_par.SetObs(trk_par.PFront(), 0); // default to the track momentum as the key observable
    trk_par.SetObs(trk_par.TFront(), 1); // default to the track time as the secondary key observable

    // Set initial track ID info
    trk_par.SetID(TrackID(&trk_par), 0);
  }

  //------------------------------------------------------------------------------------
  // Track selection
  int Mu2eEvtAna::TrackID(Track_t* track) {
    if(!track || !track->track_) return 0;
    int ID(0);
    if(track->PFront() < 70. || track->PFront() > 130.)          ID += 1 << kP;
    if(track->RMaxFront() < 430. || track->RMaxFront() > 650.)   ID += 1 << kRMax;
    if(track->TrkQual() > -10. && track->TrkQual() < 0.2)        ID += 1 << kTrkQual;
    if(track->TFront() < 600. || track->TFront() > 1650.)        ID += 1 << kT0;
    if(track->FitCon() < 1.e-5)                                  ID += 1 << kFitCon;
    if(track->ECluster() <= 0.)                                  ID += 1 << kClusterE; //require a cluster
    const float d0_sign = track->D0Front() * track->Charge();
    if(d0_sign < -100. || d0_sign > 60.)                         ID += 1 << kD0; //consistent with ST
    if(track->TanDipFront() < 0.5 || track->TanDipFront() > 1.5) ID += 1 << kTDip;
    if(track->TFront() < 500. || track->TFront() > 1650.)        ID += 1 << kT0Loose; //for control regions

    // // upstream track rejection
    // auto us_trk = trkpar.fUpstreamTrack;
    // if(us_trk) {
    //   // minimal selection on the upstream track quality
    //   if(us_trk->fFitCons > 1.e-5 &&
    //      (us_trk->TrkQual() < -10. || us_trk->TrkQual() > 0.01)) ID += 1 << kUpstream;
    // }

    // PID rejection FIXME: Add MVA identification
    if(track->EPFront() < 0.65)                                   ID += 1 << kPID;

    // // CRV rejection
    // if(trkpar.fCRVStubPar) {
    //   auto stub_par = trkpar.fCRVStubPar;
    //   const float deltat_calo = fTrack->fT0 - stub_par->fApproxTimeCaloToFront; //consider both trajectory hypotheses
    //   const float deltat_st   = fTrack->fT0 - stub_par->fApproxTimeSTToFront  ;
    //   const float deltatcrv = fTrack->fT0 - stub_par->fTime; //effective for upstream misreconstruction
    //   const float min_extrap_dt(-50.f), max_extrap_dt(60.f);
    //   if((deltat_st   > min_extrap_dt && deltat_st   < max_extrap_dt) ||
    //      (deltat_calo > min_extrap_dt && deltat_calo < max_extrap_dt) ||
    //      (deltatcrv > -25.f && deltatcrv < 0.f))               ID += 1 << kCRV;
    // }

    // // MC selection to match MC cuts applied to the Mock Data Samples
    // bool mc_selection = true;
    // if(fPrimary
    //    && fPrimary->CreationCode() != mu2e::ProcessCode::mu2eFlatPhoton  // FIXME: don't cut on RMC photons for now
    //    && trkpar.fTrack->Charge() < 0) {  // Only apply this to negative tracks
    //   mc_selection &= trkpar.fGenE <= 0. || trkpar.fGenE > 95.f; // only included DIO with E > 95 MeV
    // }

    // //MC-truth cut, don't include pileup in background-specific histogramming
    // if(fDataType == kBackground) {
    //   bool find_match(false); //FIXME: Add MC relation field to TSimParticle to check track particle relation
    //   switch(fBackgroundType) {
    //   case kInternalRPC:
    //   case kExternalRPC:
    //   case kInternalRMC:
    //   case kExternalRMC:
    //   case kDIO: find_match = true; break;
    //   default: break;
    //   }
    //   if(find_match) {
    //     bool found(false);
    //     if(fSimpBlock) {
    //       for(int isim = 0; isim < fSimpBlock->NParticles(); ++isim) {
    //         const auto sim = fSimpBlock->Particle(isim);
    //         found |= track->SimID() == int(sim->GetUniqueID());
    //         if(found) break;
    //       }
    //     }
    //     mc_selection &= found;
    //   }
    // }
    // if(!mc_selection)                                          ID += 1 << kMC;

    return ID;
  }

  //------------------------------------------------------------------------------------
  // Main event-by-event processing
  bool Mu2eEvtAna::ProcessEvent() {
    timers_["Event"].Increment();
    if(verbose_ > 4) {
      printf("Mu2eEvtAna::%s: Printing event information:\n", __func__);
    }
    FillEventHist(evt_hists_[0]); //all events with well defined inputs
    for(int itrk = 0; itrk < evt_.ntracks_; ++itrk) { // all tracks
      FillTrackHist(trk_hists_[0], &tracks_[itrk]);
      bool trk_id = tracks_[itrk].IsGood() && tracks_[itrk].FitPDG() == 11;
      trk_id &= tracks_[itrk].PZFront() > 0.f;
      if(trk_id) {
        FillTrackHist(trk_hists_[1], &tracks_[itrk]);
        const int ID = tracks_[itrk].ID();
        if(ID == 0) FillTrackHist(trk_hists_[2], &tracks_[itrk]);
      }
    }

    return false; // default to not writing output trees
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
        printf("Mu2eEvtAna::%s: Processing event %7lld (entry %8lld, event %6i/%7i/%7i): N(accept) = %7lld (%6.2f%%) (%.0f Hz)\n", __func__, nseen, entry,
               event_->evtinfo->run,event_->evtinfo->subrun,event_->evtinfo->event, naccepted, naccepted*100./((nseen <= 0) ? 1 : nseen),
               Timer("Event").AvgRate());
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
