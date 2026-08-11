//
// Mu2eEvtAna: ntuple and analysis package
// Michael MacKenzie (2025)

#ifndef MU2EEVTANA_MU2EEVTANA_HH
#define MU2EEVTANA_MU2EEVTANA_HH

// standard includes
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <bitset>

// ROOT includes
#include "Rtypes.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#endif

// Mu2e Offline includes
#include "Offline/Mu2eUtilities/inc/StopWatch.hh"

// Mu2e EventNtuple includes
#include "EventNtuple/inc/HitCount.hh"
#include "EventNtuple/inc/TrkCount.hh"
#include "EventNtuple/inc/EventInfo.hh"
#include "EventNtuple/inc/EventInfoMC.hh"
#include "EventNtuple/inc/TrkInfo.hh"
#include "EventNtuple/inc/TrkCaloHitInfo.hh"
#include "EventNtuple/inc/MVAResultInfo.hh"

// Mu2e EventNtuple RooUtil includes
#include "EventNtuple/rooutil/inc/Event.hh"
#include "EventNtuple/rooutil/inc/Track.hh"
#include "EventNtuple/rooutil/inc/Trigger.hh"
#include "EventNtuple/rooutil/inc/RooUtil.hh"

// local includes
#include "Mu2eEvtAna/inc/GlobalConstants.h"
#include "Mu2eEvtAna/inc/CutID.hh"
#include "Mu2eEvtAna/inc/Norm_t.hh"
#include "Mu2eEvtAna/inc/Event_t.hh"
#include "Mu2eEvtAna/inc/Track_t.hh"
#include "Mu2eEvtAna/inc/CRVCluster_t.hh"
#include "Mu2eEvtAna/inc/CaloCluster_t.hh"
#include "Mu2eEvtAna/inc/SimParticle_t.hh"
#include "Mu2eEvtAna/inc/Trigger_t.hh"
#include "Mu2eEvtAna/inc/EventHist_t.hh"
#include "Mu2eEvtAna/inc/TrackHist_t.hh"
#include "Mu2eEvtAna/inc/CaloClusterHist_t.hh"
#include "Mu2eEvtAna/inc/CRVHist_t.hh"

using namespace mu2e;
namespace Mu2eEvtAna {
  class Mu2eEvtAna {
  public:

    Mu2eEvtAna(int verbose = 0);
    virtual ~Mu2eEvtAna() {};

    int Process(Long64_t nentries = -1, Long64_t first = 0);
    virtual void InitHistSelections();
    virtual void BookHistograms(TDirectory* dir);

    int  AddFile(TString file_name, Long64_t max_entries = -1, Long64_t first_entry = 0);
    void SetInput(TChain* tree) { ntuple_ = tree; }
    void SetName(TString name) { name_ = name; }

    virtual bool ProcessEvent();
    virtual void InitializeEvent();
    virtual void InitEvent(Event_t& evt);
    virtual void InitTrack(const rooutil::Track* track, Track_t& track_par);
    virtual void InitCaloCluster(const rooutil::CaloCluster* calo, CaloCluster_t& cls_par);
    virtual void InitCRVCluster(const rooutil::CrvCoinc* stub, CRVCluster_t& stub_par);
    virtual void FillOutput();
    virtual void EndJob();

    virtual int InitializeInput();
    virtual int InitializeOutput();
    virtual void AddOutputBranches(TTree* t);
    virtual void BookEventHist(EventHist_t* Hist, const char* Folder);
    virtual void BookTrackHist(TrackHist_t* Hist, const char* Folder);
    virtual void BookCaloClusterHist(CaloClusterHist_t* Hist, const char* Folder);
    virtual void BookCRVHist(CRVHist_t* Hist, const char* Folder);
    virtual void FillEventHist(EventHist_t* Hist);
    virtual void FillTrackHist(TrackHist_t* Hist, Track_t* Track);
    virtual void FillCaloClusterHist(CaloClusterHist_t* Hist, CaloCluster_t* Cluster);
    virtual void FillCRVHist(CRVHist_t* Hist, CRVCluster_t* Stub);

    virtual CutID TrackID(Track_t* track);
    static TString TrackIDBitName(const int bit) {
      switch(bit) {
      case kP: return "P";
      case kRMax: return "R(max)";
      case kTrkQual: return "TrkQual";
      case kT0: return "T_case 0}";
      case kFitCon: return "p(#chi^2)";
      case kClusterE: return "E_case CL}";
      case kD0: return "D_case 0}";
      case kTDip: return "tan(dip)";
      case kT0Loose: return "T_case 0} (loose)";
      case kUpstream: return "Upstream";
      case kPID: return "PID";
      case kFitHyp: return "Fit hyp.";
      case kCosmicID: return "Cosmic ID";
      case kCRV: return "CRV";
      case kMC: return "MC";
      }
      return Form("Unknown-%i", bit);
    }

    virtual TString OutputFileName() { return "EvtAna." + name_ + ".root"; }

    void ValidateTracks() {
      const auto& tracks = event_->GetTracks();
      if(evt_.ntracks_ != int(tracks.size())) {
        throw std::runtime_error(Form("Track size changed from %i to %zu!", evt_.ntracks_, tracks.size()));
      }
      for(int itrk = 0; itrk < evt_.ntracks_; ++itrk) {
        if(&tracks[itrk] != tracks_[itrk].track_) {
          throw std::runtime_error(Form("Track %i pointer changed from %p to %p!", itrk, (void*) tracks_[itrk].track_, (void*) &tracks[itrk]));
        }
      }
    }

    rooutil::Event*  event_ ; //input TChain wrapper
    TChain* ntuple_; //input ntuple
    TTree*  tree_  = nullptr; //current tree in the TChain

    TFile* fout_; //output file
    TTree* tout_; //output ntuple
    TTree* tnorm_; //output normalization information
    Norm_t norm_; //normalization info
    TDirectory*        top_dir_;
    TDirectory*        hist_dir_;
    TDirectory*        data_dir_;
    TDirectory*        evt_dirs_  [kMaxHists];
    TDirectory*        trk_dirs_  [kMaxHists];
    TDirectory*        cls_dirs_  [kMaxHists];
    TDirectory*        crv_dirs_  [kMaxHists];
    EventHist_t*       evt_hists_ [kMaxHists];
    TrackHist_t*       trk_hists_ [kMaxHists];
    CaloClusterHist_t* cls_hists_ [kMaxHists];
    CRVHist_t*         crv_hists_ [kMaxHists];

    TString name_; //name for output file

    Long64_t report_rate_; //reporting rate for events processed
    Int_t verbose_; //verbose level

    // useful fields during processing
    Event_t evt_; //event information
    SimParticle_t simps_ [kMaxSimps ]; //Relevant sim particles
    Track_t tracks_      [kMaxTracks]; //tracks identified
    Track_t* de_tracks_  [kMaxTracks]; //Downstream electrons
    Track_t* ue_tracks_  [kMaxTracks]; //Upstream electrons
    Track_t* dmu_tracks_ [kMaxTracks]; //Downstream muons
    Track_t* umu_tracks_ [kMaxTracks]; //Upstream muons
    Trigger_t trigger_; //trigger information

    CRVCluster_t crv_clusters_[kMaxCRVClusters]; //CRV coincidence clusters
    CaloCluster_t calo_clusters_[kMaxCaloClusters]; //Calo clusters
    Long64_t entry_; //current entry

    Long64_t        cache_size_   = 200000000U; //200MB cache by default
    Bool_t          load_baskets_ = true;
    Bool_t          use_xrootd_   = true;

    Long64_t        tree_entries_ = 0;

    mu2e::StopWatch* watch_; // track processing times
  };
}

#endif
