//
// Mu2eEvtAna: ntuple and analysis package
// Michael MacKenzie (2024)

#ifndef MU2EEVTANA_MU2EEVTANA_HH
#define MU2EEVTANA_MU2EEVTANA_HH

// standard includes
#include <vector>
#include <utility>

// ROOT includes
#include "Rtypes.h"
#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"

// Mu2e Offline includes

// Mu2e EventNtuple includes
#include "EventNtuple/inc/HitCount.hh"
#include "EventNtuple/inc/TrkCount.hh"
#include "EventNtuple/inc/EventInfo.hh"
#include "EventNtuple/inc/EventInfoMC.hh"
#include "EventNtuple/inc/TrkInfo.hh"
// #include "EventNtuple/inc/TrkInfoMC.hh"
// #include "EventNtuple/inc/TrkSegInfo.hh"
// #include "EventNtuple/inc/LoopHelixInfo.hh"
// #include "EventNtuple/inc/CentralHelixInfo.hh"
// #include "EventNtuple/inc/KinematicLineInfo.hh"
// #include "EventNtuple/inc/SimInfo.hh"
// #include "EventNtuple/inc/EventWeightInfo.hh"
// #include "EventNtuple/inc/TrkStrawHitInfo.hh"
// #include "EventNtuple/inc/TrkStrawHitInfoMC.hh"
#include "EventNtuple/inc/TrkCaloHitInfo.hh"
// #include "EventNtuple/inc/CaloClusterInfoMC.hh"
// #include "EventNtuple/inc/TrkPIDInfo.hh"
// #include "EventNtuple/inc/HelixInfo.hh"
// #include "EventNtuple/inc/InfoStructHelper.hh"
// #include "EventNtuple/inc/CrvInfoHelper.hh"
// #include "EventNtuple/inc/InfoMCStructHelper.hh"
// #include "EventNtuple/inc/RecoQualInfo.hh"
#include "EventNtuple/inc/MVAResultInfo.hh"
// #include "EventNtuple/inc/MCStepInfo.hh"
// #include "EventNtuple/inc/SurfaceStepInfo.hh"
// #include "EventNtuple/inc/MCStepSummaryInfo.hh"

// local includes
#include "Mu2eEvtAna/inc/GlobalConstants.h"
#include "Mu2eEvtAna/inc/Norm_t.hh"
#include "Mu2eEvtAna/inc/Event_t.hh"
#include "Mu2eEvtAna/inc/Track_t.hh"
#include "Mu2eEvtAna/inc/Tracks_t.hh"
#include "Mu2eEvtAna/inc/TrackCaloHit_t.hh"
#include "Mu2eEvtAna/inc/TrackCaloHits_t.hh"
#include "Mu2eEvtAna/inc/EventHist_t.hh"

using namespace mu2e;
namespace Mu2eEvtAna {
  class Mu2eEvtAna {
  public:
    Mu2eEvtAna(int verbose = 0);
    virtual ~Mu2eEvtAna() {};

    int Process(Long64_t nentries = -1, Long64_t first = 0);
    void AddTrack(Tracks_t& trks, Int_t index);
    void AddTrackCollection(TTree* t, Tracks_t& trk, const char* name);
    void BookHistograms(TDirectory* dir);

    void SetInput(TTree* tree) { ntuple_ = tree; }
    void SetName(TString name) { name_ = name; }

    virtual bool ProcessEvent();
    virtual void InitializeEvent();
    virtual void FillOutput();

    virtual int InitializeInput();
    virtual int InitializeOutput();
    virtual void AddOutputBranches(TTree* t);
    virtual void BookEventHist(EventHist_t* Hist);
    virtual void FillEventHist(EventHist_t* Hist);

    virtual TString OutputFileName() { return "evtana_" + name_ + ".root"; }

    TTree* ntuple_; //input ntuple
    EventInfo* evt_info_;
    EventInfoMC* evt_info_mc_;
    std::vector<TrkInfo>* trk_;
    std::vector<TrkCaloHitInfo>* trk_calo_hit_;
    std::vector<MVAResultInfo>* trk_qual_;
    HitCount* hit_cnt_;
    TrkCount* trk_cnt_;

    TFile* fout_; //output file
    TTree* tout_; //output ntuple
    TTree* tnorm_; //output normalization information
    Norm_t norm_; //normalization info
    TDirectory* top_dir_;
    TDirectory* evt_dirs_[kMaxEventHists];
    EventHist_t* evt_hists_[kMaxEventHists];

    TString name_; //name for output file

    Long64_t report_rate_; //reporting rate for events processed
    Int_t verbose_; //verbose level

    // useful fields during processing
    Event_t evt_; //event information
    Tracks_t tracks_; //tracks identified
    Tracks_t electron_; //electron tracks
    Tracks_t muon_; //muon tracks
    Tracks_t proton_; //proton tracks
    Long64_t entry_; //current entry
  };
}

#endif
