//
// Mu2eEvtAna: ntuple and analysis package
// Michael MacKenzie (2025)

#ifndef MU2EEVTANA_MU2EEVTANA_HH
#define MU2EEVTANA_MU2EEVTANA_HH

// standard includes
#include <vector>
#include <utility>
#include <fstream>

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

// Mu2e Offline includes

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
#include "EventNtuple/rooutil/inc/RooUtil.hh"

// local includes
#include "Mu2eEvtAna/inc/GlobalConstants.h"
#include "Mu2eEvtAna/inc/Norm_t.hh"
#include "Mu2eEvtAna/inc/Event_t.hh"
#include "Mu2eEvtAna/inc/Track_t.hh"
#include "Mu2eEvtAna/inc/CRVCluster_t.hh"
#include "Mu2eEvtAna/inc/SimParticle_t.hh"
#include "Mu2eEvtAna/inc/EventHist_t.hh"
#include "Mu2eEvtAna/inc/TrackHist_t.hh"
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
    virtual void InitTrack(Track* track, Track_t& track_par);
    virtual void InitCRVCluster(CrvCoinc* stub, CRVCluster_t& stub_par);
    virtual void FillOutput();

    virtual int InitializeInput();
    virtual int InitializeOutput();
    virtual void AddOutputBranches(TTree* t);
    virtual void BookEventHist(EventHist_t* Hist, const char* Folder);
    virtual void BookTrackHist(TrackHist_t* Hist, const char* Folder);
    virtual void BookCRVHist(CRVHist_t* Hist, const char* Folder);
    virtual void FillEventHist(EventHist_t* Hist);
    virtual void FillTrackHist(TrackHist_t* Hist, Track_t* Track);
    virtual void FillCRVHist(CRVHist_t* Hist, CRVCluster_t* Stub);

    virtual int TrackID(Track_t* track);
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
      case kCRV: return "CRV";
      case kMC: return "MC";
      }
      return Form("Unknown-%i", bit);
    }

    virtual TString OutputFileName() { return "evtana_" + name_ + ".root"; }

    Event*  event_ ; //input TChain wrapper
    TChain* ntuple_; //input ntuple

    TFile* fout_; //output file
    TTree* tout_; //output ntuple
    TTree* tnorm_; //output normalization information
    Norm_t norm_; //normalization info
    TDirectory*  top_dir_;
    TDirectory*  evt_dirs_ [kMaxHists];
    TDirectory*  trk_dirs_ [kMaxHists];
    TDirectory*  crv_dirs_ [kMaxHists];
    EventHist_t* evt_hists_[kMaxHists];
    TrackHist_t* trk_hists_[kMaxHists];
    CRVHist_t*   crv_hists_[kMaxHists];

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

    CRVCluster_t crv_clusters_[kMaxCRVClusters]; //CRV coincidence clusters
    Long64_t entry_; //current entry

    Long64_t        cache_size_   = 200000000U; //200MB cache by default
    Bool_t          load_baskets_ = true;

    // Timer info
    struct Time_t {
      TString name;
      std::chrono::steady_clock::time_point last_time;
      double duration ;
      unsigned count;
      Time_t(TString name = "default") : name(name), last_time(std::chrono::steady_clock::now()), duration(0.), count(0) {}
      void Increment() {
        const auto time_now = std::chrono::steady_clock::now();
        duration += std::chrono::duration_cast<std::chrono::microseconds>(time_now-last_time).count();
        ++count;
        last_time = time_now;
      }
      void SetTime() { last_time = std::chrono::steady_clock::now(); }
      void Reset()   { duration = 0.; count = 0; }
      unsigned Count() { return count; }
      double Time()    { return duration / 1.e6; }
      double AvgTime() { return (count > 0) ? Time() / count : 0.; }
      double AvgRate() { return (duration > 0.) ? count / Time() : 0.; }
    };
    std::map<TString, Time_t> timers_; // for tracking processing time
    Time_t& Timer(TString name) {
      if(!timers_.count(name)) timers_[name] = Time_t(name);
      return timers_[name];
    }
  };
}

#endif
