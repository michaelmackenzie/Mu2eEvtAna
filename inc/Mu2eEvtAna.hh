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
#include "EventNtuple/utils/rooutil/inc/Event.hh"
#include "EventNtuple/utils/rooutil/inc/Track.hh"
#include "EventNtuple/utils/rooutil/inc/RooUtil.hh"

// local includes
#include "Mu2eEvtAna/inc/GlobalConstants.h"
#include "Mu2eEvtAna/inc/Norm_t.hh"
#include "Mu2eEvtAna/inc/Event_t.hh"
#include "Mu2eEvtAna/inc/Track_t.hh"
#include "Mu2eEvtAna/inc/EventHist_t.hh"
#include "Mu2eEvtAna/inc/TrackHist_t.hh"

using namespace mu2e;
namespace Mu2eEvtAna {
  class Mu2eEvtAna {
  public:
    Mu2eEvtAna(int verbose = 0);
    virtual ~Mu2eEvtAna() {};

    int Process(Long64_t nentries = -1, Long64_t first = 0);
    void BookHistograms(TDirectory* dir);

    int  AddFile(TString file_name, Long64_t max_entries = -1);
    void SetInput(TChain* tree) { ntuple_ = tree; }
    void SetName(TString name) { name_ = name; }

    virtual bool ProcessEvent();
    virtual void InitializeEvent();
    virtual void InitTrack(Track* track, Track_t& track_par);
    virtual void FillOutput();

    virtual int InitializeInput();
    virtual int InitializeOutput();
    virtual void AddOutputBranches(TTree* t);
    virtual void BookEventHist(EventHist_t* Hist);
    virtual void BookTrackHist(TrackHist_t* Hist, const char* Folder);
    virtual void FillEventHist(EventHist_t* Hist);
    virtual void FillTrackHist(TrackHist_t* Hist, Track_t* Track);

    virtual TString OutputFileName() { return "evtana_" + name_ + ".root"; }

    Event*  event_ ; //input TChain wrapper
    TChain* ntuple_; //input ntuple

    TFile* fout_; //output file
    TTree* tout_; //output ntuple
    TTree* tnorm_; //output normalization information
    Norm_t norm_; //normalization info
    TDirectory*  top_dir_;
    TDirectory*  evt_dirs_ [kMaxEventHists];
    EventHist_t* evt_hists_[kMaxEventHists];
    TDirectory*  trk_dirs_ [kMaxEventHists];
    TrackHist_t* trk_hists_[kMaxEventHists];

    TString name_; //name for output file

    Long64_t report_rate_; //reporting rate for events processed
    Int_t verbose_; //verbose level

    // useful fields during processing
    Event_t evt_; //event information
    Track_t tracks_  [kMaxTracks]; //tracks identified
    Track_t electron_[kMaxTracks]; //electron tracks
    Track_t muon_    [kMaxTracks]; //muon tracks
    Track_t proton_  [kMaxTracks]; //proton tracks
    Long64_t entry_; //current entry

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
      double AvgTime() { return (count > 0) ? (duration/1.e6) / count : 0.; }
      double AvgRate() { return (duration > 0.) ? count / (duration/1.e6) : 0.; }
    };
    std::map<TString, Time_t> timers_; // for tracking processing time
    Time_t& Timer(TString name) {
      if(!timers_.count(name)) timers_[name] = Time_t(name);
      return timers_[name];
    }
  };
}

#endif
