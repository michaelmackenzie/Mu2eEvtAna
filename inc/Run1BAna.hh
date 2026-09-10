//
// Run1BAna: Run1B (straight-line, field-off) analysis ntupling/histogramming
// Michael MacKenzie (2026)

#ifndef MU2EEVTANA_RUN1BANA_HH
#define MU2EEVTANA_RUN1BANA_HH

// standard includes
#include <vector>

// ROOT includes

// Mu2e Offline includes

// Mu2e EventNtuple includes

// local includes
#include "Mu2eEvtAna/inc/Mu2eEvtAna.hh"
#include "Mu2eEvtAna/inc/TimeCluster_t.hh"
#include "Mu2eEvtAna/inc/LineSeed_t.hh"
#include "Mu2eEvtAna/inc/TimeClusterHist_t.hh"
#include "Mu2eEvtAna/inc/LineSeedHist_t.hh"
#include "Mu2eEvtAna/inc/CutFlow.hh"

using namespace mu2e;
namespace Mu2eEvtAna {
  // Run1BAna reads the TimeClusterCollection/LineSeed branches (read automatically, however many
  // collections and whatever they are named -- see rooutil::Event::TimeClusterCollectionNames()/
  // LineSeedCollectionNames()) and handles tracks fit as KinematicLine (straight-line) segments
  // instead of LoopHelix. It changes nothing about the base Mu2eEvtAna/ConvAna default behavior.
  class Run1BAna : public Mu2eEvtAna {
  public:
    Run1BAna(int verbose = 0);
    ~Run1BAna() {};

    void InitHistSelections();
    void BookHistograms(TDirectory* dir);
    void BookTimeClusterHist(TimeClusterHist_t* Hist, const char* Folder);
    void BookLineSeedHist(LineSeedHist_t* Hist, const char* Folder);
    void FillTimeClusterHist(TimeClusterHist_t* Hist, const TimeCluster_t* Cluster);
    void FillLineSeedHist(LineSeedHist_t* Hist, const LineSeed_t* Seed);

    bool ProcessEvent();
    void InitializeEvent();

    // Straight-line-appropriate track selection: replaces the LoopHelix RMax/tan(dip) windows
    // with RTrackerMax()/RTrackerMin() and drops the tan(dip) requirement (meaningless for a
    // field-off track). Cut values below are placeholders -- initial, wide windows to be tuned
    // once real Run1B statistics/efficiency studies are available (see the developer plan).
    CutID TrackID(Track_t* track);

    int InitializeInput();
    int InitializeOutput();

    void EndJob();

    TString OutputFileName() { return "Run1BAna." + name_ + ".root"; }

    // Best-matched line seed for each track/time cluster, found in InitializeEvent()
    int BestLineSeed(const Track_t* track, int icoll) const;
    int BestTimeCluster(const LineSeed_t* seed, int icoll) const;

    // Populate each calo_clusters_[icls]'s best-matched line/line seed/time cluster/CRV cluster
    // (CaloCluster_t::line_/line_seed_/time_cluster_/crv_cluster_). Must run after tracks_,
    // crv_clusters_, time_clusters_, and line_seeds_ are all populated for the event -- i.e. at
    // the end of InitializeEvent(), the same way the base class' upstream-track matching runs
    // only after every track has been through InitTrack().
    void MatchCaloClusters();

    // Discovered collection names (from the input ntuple; may differ file to file)
    std::vector<TString> tc_names_;
    std::vector<TString> ls_names_;

    // Per-collection data, indexed to match tc_names_/ls_names_
    std::vector<std::vector<TimeCluster_t>> time_clusters_; // [collection][cluster]
    std::vector<std::vector<LineSeed_t>>    line_seeds_;    // [collection][seed]

    std::vector<TDirectory*>        tc_dirs_;
    std::vector<TDirectory*>        ls_dirs_;
    std::vector<TimeClusterHist_t*> tc_hists_;
    std::vector<LineSeedHist_t*>    ls_hists_;

    CutFlow cut_flow_; // straight-line track selection cut flow
  };
}

#endif
