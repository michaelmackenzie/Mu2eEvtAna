//
// BNVAna: histogramming for the search for baryon number violating muon capture into light dark
// states (the BNVAna analysis): e+- tracks above the DIO endpoint, and in-time e+e- track pairs.
// Derives from ConvAna to share its MVA evaluation, event weighting, and histogram booking.
// Michael MacKenzie (2026)

#ifndef MU2EEVTANA_BNVANA_HH
#define MU2EEVTANA_BNVANA_HH

// local includes
#include "Mu2eEvtAna/inc/ConvAna.hh"

using namespace mu2e;
namespace Mu2eEvtAna {
  class BNVAna : public ConvAna {
  public:
    BNVAna(int verbose = 0);
    ~BNVAna() {};

    enum { kBNVTrackID = 3 }; // Track_t ID index of the BNV track ID (0/1: standard e-/e+, 2: ConvAna Run 1A)

    void InitHistSelections();
    bool ProcessEvent();
    void EndJob();

    // Standard track IDs plus the BNV track ID
    void SetTrackIDs(Track_t* track);
    CutID BNVTrackID(Track_t* track);

    TString OutputFileName() { return "BNVAna." + name_ + ".root"; }

    CutFlow bnv_cut_flow_; // BNV selection
  };
}

#endif
