//
// Event histograms
// Michael MacKenzie (2025)

#ifndef MU2EEVTANA_EVENTHIST_T_HH
#define MU2EEVTANA_EVENTHIST_T_HH

// ROOT includes
#include "TH1.h"
#include "TH2.h"

namespace Mu2eEvtAna {
  struct EventHist_t {
    TH1*    fInstLumi;
    TH1*    fInstLumiApr; // only filled if apr triggered event
    TH1*    fInstLumiCpr; // only filled if cpr triggered event
    TH1*    fInstLumiAprCpr; // filled if either apr or cpr triggered event
    TH1*    fEventWeight[2];
    TH1*    fNAprHelices;
    TH1*    fNCprHelices;
    TH1*    fNHelices;
    TH1*    fNMatchedHelices;
    TH1*    fNAprTracks;
    TH1*    fNCprTracks;
    TH1*    fNTracks;
    TH1*    fNUeTracks;
    TH1*    fNDmuTracks;
    TH1*    fNUmuTracks;
    TH1*    fNCRVClusters;
    TH1*    fNGoodCRVClusters;
    TH1*    fNonCRVVetoID;
    TH1*    fNGoodTrks;
    TH1*    fNIDTrks;
    TH1*    fNDigis;
    TH1*    fNClusters;

    TH1*    fTrackerHits;
    TH1*    fCaloHits;

    TH1*    fTrigBits[2];
    TH1*    fTrigPaths[2];
    TH2*    fTrigOverlap[2];
    TH1*    fNTriggerable;
    TH1D*   fTriggered;

    // Primary process info
    TH1*    fPrimaryCode;
    TH1*    fPrimaryType;
    TH1*    fPrimaryGenE;
    TH1*    fRMCEnergy;
  };
}
#endif
