//
// Line seed histograms
// Michael MacKenzie (2026)

#ifndef MU2EEVTANA_LINESEEDHIST_T_HH
#define MU2EEVTANA_LINESEEDHIST_T_HH

// ROOT includes
#include "TH1.h"
#include "TH2.h"

namespace Mu2eEvtAna {

  struct LineSeedHist_t {
    TH1F* fStatus;
    TH1F* fNHits;
    TH1F* fNStrawHits;
    TH1F* fT0;
    TH1F* fEDep;
    TH1F* fD0;
    TH1F* fPhi0;
    TH1F* fZ0;
    TH1F* fCos;
    TH1F* fA0;
    TH1F* fB0;
    TH1F* fA1;
    TH1F* fB1;
    TH1F* fECalo;
    TH1F* fTCalo;
    TH1F* fHasCalo;

    // Combo hit information, only filled for collections whose hit lists were stored
    TH1F* fHasHits;
    TH1F* fNComboHits;
    TH1F* fHitZ;
    TH1F* fHitR;
    TH1F* fHitTime;
    TH1F* fHitEDep;

    // Matching to a reconstructed track / time cluster (filled by Run1BAna)
    TH1F* fMatchedTrackDt;
    TH1F* fMatchedTrackDD0;
    TH1F* fMatchedTCDt;
  };

}
#endif
