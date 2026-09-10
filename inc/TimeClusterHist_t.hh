//
// Time cluster histograms
// Michael MacKenzie (2026)

#ifndef MU2EEVTANA_TIMECLUSTERHIST_T_HH
#define MU2EEVTANA_TIMECLUSTERHIST_T_HH

// ROOT includes
#include "TH1.h"
#include "TH2.h"

namespace Mu2eEvtAna {

  struct TimeClusterHist_t {
    TH1F* fNHits;
    TH1F* fNStrawHits;
    TH1F* fT0;
    TH1F* fX;
    TH1F* fY;
    TH1F* fZ;
    TH1F* fR;
    TH1F* fECalo;
    TH1F* fTCalo;
    TH1F* fHasCalo;
  };

}
#endif
