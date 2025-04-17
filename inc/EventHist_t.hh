//
// Event histograms
// Michael MacKenzie (2024)

#ifndef MU2EEVTANA_EVENTHIST_T_HH
#define MU2EEVTANA_EVENTHIST_T_HH

// ROOT includes
#include "TH1.h"
#include "TH2.h"

namespace Mu2eEvtAna {
  struct EventHist_t {
    TH1* hNTrks;
    TH1* hNDigis;
    TH1* hNClusters;
  };
}
#endif
