//
// Mu2eEvtAna systematic histograms
// Michael MacKenzie (2025)

#ifndef MU2EEVTANA_SYSHIST_T_HH
#define MU2EEVTANA_SYSHIST_T_HH

// ROOT includes
#include "TH1.h"
#include "TH2.h"

// local includes
#include "Mu2eEvtAna/inc/GlobalConstants.h"

namespace Mu2eEvtAna {

  struct SysHist_t {
    TH1* fObs        [kMaxSystematics];
    TH1* fDeltaObs   [kMaxSystematics];
    TH1* fWeight     [kMaxSystematics];
    TH1* fDeltaWeight[kMaxSystematics];

    SysHist_t() {
      for(int i = 0; i < kMaxSystematics; ++i) {
        fObs        [i] = nullptr;
        fDeltaObs   [i] = nullptr;
        fWeight     [i] = nullptr;
        fDeltaWeight[i] = nullptr;
      }
    }
  };
}
#endif
