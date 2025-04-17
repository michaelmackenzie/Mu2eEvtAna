//
// Event information
// Michael MacKenzie (2024)

#ifndef MU2EEVTANA_EVENT_T_HH
#define MU2EEVTANA_EVENT_T_HH

// ROOT includes
#include "Rtypes.h"

namespace Mu2eEvtAna {
  struct Event_t {
    Double_t weight_ = 1.;
    Int_t ndigis_ = 0;
    Int_t ntrks_ = 0;
    Int_t nclusters_ = 0;

    Int_t nelectrons_ = 0;
    Int_t npositrons_ = 0;
    Int_t nmu_minus_  = 0;
    Int_t nmu_plus_   = 0;
  };
}
#endif
