//
// Event information
// Michael MacKenzie (2025)

#ifndef MU2EEVTANA_EVENT_T_HH
#define MU2EEVTANA_EVENT_T_HH

// ROOT includes
#include "Rtypes.h"

namespace Mu2eEvtAna {
  struct Event_t {
    // Event info
    Double_t weight_ = 1.;
    Int_t ndigis_    = 0;
    Int_t ngoodtrks_ = 0;
    Int_t ntrks_id_  = 0;
    Int_t nclusters_ = 0;

    // Track counters
    Int_t ntracks_    = 0;
    Int_t nelectrons_ = 0;
    Int_t nmuons_     = 0;
    Int_t nprotons_   = 0;
  };
}
#endif
