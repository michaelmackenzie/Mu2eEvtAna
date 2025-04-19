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
    Double_t weight_;
    Int_t ndigis_   ;
    Int_t ngoodtrks_;
    Int_t ntrks_id_ ;
    Int_t nclusters_;

    // Track counters
    Int_t ntracks_   ;
    Int_t nelectrons_;
    Int_t nmuons_    ;
    Int_t nprotons_  ;

    // MC info
    Int_t nsimps_    ;

    void Reset() {
      weight_    = 1.;
      ndigis_    = 0;
      ngoodtrks_ = 0;
      ntrks_id_  = 0;
      nclusters_ = 0;
      ntracks_    = 0;
      nelectrons_ = 0;
      nmuons_     = 0;
      nprotons_   = 0;
      nsimps_     = 0;
    }

    Event_t() { Reset(); }
  };
}
#endif
