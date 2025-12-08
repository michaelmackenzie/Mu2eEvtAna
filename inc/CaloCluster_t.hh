//
// Calo cluster information
// Michael MacKenzie (2025)

#ifndef MU2EEVTANA_CALOLUSTER_T_HH
#define MU2EEVTANA_CALOLUSTER_T_HH

// ROOT includes
#include "Rtypes.h"

// EventNtuple includes
#include "EventNtuple/inc/CaloClusterInfo.hh"
#include "EventNtuple/inc/CaloHitInfo.hh"

namespace Mu2eEvtAna {
  struct CaloCluster_t {
    mu2e::CaloClusterInfo* cluster_;


    //-------------------------------------------------
    // Accessors


    //-------------------------------------------------
    // Additional functions


    void Reset() {
      cluster_ = nullptr;
    }

    CaloCluster_t() { Reset(); }
  };
}
#endif
