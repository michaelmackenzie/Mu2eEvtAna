//
// Track+calorimeter hit
// Michael MacKenzie (2024)

#ifndef MU2EEVTANA_TRACKCALOHIT_T_HH
#define MU2EEVTANA_TRACKCALOHIT_T_HH

// ROOT includes
#include "Rtypes.h"

// local includes
#include "Mu2eEvtAna/inc/GlobalConstants.h"

namespace Mu2eEvtAna {
  struct TrackCaloHit_t {
    Int_t    index_          = -1; // original collection index (for matching)
    Int_t    disk_id_        = -1; // cluster information
    Float_t  cl_t_           = -1;
    Float_t  cl_t_err_       = -1;
    Float_t  cl_e_           = -1;
    Float_t  cl_e_err_       = -1;
    Float_t  cl_size_        = -1;
    Int_t    cl_ncr_         = -1;
    Float_t  x_              = -1; // position of closest approach (POCA)
    Float_t  y_              = -1;
    Float_t  z_              = -1;
    Float_t  px_             = -1; // track momentum at POCA
    Float_t  py_             = -1;
    Float_t  pz_             = -1;
    Float_t  doca_           = -1; // distance of closest approach
    Float_t  dt_             = -1; // track vs. cluster info
  };
}
#endif
