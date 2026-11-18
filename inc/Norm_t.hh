//
// Normalization information
// Michael MacKenzie (2025)

#ifndef MU2EEVTANA_NORM_T_HH
#define MU2EEVTANA_NORM_T_HH

// ROOT includes
#include "Rtypes.h"

namespace Mu2eEvtAna {
  struct Norm_t {
    Long64_t ngen_    = 0; //N(generated events)
    Long64_t nntuple_ = 0; //N(events in the input ntuple)
    Long64_t nseen_   = 0; //N(input entries in the processed range, including those failing the preselection)
    Long64_t npresel_ = 0; //N(entries passing the preselection, i.e. actually read and processed; = nseen_ without one)
    Long64_t naccept_ = 0; //N(accepted events)
    Long64_t nneg_    = 0; //N(negative weight events)
  };
}
#endif
