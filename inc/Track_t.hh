//
// A signle Track
// Michael MacKenzie (2024)

#ifndef MU2EEVTANA_TRACK_T_HH
#define MU2EEVTANA_TRACK_T_HH

// ROOT includes
#include "Rtypes.h"

// local includes
#include "Mu2eEvtAna/inc/GlobalConstants.h"
#include "Mu2eEvtAna/inc/TrackCaloHit_t.hh"

namespace Mu2eEvtAna {
  struct Track_t {
    Int_t    index_           = -1; // index in the original track collection (for matching)
    Int_t    status_          = -1; // track fit information
    Int_t    pdg_             = -1;
    Int_t    nhits_           = -1; // track hit information
    Int_t    nactive_         = -1;
    Int_t    ndouble_         = -1;
    Int_t    nplanes_         = -1;
    Int_t    nnull_           = -1;
    Int_t    first_hit_       = -1;
    Int_t    last_hit_        = -1;
    Int_t    ndof_            = -1; // fit quality information
    Float_t  chisq_           = -1;
    Float_t  fitcon_          = -1;
    Int_t    reco_id_         = -1; // quality and selection IDs
    Float_t  trk_qual_        = -1;
    Float_t  p_corr_          = -1; // corrections
    Float_t  id_weight_       = -1;
    Int_t    gen_trk_index_   = -1; // pointers to associated indices in other collections
    Int_t    seg_first_index_ = -1;
    Int_t    trk_calo_index_  = -1;

    TrackCaloHit_t calo_hit_; //list of associated track-calo hits
  };
}
#endif
