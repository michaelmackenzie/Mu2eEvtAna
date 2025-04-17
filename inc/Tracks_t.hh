//
// Track collection
// Michael MacKenzie (2024)

#ifndef MU2EEVTANA_TRACKS_T_HH
#define MU2EEVTANA_TRACKS_T_HH

// ROOT includes
#include "Rtypes.h"

// local includes
#include "Mu2eEvtAna/inc/GlobalConstants.h"
#include "Mu2eEvtAna/inc/Track_t.hh"
#include "Mu2eEvtAna/inc/TrackCaloHit_t.hh"
#include "Mu2eEvtAna/inc/TrackCaloHits_t.hh"

namespace Mu2eEvtAna {
  struct Tracks_t {
    Int_t    size_ = 0                   ; // size of the collection
    Int_t    keep_ = true                ;
    Int_t    index_          [kMaxTracks]; // index in the original track collection (for matching)
    Int_t    status_         [kMaxTracks]; // track fit information
    Int_t    pdg_            [kMaxTracks];
    Int_t    nhits_          [kMaxTracks]; // track hit information
    Int_t    nactive_        [kMaxTracks];
    Int_t    ndouble_        [kMaxTracks];
    Int_t    nplanes_        [kMaxTracks];
    Int_t    nnull_          [kMaxTracks];
    Int_t    first_hit_      [kMaxTracks];
    Int_t    last_hit_       [kMaxTracks];
    Int_t    ndof_           [kMaxTracks]; // fit quality information
    Float_t  chisq_          [kMaxTracks];
    Float_t  fitcon_         [kMaxTracks];
    Int_t    reco_id_        [kMaxTracks]; // quality and selection IDs
    Float_t  trk_qual_       [kMaxTracks];
    Float_t  p_corr_         [kMaxTracks]; // corrections
    Float_t  id_weight_      [kMaxTracks];
    Int_t    gen_trk_index_  [kMaxTracks]; // pointers to associated indices in other collections
    Int_t    seg_first_index_[kMaxTracks];
    Int_t    trk_calo_index_ [kMaxTracks];

    TrackCaloHits_t calo_hits_; //list of associated track-calo hits

    void reset() {size_ = 0; calo_hits_.reset();}

    Track_t track(int index) {
      if(index >= size_)
        throw std::runtime_error(Form("Attempting to access a track with an index (%i) beyond the collection size (%i)",
                                      index, size_));
      if(size_ >= kMaxTracks)
        throw std::runtime_error(Form("Track collection is too large with collection size = %i", size_));
      if(index < 0) return Track_t();

      // Initialize a track with the parameters from the collection
      Track_t t;
      t.index_           = index_          [index]; // index in the original track collection (for matching)
      t.status_          = status_         [index]; // track fit information
      t.pdg_             = pdg_            [index];
      t.nhits_           = nhits_          [index]; // track hit information
      t.nactive_         = nactive_        [index];
      t.ndouble_         = ndouble_        [index];
      t.nplanes_         = nplanes_        [index];
      t.nnull_           = nnull_          [index];
      t.first_hit_       = first_hit_      [index];
      t.last_hit_        = last_hit_       [index];
      t.ndof_            = ndof_           [index]; // fit quality information
      t.chisq_           = chisq_          [index];
      t.fitcon_          = fitcon_         [index];
      t.reco_id_         = reco_id_        [index]; // quality and selection IDs
      t.trk_qual_        = trk_qual_       [index];
      t.p_corr_          = p_corr_         [index]; // corrections
      t.id_weight_       = id_weight_      [index];
      t.gen_trk_index_   = gen_trk_index_  [index]; // pointers to associated indices in other collections
      t.seg_first_index_ = seg_first_index_[index];
      t.trk_calo_index_  = trk_calo_index_ [index];
      t.calo_hit_ = calo_hits_.trkcalohit(t.trk_calo_index_);
      return t;
    }
  };
}
#endif
