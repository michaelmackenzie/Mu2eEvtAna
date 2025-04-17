//
// Track+calorimeter hit collection
// Michael MacKenzie (2024)

#ifndef MU2EEVTANA_TRACKCALOHITS_T_HH
#define MU2EEVTANA_TRACKCALOHITS_T_HH

// ROOT includes
#include "Rtypes.h"

// local includes
#include "Mu2eEvtAna/inc/GlobalConstants.h"
#include "Mu2eEvtAna/inc/TrackCaloHit_t.hh"

namespace Mu2eEvtAna {
  struct TrackCaloHits_t {
    Int_t    size_                       ; // size of the collection
    Int_t    index_          [kMaxTracks]; // original collection index (for matching)
    Int_t    disk_id_        [kMaxTracks]; // cluster information
    Float_t  cl_t_           [kMaxTracks];
    Float_t  cl_t_err_       [kMaxTracks];
    Float_t  cl_e_           [kMaxTracks];
    Float_t  cl_e_err_       [kMaxTracks];
    Float_t  cl_size_        [kMaxTracks];
    Int_t    cl_ncr_         [kMaxTracks];
    Float_t  x_              [kMaxTracks]; // position of closest approach (POCA)
    Float_t  y_              [kMaxTracks];
    Float_t  z_              [kMaxTracks];
    Float_t  px_             [kMaxTracks]; // track momentum at POCA
    Float_t  py_             [kMaxTracks];
    Float_t  pz_             [kMaxTracks];
    Float_t  doca_           [kMaxTracks]; // distance of closest approach
    Float_t  dt_             [kMaxTracks]; // track vs. cluster info

    void reset() {size_ = 0;}

    TrackCaloHit_t trkcalohit(int index) {
      if(index >= size_)
        throw std::runtime_error(Form("Attempting to access a TrkCaloHit with an index (%i) beyond the collection size (%i)",
                                      index, size_));
      if(size_ >= kMaxTracks)
        throw std::runtime_error(Form("TrkCaloHit collection is too large with collection size = %i", size_));
      if(index < 0) return TrackCaloHit_t();
      TrackCaloHit_t t;
      t.index_      = index_     [index]; // original collection index (for matching)
      t.disk_id_    = disk_id_   [index]; // cluster information
      t.cl_t_       = cl_t_      [index];
      t.cl_t_err_   = cl_t_err_  [index];
      t.cl_e_       = cl_e_      [index];
      t.cl_e_err_   = cl_e_err_  [index];
      t.cl_size_    = cl_size_   [index];
      t.cl_ncr_     = cl_ncr_    [index];
      t.x_          = x_         [index]; // position of closest approach (POCA)
      t.y_          = y_         [index];
      t.z_          = z_         [index];
      t.px_         = px_        [index]; // track momentum at POCA
      t.py_         = py_        [index];
      t.pz_         = pz_        [index];
      t.doca_       = doca_      [index]; // distance of closest approach
      t.dt_         = dt_        [index]; // track vs. cluster info

      return t;
    }
  };
}
#endif
