//
// PrimaryAssociation: whether an MC particle is the event primary or one of its descendants
// Michael MacKenzie (2026)

#ifndef MU2EEVTANA_PRIMARYASSOCIATION_HH
#define MU2EEVTANA_PRIMARYASSOCIATION_HH

// Mu2e EventNtuple includes
#include "EventNtuple/inc/SimInfo.hh"

namespace Mu2eEvtAna {
  // A SimInfo's prirel is its relationship to the event's primary particles (MCRelationship(sim, primary)): same for the
  // primary itself, daughter/udaughter for its (ultimate) descendants. Siblings and mothers of the primary, and particles
  // with no relationship (e.g. pileup), are not associated. prirel is filled for track sims (trkmcsim), and for calo
  // sims (calomcsim) only by EventNtuple versions that fill it there; otherwise it is none.
  inline bool IsPrimaryRelated(const mu2e::SimInfo* sim) {
    if(!sim) return false;
    const auto rel = sim->prirel.relationship();
    return rel == mu2e::MCRelationship::same || rel == mu2e::MCRelationship::daughter || rel == mu2e::MCRelationship::udaughter;
  }
}

#endif
