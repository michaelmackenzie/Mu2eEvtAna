//
// Track histograms
// Michael MacKenzie (2025)

#ifndef MU2EEVTANA_TRACKHIST_T_HH
#define MU2EEVTANA_TRACKHIST_T_HH

// ROOT includes
#include "TH1.h"
#include "TH2.h"

namespace Mu2eEvtAna {
  struct TrackHist_t {
    TH1*    fP[2];
    TH1*    fObs; // observable being fit
    TH1*    fPt;
    TH1*    fPCenter[2];
    TH1*    fPExit;
    TH1*    fPST[2];
    TH1*    fPSTDiff;
    TH1*    fPExitDiff;
    TH1*    fT0;
    TH1*    fT0Err;
    TH1*    fD0;
    TH1*    fDP;
    TH1*    fDPCenter;
    TH1*    fChi2NDof;
    TH1*    fFitCons[2];
    TH1*    fFitMomErr; //estimated uncertainty
    TH1*    fTanDip;
    TH1*    fRadius;
    TH1*    fRMax;
    TH1*    fNActive;
    TH1*    fTrkQual;
    TH1*    fClusterE;
    TH1*    fDt;
    TH1*    fEp;
    TH1*    fBestAlg;
    TH1*    fAlgMask;
    TH1*    fTrackID;
    TH1*    fExlTrackID; //flagging what is exclusively cut by a given ID bit

    // Matched CRV cluster info
    TH1*    fCRVDeltaT; //time difference after corrections
    TH1*    fCRVDeltaTCRV; //time difference from stub time
    TH1*    fCRVDeltaTST; //assuming path from ST
    TH1*    fCRVDeltaTCalo; //assuming path from Calo
    TH1*    fCRVDeltaTExtrap; //assuming path from extrapolation
    TH1*    fCRVMinDeltaT; //minimum delta T between ST and Calo paths
    TH2F*    fCRVXZ;
    TH2F*    fCRVYZ;
    TH2F*    fCRVdTZ; //using the corrected time
    TH2F*    fCRVdTZCRV; //using the stub time

    //Matched upstream track info
    TH1*    fUpstreamDt;
    TH1*    fUpstreamDp;

    //MC truth information
    TH1*    fMCPFront; //P(front of tracker)
    TH1*    fMCPStOut; //P(ST exit)
    TH1*    fMCGenE  ; //E(sim particle) at generation
    TH1*    fMCPSig  ; //error / uncertainty
    TH1*    fMCPdg[2];
    TH1*    fMCStrawHits;
    TH1*    fMCGoodHits;
    TH1*    fMCTrajectory;
    TH1*    fMCSimProc;
  };
}
#endif
