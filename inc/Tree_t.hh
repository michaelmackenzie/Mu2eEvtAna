//
// Slim tree information
// Michael MacKenzie (2026)

#ifndef MU2EEVTANA_TREE_T_HH
#define MU2EEVTANA_TREE_T_HH

#include "TTree.h"

#include <vector>

namespace Mu2eEvtAna {

  struct Tree_t {

    // Event info
    int   run;
    int   subrun;
    int   event;
    float weight;
    float train;

    // Track information
    float trk_p;
    float trk_t0;
    float trk_d0;
    float trk_tandip;
    float trk_cos;
    float trk_cluster;
    float trk_ep;
    float trk_dt;
    float trk_rmax;
    float trk_nactive;
    float trk_active_ratio;
    float trk_null_ratio;
    float trk_active_mat_ratio;
    float trk_fitcon;
    float trk_logfitcon;
    float trk_momerr;
    float trk_t0err;
    float trk_tzslope;
    float trk_tzslope_sig;
    float trk_tzslope_ratio;
    float trk_pexit_diff;
    float trk_qual;
    float trk_pid;
    float trk_trkonlypid;
    float trk_altqual;
    float trk_altpid;
    float trk_cosmicid;
    float trk_charge;
    float trk_min_crv_time;
    int   trk_nst_down; // intersections
    int   trk_nst_up;
    int   trk_stboundary;
    int   trk_nipa;
    int   trk_opa;
    int   trk_tsda;
    float trk_mc_dp; // MC info
    float trk_mc_pdg;

    // CRV information
    float crv_z;
    float crv_deltat;
    float crv_npulses;
    float crv_npe;
    float crv_npepp;

    //----------------------------------------------------------
    // Track quality info

    float trkqual_nactive;
    float trkqual_activehitsfraction;
    float trkqual_nullhitsfraction;
    float trkqual_activematsitesfraction;
    float trkqual_fitcons;
    float trkqual_momerr;
    float trkqual_t0err;

    //----------------------------------------------------------
    // Specific for TReflectionAna

    // Upstream track info
    float refl_p_us;
    float refl_pt_us;
    float refl_cos_us;
    float refl_d0_us;
    float refl_r_us;
    float refl_rmax_us;
    float refl_rmin_us;
    float refl_ep_us;
    int refl_nhits_us;
    int refl_q_us;
    int refl_stboundary_us;
    int refl_stinters_us;
    int refl_triggered_us;

    // Downstream track info
    float refl_p_ds;
    float refl_pt_ds;
    float refl_cos_ds;
    float refl_d0_ds;
    float refl_r_ds;
    float refl_rmax_ds;
    float refl_rmin_ds;
    float refl_ep_ds;
    int refl_nhits_ds;
    int refl_q_ds;
    int refl_stboundary_ds;
    int refl_stinters_ds;
    int refl_triggered_ds;

    // Comparison info
    float refl_dp;
    float refl_dt0;
    float refl_dtfront;

    // Photon CNN info
    std::vector<float> crystal_e;
    std::vector<float> crystal_t;

    // TTree
    TTree* tree = nullptr;

    Tree_t() { Reset(); }

    void Reset() {

      run = 0;
      subrun = 0;
      event = 0;
      weight = 1.f;
      train = 0.f;

      trk_p = 0.f;
      trk_t0 = 0.f;
      trk_d0 = 0.f;
      trk_tandip = 0.f;
      trk_cluster = 0.f;
      trk_ep = 0.f;
      trk_dt = 0.f;
      trk_rmax = 0.f;
      trk_nactive = 0.f;
      trk_active_ratio = 0.f;
      trk_null_ratio = 0.f;
      trk_active_mat_ratio = 0.f;
      trk_fitcon = 0.f;
      trk_logfitcon = 0.f;
      trk_momerr = 0.f;
      trk_t0err = 0.f;
      trk_tzslope = 0.f;
      trk_tzslope_sig = 0.f;
      trk_tzslope_ratio = 0.f;
      trk_mc_dp = 0.f;
      trk_mc_pdg = 0.f;

      crv_z = 0.f;
      crv_deltat = 0.f;
      crv_npulses = 0.f;
      crv_npe = 0.f;
      crv_npepp = 0.f;

      trkqual_nactive = 0.f;
      trkqual_activehitsfraction = 0.f;
      trkqual_nullhitsfraction = 0.f;
      trkqual_activematsitesfraction = 0.f;
      trkqual_fitcons = 0.f;
      trkqual_momerr = 0.f;
      trkqual_t0err = 0.f;

      refl_p_us = 0.f;
      refl_pt_us = 0.f;
      refl_cos_us = 0.f;
      refl_d0_us = 0.f;
      refl_r_us = 0.f;
      refl_rmax_us = 0.f;
      refl_rmin_us = 0.f;
      refl_ep_us = 0.f;
      refl_nhits_us = 0;
      refl_q_us = 0;
      refl_stboundary_us = 0;
      refl_stinters_us = 0;
      refl_triggered_us = 0;
      refl_p_ds = 0.f;
      refl_pt_ds = 0.f;
      refl_cos_ds = 0.f;
      refl_d0_ds = 0.f;
      refl_r_ds = 0.f;
      refl_rmax_ds = 0.f;
      refl_rmin_ds = 0.f;
      refl_ep_ds = 0.f;
      refl_nhits_ds = 0;
      refl_q_ds = 0;
      refl_stboundary_ds = 0;
      refl_stinters_ds = 0;
      refl_triggered_ds = 0;
      refl_dp = 0.f;
      refl_dt0 = 0.f;
      refl_dtfront = 0.f;

      crystal_e.clear();
      crystal_t.clear();
    }
  };
} // namespace ConvAna
#endif
