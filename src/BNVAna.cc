#include "Mu2eEvtAna/inc/BNVAna.hh"

using namespace mu2e;
namespace Mu2eEvtAna {

  //------------------------------------------------------------------------------------
  // Constructor
  BNVAna::BNVAna(int verbose) : ConvAna(verbose) {
    hist_track_id_ = kBNVTrackID; // show the BNV ID bits in the track ID histograms
  }

  //------------------------------------------------------------------------------------
  // BNV histogram sets: e+- above the DIO endpoint (charge-blind) and in-time e+e- pairs
  void BNVAna::InitHistSelections() {
    struct hist_info_t {
      TString _dsc; // description of the selection
      bool    _trk; // track histograms
      bool    _crv; // CRV histograms
      hist_info_t(TString dsc = "", bool trk = false, bool crv = false) : _dsc(dsc), _trk(trk), _crv(crv) {}
    };

    hist_info_t* hist_sets[kMaxHists];
    for(int i = 0; i < kMaxHists; ++i) hist_sets[i] = nullptr;

    //                                  description                               trk    crv
    hist_sets[  0] = new hist_info_t("All events"                            , false, false);
    hist_sets[  1] = new hist_info_t("Exactly one downstream e+-"            , false, false);
    hist_sets[200] = new hist_info_t("BNV: e+-, p > 100, no ID"              ,  true, false);
    hist_sets[201] = new hist_info_t("BNV: e+-, BNV ID, p > 105"             ,  true,  true);
    hist_sets[202] = new hist_info_t("BNV: e-, BNV ID, p > 105"              ,  true, false);
    hist_sets[203] = new hist_info_t("BNV: e+, BNV ID, p > 105"              ,  true, false);
    hist_sets[204] = new hist_info_t("BNV: e+-, BNV ID w/o cosmic ID"        ,  true, false);
    hist_sets[205] = new hist_info_t("BNV: e+-, BNV ID, p > 140"             ,  true, false);
    hist_sets[206] = new hist_info_t("BNV: e+-, BNV ID + trigger"            ,  true, false);
    hist_sets[207] = new hist_info_t("BNV: in-time e+e- track pair"          ,  true, false);
    hist_sets[208] = new hist_info_t("BNV: e+-, std ID w/o p window, p > 105",  true, false);

    for(int i = 0; i < kMaxHists; ++i) {
      if(!hist_sets[i]) continue;
      evt_hists_[i] = new EventHist_t;
      if(hist_sets[i]->_trk) trk_hists_[i] = new TrackHist_t;
      if(hist_sets[i]->_crv) crv_hists_[i] = new CRVHist_t;
    }
  }

  //------------------------------------------------------------------------------------
  // Track IDs: the standard e-/e+ IDs (indices 0/1) and the BNV ID (kBNVTrackID)
  void BNVAna::SetTrackIDs(Track_t* track) {
    Mu2eEvtAna::SetTrackIDs(track);
    if(!track || !track->track_) return;
    track->SetID(BNVTrackID(track), kBNVTrackID);
  }

  //------------------------------------------------------------------------------------
  // BNV track ID, a charge-blind loosening of the standard IDs for e+- signals above the DIO endpoint
  // (single e- lines up to ~180 MeV/c and dark cascade e+e- pairs, see BNVAna/doc/BNVMuonCapture.md).
  // Changes from the standard IDs, motivated by the signal samples (ConvAna set 200, p > 100 MeV/c):
  //  - no momentum window: the BNV sets apply their own momentum requirements
  //  - R(max) 400-680 mm for both charges, and no OPA/TSDA intersection requirement: these cost ~20% of the
  //    e- line signal, and energy loss in them can only lower the momentum, so it cannot push DIO above 105 MeV/c
  //  - p(chi2) > 1e-8 for both charges (the e- value 1e-5 costs ~10%)
  //  - calorimeter cluster not required: the tracker-only PID is used without one (as in the e+ ID)
  //  - tan(dip) 0.5-2.5: signals at lower dark-state masses are more forward
  //  - kinematic cosmic ID > 0.5 for both charges (the e- value 0.85 costs ~30%; its MVA was trained on
  //    conversion-like kinematics), to be revisited with the cosmic-ray background
  //  - no in-time multi-track veto (it would reject reconstructed signal e+e- pairs); the upstream
  //    reflection veto is kept
  // Charge-specific track quality, PID, and timing requirements are those of the standard IDs.
  CutID BNVAna::BNVTrackID(Track_t* track) {
    if(!track || !track->track_) return 0;
    CutID ID;

    const bool electron = track->Charge() < 0;

    if(track->RMaxFront() < 400. || track->RMaxFront() > 680.)     ID.SetBit(kRMax);
    if(track->FitCon() < 1.e-8)                                    ID.SetBit(kFitCon);
    if(track->TanDipFront() < 0.5 || track->TanDipFront() > 2.5)   ID.SetBit(kTDip);
    if(track->CosmicID() > -100.f && track->CosmicID() < 0.5f)     ID.SetBit(kCosmicID);
    if(!track->STBoundary())                                       ID.SetBit(kD0); // consistent with stopping target
    if(track->TFront() < 475. || track->TFront() > 1650.)          ID.SetBit(kT0Loose); // for control regions

    if(electron) {
      if(track->AltTrkQual() > -10. && track->AltTrkQual() < 0.2)  ID.SetBit(kTrkQual);
      if(track->TFront() < 540. || track->TFront() > 1650.)        ID.SetBit(kT0);
    } else {
      if(track->TrkQual() > -10. && track->TrkQual() < 0.015)      ID.SetBit(kTrkQual);
      if(track->TFront() < 500. || track->TFront() > 1650.)        ID.SetBit(kT0);
    }

    // PID: full PID with a calorimeter cluster, tracker-only PID without one
    if(track->ECluster() <= 0.) {
      if(track->TrkPID() < -100.f)                                 ID.SetBit(kClusterE); // no score --> fail it
      else if(track->TrkPID() < 0.15f)                             ID.SetBit(kPID);
    } else if(track->AltPID() < ((electron) ? 0.5f : 0.10f))       ID.SetBit(kPID);

    // Upstream reflection veto (no in-time multi-track veto)
    for(int i = 0; i < evt_.ntracks_; ++i) {
      const auto alt_trk = &tracks_[i];
      if(alt_trk == track || !alt_trk->IsGood()) continue;
      if(alt_trk->PZFront() < 0.f) {
        const float dt = track->TFront() - alt_trk->TFront();
        if(dt > 40.f && dt < 110.f) {
          ID.SetBit(kUpstream);
          break;
        }
      }
    }

    // CRV veto, as in the standard IDs
    if(track->stub_) {
      auto stub = track->stub_;
      const float deltat_st     = track->TFront() - stub->TimeViaSTBack();
      const float deltat_calo   = track->TFront() - stub->TimeViaCaloFront();
      const float deltat_crv    = track->TFront() - stub->Time();
      const float min_extrap_dt(-50.f), max_extrap_dt(60.f);
      if((deltat_st   > min_extrap_dt && deltat_st   < max_extrap_dt) ||
         (deltat_calo > min_extrap_dt && deltat_calo < max_extrap_dt) ||
         (deltat_crv > -25.f && deltat_crv < 0.f))                 ID.SetBit(kCRV);
    }
    return ID;
  }

  //------------------------------------------------------------------------------------
  // Main event-by-event processing
  bool BNVAna::ProcessEvent() {
    ValidateTracks();
    SetEventWeight();

    // Recompute the BNV ID now that all tracks in the event are initialized (its upstream reflection veto
    // looks at the other tracks, which may not be set yet when SetTrackIDs runs from InitTrack)
    for(int itrk = 0; itrk < evt_.ntracks_; ++itrk) {
      if(tracks_[itrk].IsGood()) tracks_[itrk].SetID(BNVTrackID(&tracks_[itrk]), kBNVTrackID);
    }

    FillEventHist(evt_hists_[0]); // all events with well defined inputs

    bnv_cut_flow_.ResetEvent();
    bnv_cut_flow_.Increment("All");

    // An online helix trigger fired
    const bool triggered = trigger_.FiredAPR() || trigger_.FiredCPR();

    for(int itrk = 0; itrk < evt_.ntracks_; ++itrk) {
      track_ = &tracks_[itrk];
      if(!track_->IsGood()) continue; // if not a properly fit track, skip it
      if(track_->PFront() <= 0.) continue;
      if(std::abs(track_->FitPDG()) != 11 || track_->PZFront() <= 0.f) continue; // downstream e+- fits
      bnv_cut_flow_.Increment("downstream_e");

      const float p_trk = track_->PFront();
      if(p_trk > 100.f) {
        bnv_cut_flow_.Increment("p_100");
        FillAllHistograms(200);
      }

      const auto bnv_id = track_->ID(kBNVTrackID);
      if(p_trk > 105.f) {
        bnv_cut_flow_.Increment("p_105");
        if(bnv_id.Passes()) {
          bnv_cut_flow_.Increment("bnv_id");
          FillAllHistograms(201);
          FillAllHistograms((track_->Charge() < 0) ? 202 : 203);
          if(p_trk > 140.f) FillAllHistograms(205); // above the RPC endpoint
          if(triggered)     FillAllHistograms(206);
        }
        if(bnv_id.ID(~(1 << kCosmicID)) == 0) FillAllHistograms(204); // for cosmic ID studies

        // Charge-appropriate standard ID without its momentum window or multi-track veto, for comparison
        bool upstream_veto(true);
        for(int i = 0; i < evt_.ntracks_; ++i) {
          const auto alt_trk = &tracks_[i];
          if(i == itrk || !alt_trk->IsGood() || alt_trk->PZFront() >= 0.f) continue;
          const float dt = track_->TFront() - alt_trk->TFront();
          upstream_veto &= dt < 40.f || dt > 110.f;
        }
        const auto std_id = track_->ID((track_->Charge() < 0) ? 0 : 1);
        if(upstream_veto && std_id.ID(~((1 << kP) | (1 << kUpstream))) == 0) FillAllHistograms(208);
      }

      // In-time e+e- pair (no ID): fill once per pair, from the positron
      if(track_->Charge() > 0) {
        for(int i = 0; i < evt_.ntracks_; ++i) {
          if(i == itrk) continue;
          const auto alt_trk = &tracks_[i];
          if(!alt_trk->IsGood() || std::abs(alt_trk->FitPDG()) != 11 || alt_trk->PZFront() <= 0.f) continue;
          if(alt_trk->Charge() >= 0 || alt_trk->PFront() <= 0.f) continue;
          if(std::fabs(track_->TFront() - alt_trk->TFront()) < 25.f) {
            bnv_cut_flow_.Increment("ee_pair");
            FillAllHistograms(207);
            break;
          }
        }
      }
    }

    if(evt_.nde_tracks_ != 1) return false; // exactly one downstream e+-
    FillEventHist(evt_hists_[1]);
    return true;
  }

  //------------------------------------------------------------------------------------
  void BNVAna::EndJob() {
    printf("BNVAna::%s\n", __func__);
    bnv_cut_flow_.Print();
  }
}
