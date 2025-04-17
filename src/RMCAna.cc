#include "Mu2eEvtAna/inc/RMCAna.hh"

using namespace mu2e;
namespace Mu2eEvtAna {

  //------------------------------------------------------------------------------------
  // Constructor
  RMCAna::RMCAna(int verbose) : Mu2eEvtAna(verbose) {
    //Default to one histogram set/selection
    if(!evt_hists_[0]) evt_hists_[0] = new EventHist_t; //check if it already exists from the defaul constructor
    evt_hists_[1] = new EventHist_t; //basic RMC selection

    // Track collections to keep
    tracks_  .keep_ = false;
    electron_.keep_ =  true;
    muon_    .keep_ = false;
    proton_  .keep_ = false;
  }

  //------------------------------------------------------------------------------------
  // Initialize the input ntuple information
  int RMCAna::InitializeInput() {
    Mu2eEvtAna::InitializeInput();
    return 0;
  }

  //------------------------------------------------------------------------------------

  //------------------------------------------------------------------------------------
  // Initialize the input ntuple information
  int RMCAna::InitializeOutput() {
    Mu2eEvtAna::InitializeOutput();
    return 0;
  }


  //------------------------------------------------------------------------------------
  // Initialize event information
  void RMCAna::InitializeEvent() {
    Mu2eEvtAna::InitializeEvent();
  }

  //------------------------------------------------------------------------------------
  // Main event-by-event processing
  bool RMCAna::ProcessEvent() {
    FillEventHist(evt_hists_[0]); //all events with well defined inputs
    if(evt_.nelectrons_ + evt_.npositrons_ != 1) return false; //exactly one positron or electron
    if(electron_.size_ == 0) {
      throw std::runtime_error(Form("No electron/positron tracks, but requiring a track (size = %i, fields = %i, %i)",
                                    electron_.size_, evt_.nelectrons_, evt_.npositrons_));
    }
    //kinematic selections
    const auto track = electron_.track(0);
    //require a trk+calo hit
    if(track.trk_calo_index_ < 0) return false;

    FillEventHist(evt_hists_[1]);
    return true;
  }
}
