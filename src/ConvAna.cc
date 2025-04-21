#include "Mu2eEvtAna/inc/ConvAna.hh"

using namespace mu2e;
namespace Mu2eEvtAna {

  //------------------------------------------------------------------------------------
  // Constructor
  ConvAna::ConvAna(int verbose) : Mu2eEvtAna(verbose) {
  }


  //------------------------------------------------------------------------------------
  // Define the histogram selections
  void ConvAna::InitHistSelections() {
    //Default histogram selections
    evt_hists_[0] = new EventHist_t;
    trk_hists_[0] = new TrackHist_t;
    crv_hists_[0] = new CRVHist_t  ;

    // Basic selection
    evt_hists_[1] = new EventHist_t;
    trk_hists_[1] = new TrackHist_t;
    crv_hists_[1] = new CRVHist_t  ;
  }

  //------------------------------------------------------------------------------------
  // Initialize the input ntuple information
  int ConvAna::InitializeInput() {
    Mu2eEvtAna::InitializeInput();
    return 0;
  }

  //------------------------------------------------------------------------------------
  // Initialize the output ntuple information
  int ConvAna::InitializeOutput() {
    Mu2eEvtAna::InitializeOutput();
    return 0;
  }


  //------------------------------------------------------------------------------------
  // Initialize event information
  void ConvAna::InitializeEvent() {
    Mu2eEvtAna::InitializeEvent();
  }

  //------------------------------------------------------------------------------------
  // Main event-by-event processing
  bool ConvAna::ProcessEvent() {
    FillEventHist(evt_hists_[0]); //all events with well defined inputs
    if(evt_.nde_tracks_ != 1) return false; //exactly one positron or electron
    FillEventHist(evt_hists_[1]);
    return true;
  }
}
