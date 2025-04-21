//
// ConvAna: Conversion search analysis ntupling/histogramming
// Michael MacKenzie (2025)

#ifndef MU2EEVTANA_CONVANA_HH
#define MU2EEVTANA_CONVANA_HH

// standard includes

// ROOT includes

// Mu2e Offline includes

// Mu2e EventNtuple includes

// local includes
#include "Mu2eEvtAna/inc/Mu2eEvtAna.hh"

using namespace mu2e;
namespace Mu2eEvtAna {
  class ConvAna : public Mu2eEvtAna {
  public:
    ConvAna(int verbose = 0);
    ~ConvAna() {};

    void InitHistSelections();
    bool ProcessEvent();
    void InitializeEvent();

    int InitializeInput();
    int InitializeOutput();

    TString OutputFileName() { return "convana_" + name_ + ".root"; }

  };
}

#endif
