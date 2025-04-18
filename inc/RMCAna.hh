//
// RMCAna: RMC analysis ntupling
// Michael MacKenzie (2025)

#ifndef MU2EEVTANA_RMCANA_HH
#define MU2EEVTANA_RMCANA_HH

// standard includes

// ROOT includes

// Mu2e Offline includes

// Mu2e EventNtuple includes

// local includes
#include "Mu2eEvtAna/inc/Mu2eEvtAna.hh"

using namespace mu2e;
namespace Mu2eEvtAna {
  class RMCAna : public Mu2eEvtAna {
  public:
    RMCAna(int verbose = 0);
    ~RMCAna() {};

    virtual bool ProcessEvent();
    void InitializeEvent();

    int InitializeInput();
    int InitializeOutput();

    TString OutputFileName() { return "rmcana_" + name_ + ".root"; }

  };
}

#endif
