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
#include "Mu2eEvtAna/inc/SysHist_t.hh"
#include "Mu2eEvtAna/inc/Systematics.hh"
#include "Mu2eEvtAna/inc/Tree_t.hh"
#include "Mu2eEvtAna/inc/MVATools.hh"

using namespace mu2e;
namespace Mu2eEvtAna {
  class ConvAna : public Mu2eEvtAna {
  public:
    ConvAna(int verbose = 0);
    ~ConvAna() {};

    void InitHistSelections();
    void BookSystematicHist(SysHist_t* Hist, const char* Folder);
    void BookHistograms(TDirectory* dir);
    bool ProcessEvent();
    void InitializeEvent();
    void InitTrack(rooutil::Track* track, Track_t& trk_par);

    int InitializeInput();
    int InitializeOutput();

    void InitTreeData();
    int  Run1ATrackID(Track_t* track);

    void ValidateVariable(float var, const char* name) {
      if(!std::isfinite(var)) {
        printf(">>> Event %5i/%5i/%6i: Variable %s is non-finite = %f\n", evt_.run_, evt_.subrun_, evt_.event_, name, var);
      }
    }

    TString OutputFileName() { return "ConvAna." + name_ + ".root"; }

    Bool_t             fill_verbose_sys_ = false        ; // add additional info with each systematic

    SysHist_t*         sys_hists_[kMaxHists]            ; // systematic histograms
    TDirectory*        sys_dirs_ [kMaxHists]            ;
    Systematics        systematics_                     ; // systematic information

    Tree_t             tree_                            ; // selected data
    Track_t*           track_ = nullptr                 ;
    CaloCluster_t*     calo_cluster_ = nullptr          ;
    CRVCluster_t*      crv_cluster_ = nullptr           ;

    // MVA info
    int evaluate_mvas_ = 1;
    TMVA::Reader* trkqual_ = nullptr;
    int trkqual_version_ = 0;
    TMVA::Reader* pid_ = nullptr;
    int pid_version_ = 0;
    TMVA::Reader* trkpid_ = nullptr;
    int trkpid_version_ = 0;
    TMVA::Reader* cosmic_id_ = nullptr;
    int cosmic_id_version_ = 0;
  };
}

#endif
