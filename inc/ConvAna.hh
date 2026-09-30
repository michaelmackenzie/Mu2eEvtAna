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
#include "Mu2eEvtAna/inc/CutID.hh"
#include "Mu2eEvtAna/inc/CutFlow.hh"

using namespace mu2e;
namespace Mu2eEvtAna {
  class ConvAna : public Mu2eEvtAna {
  public:
    // Sample type, inferred from the analysis name (name_)
    enum Dataset_t {
      kUnknown = 0,
      kCeMinus,  // cele: CeM leading log
      kCePlus,   // cpos: CeP leading log
      kCosmic,   // cry4a: cosmic signal
      kDIO,      // dio00: DIO tail
      kRMCE0,    // rmce0: RMC external conversion, 0 neutron knockout
      kRMCE1,    // rmce1: RMC external conversion, 1 neutron knockout
      kRMCI0,    // rmci0: RMC internal conversion, 0 neutron knockout
      kRMCI1,    // rmci1: RMC internal conversion, 1 neutron knockout
      kRPCE,     // rpce: RPC external conversion
      kRPCI,     // rpci: RPC internal conversion
      kPbar,     // pbar: antiproton resampling
      kFlatPlus, // fpos: flat e+
      kEnsemble  // mds: mock data ensemble
    };

    ConvAna(int verbose = 0);
    ~ConvAna() {};

    void InitHistSelections() override;
    void BookSystematicHist(SysHist_t* Hist, const char* Folder);
    void BookHistograms(TDirectory* dir) override;
    void FillSystematicHist(SysHist_t* Hist);
    bool ProcessEvent() override;
    void InitializeEvent() override;
    void SetEventWeight(); // apply per-event sample weights (antiproton reweighting)
    void InitTrack(const rooutil::Track* track, Track_t& trk_par) override;

    int InitializeInput() override;
    void InitDataset(); // set dataset_ from name_
    static const char* DatasetName(Dataset_t dataset);
    bool IsRMC() const { return dataset_ == kRMCE0 || dataset_ == kRMCE1 || dataset_ == kRMCI0 || dataset_ == kRMCI1; }
    bool IsRMCInternal() const { return dataset_ == kRMCI0 || dataset_ == kRMCI1; }
    bool IsRMCExternal() const { return dataset_ == kRMCE0 || dataset_ == kRMCE1; }
    int InitializeOutput() override;

    void EndJob() override;

    void FillAllHistograms(const int index);
    void InitTreeData();
    CutID Run1ATrackID(Track_t* track);
    bool Run1ACutFlow();
    bool StandardCutFlow();

    bool ValidateVariable(float var, const char* name) {
      if(!std::isfinite(var)) {
        printf(">>> Event %5i/%5i/%6i: Variable %s is non-finite = %f\n", evt_.run_, evt_.subrun_, evt_.event_, name, var);
        return false;
      }
      return true;
    }

    double phase_space_cdf(double k, double kmax, double power) {
      if(kmax <= 0.) return 0.;
      if(power < 0) return 0.;
      k = std::max(0., std::min(kmax, k));
      const double x_1 = k / kmax;
      const double x_2 = 1.;
      const double val_1 = (x_1 - 1.)*std::pow(1-x_1, power)*(power*x_1+x_1+1.);
      const double val_2 = (x_2 - 1.)*std::pow(1-x_2, power)*(power*x_2+x_2+1.);
      const double integral = val_2 - val_1;
      return integral;
    }

    TString OutputFileName() override { return "ConvAna." + name_ + ".root"; }

    CutFlow            cut_flow_                        ; // standard selection
    CutFlow            run1a_cut_flow_                  ; // Run 1A paper selection
    CutFlow            dev_cut_flow_                    ; // For cut-set testing

    Dataset_t          dataset_ = kUnknown              ; // input sample type
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
