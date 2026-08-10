// Dataset information
#ifndef __MU2E_EVT_ANA_DATASETS_C__
#define __MU2E_EVT_ANA_DATASETS_C__

namespace DATA {
  struct Dataset_t {
    Bool_t process_;
    TString name_;
    TString full_name_;
    Int_t    n_events_;
    Int_t    n_gen_events_;

    Dataset_t(Bool_t process, TString name, TString full_name, Int_t n_events, Int_t n_gen_events)
      : process_(process), name_(name), full_name_(full_name), n_events_(n_events), n_gen_events_(n_gen_events) {}
  };

  vector<Dataset_t> datasets() {
    vector<Dataset_t> datasets;
    datasets.emplace_back(false,  "cele1b0s5r0100", "nts.mu2e.CeMLeadingLogOnSpill-reco-ntuple.MDC2025-002.root"     , 4009075,   10000000);
    datasets.emplace_back(false,  "cry4ab0s5r0100", "nts.mu2e.CosmicSignalOnSpill-reco-ntuple.MDC2025-002.root"      , 4120241,    4437500. /*livetime*/);
    datasets.emplace_back(false,  "dio00b0s5r0100", "nts.mu2e.DIOtail95OnSpill-reco-ntuple.MDC2025-002.root"         , 8780533,   25000000);

    datasets.emplace_back(true ,  "cele1b1s5r0100", "nts.mu2e.CeMLeadingLogMix1BB.MDC2025ar_best_v1_1.root"          , 4009075,   10000000);
    datasets.emplace_back(false,  "cpos1b1s5r0100", "nts.mu2e.CePLeadingLogMix1BB.MDC2025ar_best_v1_1.root"          , 3235878,   10000000);
    datasets.emplace_back(true ,  "cry4ab1s5r0100", "nts.mu2e.CosmicSignalMix1BB.MDC2025ar_best_v1_1.root"           , 4120241,    4437500. /*livetime*/);
    datasets.emplace_back(true ,  "dio00b1s5r0100", "nts.mu2e.DIOtail95Mix1BB.MDC2025ar_best_v1_1.root"              , 9348862,   25000000);
    datasets.emplace_back(true ,  "rmce0b1s5r0100", "nts.mu2e.RMCPhaseSpace0NExternalMix1BB.MDC2025ar_best_v1_1.root", 4952890, 7000000000);
    datasets.emplace_back(false,  "rmce1b1s5r0100", "nts.mu2e.RMCPhaseSpace1NExternalMix1BB.MDC2025ar_best_v1_1.root", 2965072, 7000000000);
    datasets.emplace_back(true ,  "rmci0b1s5r0100", "nts.mu2e.RMCPhaseSpace0NInternalMix1BB.MDC2025ar_best_v1_1.root", 1211753,   50000000);
    datasets.emplace_back(false,  "rmci1b1s5r0100", "nts.mu2e.RMCPhaseSpace1NInternalMix1BB.MDC2025ar_best_v1_1.root",  507641,   50000000);
    datasets.emplace_back(true ,  "rpce1b1s5r0100", "nts.mu2e.RPCExternalPhysicalMix1BB.MDC2025au_best_v1_1.root"    ,  458818, 5000000000);
    datasets.emplace_back(true ,  "rpci1b1s5r0100", "nts.mu2e.RPCInternalPhysicalMix1BB.MDC2025au_best_v1_1.root"    , 1899806,  125000000);

    datasets.emplace_back(true ,  "mds3cb1s5r0100", "nts.mu2e.ensembleMDS3cMix1BB.MDC2025ar_best_v1_1.root"          , 4863968, 1. /*data-ish*/);
    return datasets;
  }
}
#endif // __MU2E_EVT_ANA_DATASETS_C__
