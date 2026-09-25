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

    const bool ele = false; // datasets for electron analysis
    const bool pos = true ; // datasets for positron analysis

    datasets.emplace_back(false,  "cele1b0s5r0100", "nts.mu2e.CeMLeadingLogOnSpill-reco-ntuple.MDC2025-002.root"     , 4009075,   10000000);
    datasets.emplace_back(false,  "cry4ab0s5r0100", "nts.mu2e.CosmicSignalOnSpill-reco-ntuple.MDC2025-002.root"      , 4120241,    4437500. /*livetime*/);
    datasets.emplace_back(false,  "dio00b0s5r0100", "nts.mu2e.DIOtail95OnSpill-reco-ntuple.MDC2025-002.root"         , 8780533,   25000000);

    datasets.emplace_back(ele  ,  "cele1b1s5r0100", "nts.mu2e.CeMLeadingLogMix1BB.MDC2025au_best_v1_1.root"          , 4141125,   10000000);
    datasets.emplace_back(pos  ,  "cpos1b1s5r0100", "nts.mu2e.CePLeadingLogMix1BB.MDC2025au_best_v1_1.root"          , 3235878,   10000000);
    datasets.emplace_back(true ,  "cry4ab1s5r0100", "nts.mu2e.CosmicSignalMix1BB.MDC2025au_best_v1_1.root"           , 4155435,    4437713. /*livetime*/);
    datasets.emplace_back(ele  ,  "dio00b1s5r0100", "nts.mu2e.DIOtail95Mix1BB.MDC2025au_best_v1_1.root"              , 9368976,   25000000);
    datasets.emplace_back(true ,  "rmce0b1s5r0100", "nts.mu2e.RMCPhaseSpace0NExternalMix1BB.MDC2025au_best_v1_1.root", 4967393, 7000000000);
    datasets.emplace_back(pos  ,  "rmce1b1s5r0100", "nts.mu2e.RMCPhaseSpace1NExternalMix1BB.MDC2025au_best_v1_1.root", 2974188, 7000000000);
    datasets.emplace_back(pos  ,  "rmci0b1s5r0100", "nts.mu2e.RMCPhaseSpace0NInternalMix1BB.MDC2025au_best_v1_1.root", 1229006,   50000000);
    datasets.emplace_back(pos  ,  "rmci1b1s5r0100", "nts.mu2e.RMCPhaseSpace1NInternalMix1BB.MDC2025au_best_v1_1.root",  517029,   50000000);
    datasets.emplace_back(true ,  "rpce1b1s5r0100", "nts.mu2e.RPCExternalPhysicalMix1BB.MDC2025au_best_v1_1.root"    ,  458818, 5000000000);
    datasets.emplace_back(true ,  "rpci1b1s5r0100", "nts.mu2e.RPCInternalPhysicalMix1BB.MDC2025au_best_v1_1.root"    , 1899806,  125000000);
    datasets.emplace_back(pos  ,  "pbar1b1s5r0100", "nts.mu2e.PbarResamplingMix1BB.MDC2025au_best_v1_1.root"         , 6461314,   30000000);
    datasets.emplace_back(false,  "fpos0b1s5r0100", "nts.mu2e.FlatePlusMix1BB.MDC2025au_best_v1_1.root"              , 2285546,    7500000);

    datasets.emplace_back(false,  "mds3cb1s5r0100", "nts.mu2e.ensembleMDS3cMix1BB.MDC2025au_best_v1_1.root"          , 4926941, 1. /*data-ish*/);
    datasets.emplace_back(false,  "mds3db0s5r0100", "nts.mu2e.ensembleMDS3dOnSpill.MDC2025au_best_v1_1-001.root"     , 9837962, 1. /*data-ish*/);

    // Run 1B datasets
    // datasets.emplace_back(true , "mnbs0b1s51r0204", "dig.mu2e.NoPrimaryMix1BB.Run1Ban_best_v1_4-000.art"       ,  99995000,   99995000); // unbiased pileup
    datasets.emplace_back(false, "mnbs1b1s51r0204", "nts.mu2e.NoPrimaryMix1BB-KL.Run1B-011.root"                   , 344196254, 5000000000); // biased to high energy clusters
    datasets.emplace_back(false, "cele0b0s51r0204", "nts.mu2e.CeEndpoint-KL.Run1B-010.root"                        ,   1184477, 1900000000);
    datasets.emplace_back(false, "cele0b1s51r0204", "nts.mu2e.CeEndpointMix1BB-KL.Run1Baw_best_v1_5.root"          ,   1326786, 1999000000);
    datasets.emplace_back(false, "fgam0b1s51r0204", "nts.mu2e.FlatGammaMix1BB-KL.Run1Baw_best_v1_5.root"           ,   1039674, 1999000000);
    datasets.emplace_back(false, "csms0b1s51r0204", "nts.mu2e.CosmicCRYAllMix1BB-KL.Run1Baw_best_v1_5.root"        ,   2351533,     2.29e4); // cosmic N(gen) is livetime
    datasets.emplace_back(false, "dio00b1s51r0204", "dig.mu2e.DIOMix1BB.Run1Ban_best_v1_4-000.art"                 ,        1.,         1.);
    datasets.emplace_back(false, "fele0b1s51r0204", "nts.mu2e.FlateMinusMix1BB-KL.Run1Baw_best_v1_5.root"          ,    704053, 1998000000);
    datasets.emplace_back(false, "pgamcb1s51r0204", "nts.mu2e.PolyFlatGammaCaloMix1BB-KL.Run1Baw_best_v1_5.root"   ,   5249814,  100000000);
    datasets.emplace_back(false, "neut0b1s51r0204", "nts.mu2e.MuCapNeutronTailCaloMix1BB-KL.Run1Baw_best_v1_5.root",    154086,  125000000);
    datasets.emplace_back(false, "prot0b1s51r0204", "nts.mu2e.MuCapProtonTailCaloMix1BB-KL.Run1Baw_best_v1_5.root" ,    158688,  100000000);
    datasets.emplace_back(false, "rpce0b1s51r0204", "nts.mu2e.RPCExternalMix1BB-KL.Run1Baw_best_v1_5.root"         ,    643155, 4999841344);

    return datasets;
  }
}
#endif // __MU2E_EVT_ANA_DATASETS_C__
