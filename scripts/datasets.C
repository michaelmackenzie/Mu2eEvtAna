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

    //---------------------------------------
    // Run 1A datasets
    //---------------------------------------

    // MDC2025au version, old PID
    datasets.emplace_back(false,  "cele1b1s5r0100", "nts.mu2e.CeMLeadingLogMix1BB.MDC2025au_best_v1_1.root"              , 4141125,   10000000);
    datasets.emplace_back(false,  "cpos1b1s5r0100", "nts.mu2e.CePLeadingLogMix1BB.MDC2025au_best_v1_1.root"              , 3235878,   10000000);
    datasets.emplace_back(false,  "cry4ab1s5r0100", "nts.mu2e.CosmicSignalMix1BB.MDC2025au_best_v1_1.root"               , 4155435,    4437713. /*livetime*/);
    datasets.emplace_back(false,  "dio00b1s5r0100", "nts.mu2e.DIOtail95Mix1BB.MDC2025au_best_v1_1.root"                  , 9368976,   25000000);
    datasets.emplace_back(false,  "rmce0b1s5r0100", "nts.mu2e.RMCPhaseSpace0NExternalMix1BB.MDC2025au_best_v1_1.root"    , 4967393, 7000000000);
    datasets.emplace_back(false,  "rmce1b1s5r0100", "nts.mu2e.RMCPhaseSpace1NExternalMix1BB.MDC2025au_best_v1_1.root"    , 2974188, 7000000000);
    datasets.emplace_back(false,  "rmci0b1s5r0100", "nts.mu2e.RMCPhaseSpace0NInternalMix1BB.MDC2025au_best_v1_1.root"    , 1229006,   50000000);
    datasets.emplace_back(false,  "rmci1b1s5r0100", "nts.mu2e.RMCPhaseSpace1NInternalMix1BB.MDC2025au_best_v1_1.root"    ,  517029,   50000000);
    datasets.emplace_back(false,  "rpce1b1s5r0100", "nts.mu2e.RPCExternalPhysicalMix1BB.MDC2025au_best_v1_1.root"        ,  458818, 5000000000);
    datasets.emplace_back(false,  "rpci1b1s5r0100", "nts.mu2e.RPCInternalPhysicalMix1BB.MDC2025au_best_v1_1.root"        , 1899806,  125000000);
    datasets.emplace_back(false,  "pbar1b1s5r0100", "nts.mu2e.PbarResamplingMix1BB.MDC2025au_best_v1_1.root"             , 6461314,   30000000);
    datasets.emplace_back(false,  "fpos0b1s5r0100", "nts.mu2e.FlatePlusMix1BB.MDC2025au_best_v1_1.root"                  , 2285546,    7500000);
    datasets.emplace_back(false,  "mds3cb1s5r0100", "nts.mu2e.ensembleMDS3cMix1BB.MDC2025au_best_v1_1.root"              , 4926941, 1. /*data-ish*/);
    datasets.emplace_back(false,  "mds3db0s5r0100", "nts.mu2e.ensembleMDS3dOnSpill.MDC2025au_best_v1_1-001.root"         , 9837962, 1. /*data-ish*/);

    // MDC2025au-001 version, new PID (same mcs, so same counts)
    datasets.emplace_back(false,  "cele0b0s5r0101", "nts.mu2e.CeEndpointOnSpill.MDC2025ba_best_v1_5.root"                , 4247008,   10000000);
    datasets.emplace_back(false,  "dio00b0s5r0101", "nts.mu2e.DIOtail95OnSpill.MDC2025ba_best_v1_5.root"                 , 2995328,    8000000);
    datasets.emplace_back(ele  ,  "cele1b1s5r0101", "nts.mu2e.CeMLeadingLogMix1BB.MDC2025au_best_v1_1-001.root"          , 4141125,   10000000);
    datasets.emplace_back(pos  ,  "cpos1b1s5r0101", "nts.mu2e.CePLeadingLogMix1BB.MDC2025au_best_v1_1-001.root"          , 3235878,   10000000);
    datasets.emplace_back(true ,  "cry4ab1s5r0101", "nts.mu2e.CosmicSignalMix1BB.MDC2025au_best_v1_1-001.root"           , 4155435,    4437713. /*livetime*/);
    datasets.emplace_back(ele  ,  "dio00b1s5r0101", "nts.mu2e.DIOtail95Mix1BB.MDC2025au_best_v1_1-001.root"              , 9368976,   25000000);
    datasets.emplace_back(true ,  "rmce0b1s5r0101", "nts.mu2e.RMCPhaseSpace0NExternalMix1BB.MDC2025au_best_v1_1-001.root", 4967393, 7000000000);
    datasets.emplace_back(pos  ,  "rmce1b1s5r0101", "nts.mu2e.RMCPhaseSpace1NExternalMix1BB.MDC2025au_best_v1_1-001.root", 2974188, 7000000000);
    datasets.emplace_back(pos  ,  "rmci0b1s5r0101", "nts.mu2e.RMCPhaseSpace0NInternalMix1BB.MDC2025au_best_v1_1-001.root", 1229006,   50000000);
    datasets.emplace_back(pos  ,  "rmci1b1s5r0101", "nts.mu2e.RMCPhaseSpace1NInternalMix1BB.MDC2025au_best_v1_1-001.root",  517029,   50000000);
    datasets.emplace_back(true ,  "rpce1b1s5r0101", "nts.mu2e.RPCExternalPhysicalMix1BB.MDC2025au_best_v1_1-001.root"    ,  458818, 5000000000);
    datasets.emplace_back(true ,  "rpci1b1s5r0101", "nts.mu2e.RPCInternalPhysicalMix1BB.MDC2025au_best_v1_1-001.root"    , 1899806,  125000000);
    datasets.emplace_back(true ,  "pbar1b1s5r0101", "nts.mu2e.PbarResamplingMix1BB.MDC2025au_best_v1_1-001.root"         , 6461314,   30000000);
    datasets.emplace_back(false,  "fpos0b1s5r0101", "nts.mu2e.FlatePlusMix1BB.MDC2025au_best_v1_1-001.root"              , 2285546,    7500000);
    datasets.emplace_back(false,  "mds3cb1s5r0101", "nts.mu2e.ensembleMDS3cMix1BB.MDC2025au_best_v1_1-001.root"          , 4926941, 1. /*data-ish*/);
    datasets.emplace_back(false,  "mds3db0s5r0101", "nts.mu2e.ensembleMDS3dOnSpill.MDC2025au_best_v1_1-001-001.root"     , 9837962, 1. /*data-ish*/);

    // Flat-spectrum e-/e+/mu- samples for PID performance studies (TrkPID v1 in the trkpid branch)
    datasets.emplace_back(false,  "fele0b1s5r0101", "nts.mu2e.FlateMinusMix1BB.MDC2025au_best_v1_1-001.root"            , 2276505,    7500000);
    datasets.emplace_back(false,  "fele0b0s5r0101", "nts.mu2e.FlateMinusOnSpill.MDC2025au_best_v1_5-001.root"           , 2282410,    7500000);
    datasets.emplace_back(false,  "fpos0b0s5r0101", "nts.mu2e.FlatePlusOnSpill.MDC2025au_best_v1_5-001.root"            , 2291243,    7500000);
    datasets.emplace_back(false,  "fmum0b0s5r0101", "nts.mu2e.FlatMuMinusOnSpill.MDC2025au_best_v1_5-001.root"          , 2022034,    7500000);

    // Nominal-target OnSpill (unmixed) references for the alternate target samples. There is no CeEndpointOnSpill at
    // au_best_v1_5, so the signal reference is CeMLeadingLog. N(events) and N(gen) are from the mcs datasets; the
    // DIO -001 ntuple holds 122 of the 125 mcs files
    datasets.emplace_back(false,  "cele1b0s5r0101", "nts.mu2e.CeMLeadingLogOnSpill.MDC2025au_best_v1_5-001.root"       , 4144677,   10000000);
    // datasets.emplace_back(false,  "dio00b0s5r0101", "nts.mu2e.DIOtail95OnSpill.MDC2025au_best_v1_5-001.root"          , 9386180,   25000000);

    // Alternate target samples (OnSpill): 37 foils + 6 up/6 down (r0110), 49 foils with 1-5 mm (r0120) or 9-41 mm
    // (r0130) holes. N(events) and N(gen) are from the mcs datasets (sum of gen.count); the DIO samples have 8e6 generated
    datasets.emplace_back(false,  "cele0b0s5r0110", "nts.mu2e.CeEndpointOnSpill.MDC2025ba_st37plus6up6dn_best_v1_5.root", 4002040,   10000000);
    datasets.emplace_back(false,  "dio00b0s5r0110", "nts.mu2e.DIOtail95OnSpill.MDC2025ba_st37plus6up6dn_best_v1_5.root" , 2764431,    8000000);
    datasets.emplace_back(false,  "cele0b0s5r0120", "nts.mu2e.CeEndpointOnSpill.MDC2025ba_st49hole1to5_best_v1_5.root"  , 3967395,   10000000);
    datasets.emplace_back(false,  "dio00b0s5r0120", "nts.mu2e.DIOtail95OnSpill.MDC2025ba_st49hole1to5_best_v1_5.root"   , 2717026,    8000000);
    datasets.emplace_back(false,  "cele0b0s5r0130", "nts.mu2e.CeEndpointOnSpill.MDC2025ba_st49hole9to41_best_v1_5.root" , 3961096,   10000000);
    datasets.emplace_back(false,  "dio00b0s5r0130", "nts.mu2e.DIOtail95OnSpill.MDC2025ba_st49hole9to41_best_v1_5.root"  , 2700249,    8000000);

    //---------------------------------------
    // Run 1B datasets
    //---------------------------------------

    // EventNtuple Run1B-010/-011 versions (with the primary branch, needed for spectrum weights). The Run1B-010 ntuples are
    // made from the same Run1Baw_best_v1_5 reco datasets as the earlier versions, so N(events) is unchanged.
    datasets.emplace_back(false, "cele0b0s51r0204", "nts.mu2e.CeEndpoint-KL.Run1B-010.root"                        ,   1184477, 1900000000);
    datasets.emplace_back(false, "fgam0b0s51r0204", "nts.mu2e.FlatGamma-KL.Run1B-010.root"                         ,    897555, 1900000000);
    datasets.emplace_back(false, "csms0b0s51r0204", "nts.mu2e.CosmicCRYAll-KL.Run1B-010.root"                      , 115931437,     2.29e4); // livetime assumed equal to Mix1BB
    datasets.emplace_back(false, "fele0b0s51r0204", "nts.mu2e.FlateMinus-KL.Run1B-010.root"                        ,    570932, 1899000000);
    datasets.emplace_back(false, "mnbs0b1s51r0204", "nts.mu2e.NoPrimaryMix1BB-KL.Run1B-012.root"                   , 100000000,  100000000); // unbiased pileup

    datasets.emplace_back(true , "mnbs1b1s51r0204", "nts.mu2e.NoPrimaryMix1BB-KL.Run1B-011.root"                   , 344196254, 5000000000); // biased to high energy clusters
    datasets.emplace_back(true , "cele0b1s51r0204", "nts.mu2e.CeEndpointMix1BB-KL.Run1B-010.root"                  ,   1326786, 1999000000);
    datasets.emplace_back(true , "fgam0b1s51r0204", "nts.mu2e.FlatGammaMix1BB-KL.Run1B-010.root"                   ,   1039674, 1999000000);
    datasets.emplace_back(true , "csms0b1s51r0204", "nts.mu2e.CosmicCRYAllMix1BB-KL.Run1B-010.root"                ,   2351533,     2.29e4); // cosmic N(gen) is livetime
    datasets.emplace_back(true , "dio00b1s51r0204", "dig.mu2e.DIOMix1BB.Run1Ban_best_v1_4-000.art"                 ,        1.,         1.);
    datasets.emplace_back(true , "fele0b1s51r0204", "nts.mu2e.FlateMinusMix1BB-KL.Run1B-010.root"                  ,    704053, 1998000000);
    datasets.emplace_back(true , "pgamcb1s51r0204", "nts.mu2e.PolyFlatGammaCaloMix1BB-KL.Run1B-010.root"           ,   5249814,  100000000);
    datasets.emplace_back(true , "neut0b1s51r0204", "nts.mu2e.MuCapNeutronTailCaloMix1BB-KL.Run1B-010.root"        ,    154086,  125000000);
    datasets.emplace_back(true , "prot0b1s51r0204", "nts.mu2e.MuCapProtonTailCaloMix1BB-KL.Run1B-010.root"         ,    158688,  100000000);
    datasets.emplace_back(true , "rpce0b1s51r0204", "nts.mu2e.RPCExternalMix1BB-KL.Run1B-010.root"                 ,    643155, 4999841344);
    datasets.emplace_back(true , "rpci0b1s51r0204", "nts.mu2e.RPCInternalMix1BB-KL.Run1B-010.root"                 ,    230559, 4999841344);

    return datasets;
  }
}
#endif // __MU2E_EVT_ANA_DATASETS_C__
