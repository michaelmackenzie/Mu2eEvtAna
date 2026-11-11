// Dataset information
#ifndef __MU2E_EVT_ANA_DATASETS_C__
#define __MU2E_EVT_ANA_DATASETS_C__

namespace DATA {
  struct Dataset_t {
    Bool_t process_;
    TString name_;
    TString full_name_;

    Dataset_t(Bool_t process, TString name, TString full_name)
      : process_(process), name_(name), full_name_(full_name) {}
  };

  vector<Dataset_t> datasets() {
    vector<Dataset_t> datasets;

    const bool ele = false; // datasets for electron analysis
    const bool pos = true ; // datasets for positron analysis

    //---------------------------------------
    // Run 1A datasets
    //---------------------------------------

    // MDC2025au version, old PID
    datasets.emplace_back(false,  "cele1b1s5r0100", "nts.mu2e.CeMLeadingLogMix1BB.MDC2025au_best_v1_1.root"              );
    datasets.emplace_back(false,  "cpos1b1s5r0100", "nts.mu2e.CePLeadingLogMix1BB.MDC2025au_best_v1_1.root"              );
    datasets.emplace_back(false,  "cry4ab1s5r0100", "nts.mu2e.CosmicSignalMix1BB.MDC2025au_best_v1_1.root"               );
    datasets.emplace_back(false,  "dio00b1s5r0100", "nts.mu2e.DIOtail95Mix1BB.MDC2025au_best_v1_1.root"                  );
    datasets.emplace_back(false,  "rmce0b1s5r0100", "nts.mu2e.RMCPhaseSpace0NExternalMix1BB.MDC2025au_best_v1_1.root"    );
    datasets.emplace_back(false,  "rmce1b1s5r0100", "nts.mu2e.RMCPhaseSpace1NExternalMix1BB.MDC2025au_best_v1_1.root"    );
    datasets.emplace_back(false,  "rmci0b1s5r0100", "nts.mu2e.RMCPhaseSpace0NInternalMix1BB.MDC2025au_best_v1_1.root"    );
    datasets.emplace_back(false,  "rmci1b1s5r0100", "nts.mu2e.RMCPhaseSpace1NInternalMix1BB.MDC2025au_best_v1_1.root"    );
    datasets.emplace_back(false,  "rpce1b1s5r0100", "nts.mu2e.RPCExternalPhysicalMix1BB.MDC2025au_best_v1_1.root"        );
    datasets.emplace_back(false,  "rpci1b1s5r0100", "nts.mu2e.RPCInternalPhysicalMix1BB.MDC2025au_best_v1_1.root"        );
    datasets.emplace_back(false,  "pbar1b1s5r0100", "nts.mu2e.PbarResamplingMix1BB.MDC2025au_best_v1_1.root"             );
    datasets.emplace_back(false,  "fpos0b1s5r0100", "nts.mu2e.FlatePlusMix1BB.MDC2025au_best_v1_1.root"                  );
    datasets.emplace_back(false,  "mds3cb1s5r0100", "nts.mu2e.ensembleMDS3cMix1BB.MDC2025au_best_v1_1.root"              );
    datasets.emplace_back(false,  "mds3db0s5r0100", "nts.mu2e.ensembleMDS3dOnSpill.MDC2025au_best_v1_1-001.root"         );

    // MDC2025au-001 version, new PID (same mcs, so same counts)
    datasets.emplace_back(false,  "cele0b0s5r0101", "nts.mu2e.CeEndpointOnSpill.MDC2025ba_best_v1_5.root"                );
    datasets.emplace_back(false,  "dio00b0s5r0101", "nts.mu2e.DIOtail95OnSpill.MDC2025ba_best_v1_5.root"                 );
    datasets.emplace_back(ele  ,  "cele1b1s5r0101", "nts.mu2e.CeMLeadingLogMix1BB.MDC2025au_best_v1_1-001.root"          );
    datasets.emplace_back(pos  ,  "cpos1b1s5r0101", "nts.mu2e.CePLeadingLogMix1BB.MDC2025au_best_v1_1-001.root"          );
    datasets.emplace_back(true ,  "cry4ab1s5r0101", "nts.mu2e.CosmicSignalMix1BB.MDC2025au_best_v1_1-001.root"           );
    datasets.emplace_back(ele  ,  "dio00b1s5r0101", "nts.mu2e.DIOtail95Mix1BB.MDC2025au_best_v1_1-001.root"              );
    datasets.emplace_back(true ,  "rmce0b1s5r0101", "nts.mu2e.RMCPhaseSpace0NExternalMix1BB.MDC2025au_best_v1_1-001.root");
    datasets.emplace_back(pos  ,  "rmce1b1s5r0101", "nts.mu2e.RMCPhaseSpace1NExternalMix1BB.MDC2025au_best_v1_1-001.root");
    datasets.emplace_back(pos  ,  "rmci0b1s5r0101", "nts.mu2e.RMCPhaseSpace0NInternalMix1BB.MDC2025au_best_v1_1-001.root");
    datasets.emplace_back(pos  ,  "rmci1b1s5r0101", "nts.mu2e.RMCPhaseSpace1NInternalMix1BB.MDC2025au_best_v1_1-001.root");
    datasets.emplace_back(true ,  "rpce1b1s5r0101", "nts.mu2e.RPCExternalPhysicalMix1BB.MDC2025au_best_v1_1-001.root"    );
    datasets.emplace_back(true ,  "rpci1b1s5r0101", "nts.mu2e.RPCInternalPhysicalMix1BB.MDC2025au_best_v1_1-001.root"    );
    datasets.emplace_back(true ,  "pbar1b1s5r0101", "nts.mu2e.PbarResamplingMix1BB.MDC2025au_best_v1_1-001.root"         );
    datasets.emplace_back(false,  "fpos0b1s5r0101", "nts.mu2e.FlatePlusMix1BB.MDC2025au_best_v1_1-001.root"              );
    datasets.emplace_back(false,  "mds3cb1s5r0101", "nts.mu2e.ensembleMDS3cMix1BB.MDC2025au_best_v1_1-001.root"          );
    datasets.emplace_back(false,  "mds3db0s5r0101", "nts.mu2e.ensembleMDS3dOnSpill.MDC2025au_best_v1_1-001-001.root"     );

    // Flat-spectrum e-/e+/mu- samples for PID performance studies (TrkPID v1 in the trkpid branch)
    datasets.emplace_back(false,  "fele0b1s5r0101", "nts.mu2e.FlateMinusMix1BB.MDC2025au_best_v1_1-001.root"            );
    datasets.emplace_back(false,  "fele0b0s5r0101", "nts.mu2e.FlateMinusOnSpill.MDC2025au_best_v1_5-001.root"           );
    datasets.emplace_back(false,  "fpos0b0s5r0101", "nts.mu2e.FlatePlusOnSpill.MDC2025au_best_v1_5-001.root"            );
    datasets.emplace_back(false,  "fmum0b0s5r0101", "nts.mu2e.FlatMuMinusOnSpill.MDC2025au_best_v1_5-001.root"          );

    // Nominal-target OnSpill (unmixed) references for the alternate target samples. There is no CeEndpointOnSpill at
    // au_best_v1_5, so the signal reference is CeMLeadingLog. N(events) and N(gen) are from the mcs datasets; the
    // DIO -001 ntuple holds 122 of the 125 mcs files
    datasets.emplace_back(false,  "cele1b0s5r0101", "nts.mu2e.CeMLeadingLogOnSpill.MDC2025au_best_v1_5-001.root");
    // datasets.emplace_back(false,  "dio00b0s5r0101", "nts.mu2e.DIOtail95OnSpill.MDC2025au_best_v1_5-001.root");

    // Alternate target samples (OnSpill): 37 foils + 6 up/6 down (r0110), 49 foils with 1-5 mm (r0120) or 9-41 mm
    // (r0130) holes. N(events) and N(gen) are from the mcs datasets (sum of gen.count); the DIO samples have 8e6 generated
    datasets.emplace_back(false,  "cele0b0s5r0110", "nts.mu2e.CeEndpointOnSpill.MDC2025ba_st37plus6up6dn_best_v1_5.root");
    datasets.emplace_back(false,  "dio00b0s5r0110", "nts.mu2e.DIOtail95OnSpill.MDC2025ba_st37plus6up6dn_best_v1_5.root" );
    datasets.emplace_back(false,  "cele0b0s5r0120", "nts.mu2e.CeEndpointOnSpill.MDC2025ba_st49hole1to5_best_v1_5.root"  );
    datasets.emplace_back(false,  "dio00b0s5r0120", "nts.mu2e.DIOtail95OnSpill.MDC2025ba_st49hole1to5_best_v1_5.root"   );
    datasets.emplace_back(false,  "cele0b0s5r0130", "nts.mu2e.CeEndpointOnSpill.MDC2025ba_st49hole9to41_best_v1_5.root" );
    datasets.emplace_back(false,  "dio00b0s5r0130", "nts.mu2e.DIOtail95OnSpill.MDC2025ba_st49hole9to41_best_v1_5.root"  );

    //---------------------------------------
    // Run 1B datasets
    //---------------------------------------

    // EventNtuple Run1B-013/-014/-015 versions (fixed RPC pion survival weight, primary branch). The ntuples are made from the
    // same Run1Baw_best_v1_5 reco datasets as the earlier versions; N(events) and N(gen) are summed over the mcs parents of
    // the ntuple files (scripts/make_file_lists.py). As of 2026-11-03, NoPrimary -014 (biased) covers 11000 of the 20000 mcs
    // files and RPCExternal -013 covers 2452 of 5000.
    datasets.emplace_back(false, "cele0b0s51r0204", "nts.mu2e.CeEndpoint-KL.Run1B-013.root"                 );
    datasets.emplace_back(false, "fgam0b0s51r0204", "nts.mu2e.FlatGamma-KL.Run1B-013.root"                  );
    datasets.emplace_back(false, "csms0b0s51r0204", "nts.mu2e.CosmicCRYAll-KL.Run1B-013.root"               );
    datasets.emplace_back(false, "fele0b0s51r0204", "nts.mu2e.FlateMinus-KL.Run1B-013.root"                 );
    datasets.emplace_back(false, "mnbs0b1s51r0204", "nts.mu2e.NoPrimaryMix1BB-KL.Run1B-015.root"            ); // unbiased pileup

    datasets.emplace_back(true , "mnbs1b1s51r0204", "nts.mu2e.NoPrimaryMix1BB-KL.Run1B-014.root"            ); // biased to high energy clusters
    datasets.emplace_back(true , "cele0b1s51r0204", "nts.mu2e.CeEndpointMix1BB-KL.Run1B-013.root"           );
    datasets.emplace_back(true , "fgam0b1s51r0204", "nts.mu2e.FlatGammaMix1BB-KL.Run1B-013.root"            );
    datasets.emplace_back(true , "csms0b1s51r0204", "nts.mu2e.CosmicCRYAllMix1BB-KL.Run1B-013.root"         );
    datasets.emplace_back(true , "dio00b1s51r0204", "dig.mu2e.DIOMix1BB.Run1Ban_best_v1_4-000.art"          );
    datasets.emplace_back(true , "fele0b1s51r0204", "nts.mu2e.FlateMinusMix1BB-KL.Run1B-013.root"           );
    datasets.emplace_back(true , "pgamcb1s51r0204", "nts.mu2e.PolyFlatGammaCaloMix1BB-KL.Run1B-013.root"    );
    datasets.emplace_back(true , "neut0b1s51r0204", "nts.mu2e.MuCapNeutronTailCaloMix1BB-KL.Run1B-013.root" );
    datasets.emplace_back(true , "prot0b1s51r0204", "nts.mu2e.MuCapProtonTailCaloMix1BB-KL.Run1B-013.root"  );
    datasets.emplace_back(true , "rpce0b1s51r0204", "nts.mu2e.RPCExternalMix1BB-KL.Run1B-013.root"          );
    datasets.emplace_back(false, "rpci0b1s51r0204", "nts.mu2e.RPCInternalMix1BB-KL.Run1B-013.root"          );

    return datasets;
  }
}
#endif // __MU2E_EVT_ANA_DATASETS_C__
