// Test the Mu2eEvtAna tools
Mu2eEvtAna::Mu2eEvtAna* gAna;

int test(Long64_t max_entries = 1e4) {

  // TString file_list = "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CeEndpointOnSpillTriggered/MDC2020aq_best_v1_3_v06_03_00/root/ed/df/nts.mu2e.CeEndpointOnSpillTriggered.MDC2020aq_best_v1_3_v06_03_00.001210_00000000.root";
  // TString file_list = "Mu2eEvtAna/file_lists/nts.mu2e.CeEndpointOnSpillTriggered.MDC2020aq_best_v1_3_v06_03_00.root.files";
  TString file_list = "Mu2eEvtAna/file_lists/nts.mu2e.CosmicCRYSignalAllOnSpillTriggered.MDC2020aq_best_v1_3_v06_03_00.root.files";
  // TString file_list = "Mu2eEvtAna/file_lists/nts.mu2e.RPCExternalOnSpillTriggered.MDC2020aq_best_v1_3_v06_03_00.root.files";

  // Setup a Mu2eEvtAna processing
  gAna = new Mu2eEvtAna::Mu2eEvtAna((max_entries < 10) ? (max_entries == 1) ? 10 : 2 : 0);
  gAna->AddFile(file_list, max_entries);
  gAna->SetName("test");

  // Process the tree
  const int status = gAna->Process(max_entries);
  cout << "Status code = " << status << endl;

  return 0;
}
