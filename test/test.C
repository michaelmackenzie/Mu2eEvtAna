// Test the Mu2eEvtAna tools

int test(Long64_t max_entries = -1) {

  TString file_name = "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CeEndpointOnSpillTriggered/MDC2020aq_best_v1_3_v06_03_00/root/ed/df/nts.mu2e.CeEndpointOnSpillTriggered.MDC2020aq_best_v1_3_v06_03_00.001210_00000000.root";

  // Setup a Mu2eEvtAna processing
  Mu2eEvtAna::Mu2eEvtAna ana(0);
  ana.AddFile(file_name);
  ana.SetName("test");

  // Process the tree
  const int status = ana.Process(max_entries);
  cout << "Status code = " << status << endl;

  return 0;
}
