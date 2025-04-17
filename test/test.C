// Test the Mu2eEvtAna tools

int test(Long64_t max_entries = 1000) {

  // Retrieve an EventNtuple tree
  auto f_in = TFile::Open("/pnfs/mu2e/tape/phy-nts/nts/mu2e/CeEndpointOnSpillTriggered/MDC2020aq_best_v1_3_v06_03_00/root/ed/df/nts.mu2e.CeEndpointOnSpillTriggered.MDC2020aq_best_v1_3_v06_03_00.001210_00000000.root", "READ");
  if(!f_in) return -1;
  auto t_in = (TTree*) f_in->Get("EventNtuple/ntuple");
  if(!t_in) {
    cout << "Input ntuple not found!\n";
    f_in->ls();
    return 1;
  }

  // Setup a Mu2eEvtAna processing
  Mu2eEvtAna::Mu2eEvtAna ana(0);
  ana.SetInput(t_in);
  ana.SetName("test");

  // Process the tree
  const int status = ana.Process(max_entries);
  cout << "Status code = " << status << endl;

  return 0;
}
