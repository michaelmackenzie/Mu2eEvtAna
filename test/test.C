// Test the Mu2eEvtAna tools
Mu2eEvtAna::Mu2eEvtAna* gAna;

int test(Long64_t max_entries = 1e4, Long64_t first_entry = 0) {

  TString file_list = "Mu2eEvtAna/file_lists/nts.mu2e.CeMLeadingLogMix1BB.MDC2025au_best_v1_1-001.root.files";

  // Setup a Mu2eEvtAna processing
  gAna = new Mu2eEvtAna::Mu2eEvtAna((max_entries < 10) ? (max_entries == 1) ? 10 : 2 : 0);
  gAna->AddFile(file_list, max_entries, first_entry);
  gAna->SetName("test");
  gAna->cache_size_ = 2000000000U; //2 GB cache
  gAna->load_baskets_ = false; // whether or not to load tree baskets when a new file is opened

  // Process the tree
  const int status = gAna->Process(max_entries);
  cout << "Status code = " << status << endl;

  return 0;
}
