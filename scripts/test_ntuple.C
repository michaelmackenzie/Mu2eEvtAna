void test_ntuple(TString file_name) {
  auto ntuple_ = new TChain("EventNtuple/ntuple");
  Long64_t tree_entries = 0;
  // Check if the given filename contains .root at the end
  if (file_name.EndsWith(".root")) { // assume it's a single file FIXME: Allow for wildcards
    ntuple_->Add(file_name.Data());
    tree_entries = ntuple_->GetEntries();
  } else { // assume it's a file list
    std::ifstream filelist(file_name.Data());
    if (filelist.is_open()) {
      // Count the number of input files in the text file
      const int nfiles = std::count(std::istreambuf_iterator<char>(filelist),
                                    std::istreambuf_iterator<char>(), '\n');
      // Clear any error flags and reset to the beginning
      filelist.clear();
      filelist.seekg(0, std::ios::beg);

      // Add each file to the TChain
      std::string line;
      int ifile = 0;
      Long64_t accumulated_entries = 0;
      while (std::getline(filelist, line)) {
        ++ifile;
        TString line_t(line);
        ntuple_->Add(line_t.Data());
        TFile* f = TFile::Open(line_t, "READ");
        if(f) {
          TTree* t = (TTree*) f->Get("EventNtuple/ntuple");
          if(t) accumulated_entries += t->GetEntriesFast();
          else  printf("!!! File %s has no tree!\n", line_t.Data());
          f->Close();
          printf("\r%s: Loaded %i files of %i with %llu entries", __func__, ifile, nfiles, accumulated_entries);
        } else {
          printf("!!! Unable to read file %s!\n", line_t.Data());
        }
      }
      tree_entries = accumulated_entries;
      filelist.close();
    } else {
      printf("%s: Error! Unable to read input file list %s\n", __func__, file_name.Data());
    }
  }
  printf("Read in %llu tree entries\n", tree_entries);
}
