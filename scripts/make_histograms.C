// Process the input datasets and write out histogram files
#include "Mu2eEvtAna/scripts/datasets.C"
#include "Mu2eEvtAna/scripts/utils.C"
#include "Mu2eEvtAna/scripts/functions.C"

/**
   processes: Number of parallel processes to submit
   dataset  : Specific dataset to process, if empty it will process all enabled datasets.
              Can also be a single ntuple file (*.root) or a file list of ntuples, in which
              case only that input is processed
   mode     : Histogramming mode
   function : histogramming processing function, defined in ana/scripts/
   n_threads: Number of threads per process to split file lists (default 1)
   name_tag : Tag used to name the output files, defaults to the dataset name for a dataset
              or to the input file name for a file input
 **/
int make_histograms(int processes = 1, TString dataset = "", const int mode = 1,
                    const char* function = "cnv_ana", int n_threads = 1, Long64_t max_entries = 1e5,
                    TString name_tag = "") {

  if(processes > 3) {
    cout << "Requested " << processes << " parallel processes, but this exceeds the interactive maximum of about 2-3!\n";
    return 1;
  }
  if(n_threads > 20) {
    cout << "Requested " << n_threads << " parallel threads, but this exceeds the maximum of about 20!\n";
    return 1;
  }


  // Inputs that aren't a known dataset name (a single ntuple file or a file list) are processed directly
  const bool file_input = dataset != "" && !IsKnownDataset(dataset);
  const TString tag = (name_tag != "") ? name_tag : (file_input ? DefaultNameTag(dataset) : dataset);

  // if(tag.BeginsWith("cele")) max_entries = min(max_entries, Long64_t(1e6/n_threads)); // don't need as many signal events
  if(tag.BeginsWith("cpos")) max_entries = min(max_entries, Long64_t(1e6/n_threads));

  if(processes > 1 || n_threads > 1) {
    gSystem->Exec("[ ! -d log ] && mkdir log");
    gSystem->Exec("[ ! -d temp ] && mkdir temp");
    gSystem->Exec("[ ! -d output ] && mkdir output");
  }
  // Check the token is available
  if(use_xrootd_) {
    TString token_str = gSystem->GetFromPipe("timeout 10 getToken; if [ $? -ne 0 ]; then echo Failed; else echo Passed; fi");
    if(!token_str.Contains("Passed")) {
      cout << "Failed to get token!\n";
      return 1;
    }
  }

  if(file_input) { // a single ntuple file or a file list, process it directly
    if(gSystem->AccessPathName(dataset) && !dataset.Contains("://")) {
      cout << "Input " << dataset << " is neither a known dataset nor an existing file!\n";
      return 1;
    }
    printf(" Processing input %s with name tag %s (%i threads)...\n", dataset.Data(), tag.Data(), n_threads);
    if(     strcmp(function, "mu2e_ana") == 0) return mu2e_ana(dataset, mode, max_entries, 0, n_threads, tag);
    else if(strcmp(function, "rmc_ana")  == 0) return rmc_ana (dataset, mode, max_entries, 0, n_threads, tag);
    else if(strcmp(function, "cnv_ana")  == 0) return cnv_ana (dataset, mode, max_entries, 0, n_threads, tag);
    return gInterpreter->ProcessLine(Form("%s(\"%s\", %i, %lld, 0, %i, \"%s\");",
                                          function, dataset.Data(), mode, max_entries, n_threads, tag.Data()));
  }

  auto datasets = DATA::datasets();
  const bool is_run1b = TString(function) == "run1b_ana";

  vector<TString> logs;
  for(auto config : datasets) {
    if(dataset == "") {
      if(!config.process_) continue;
      const bool is_run1b_data = config.name_.EndsWith("r0204");
      if(is_run1b_data != is_run1b) continue;
    }
    if(dataset != "" && config.name_ != dataset) continue;
    // allow overriding the output name tag when a single dataset is requested
    const TString out_tag = (name_tag != "" && dataset != "") ? name_tag : config.name_;
    if(processes > 1) {
      while(CountAnalyzerProcesses() >= processes) sleep(10);
      TString command = Form("root.exe -q -b \"${MUSE_WORK_DIR}/Mu2eEvtAna/scripts/make_histograms.C(0, \\\"%s\\\", %i, \\\"%s\\\", %i)\" >| log/out_%s.log 2>&1 &",
                             config.name_.Data(), mode, function, n_threads, config.name_.Data());
      printf(" Submitting %-20s histogramming (%i threads)...\n", config.name_.Data(), n_threads);
      logs.push_back(Form("log/out_%s.log", config.name_.Data()));
      gSystem->Exec(command.Data());
    } else {
      if(n_threads > 1) {
        if(     strcmp(function, "mu2e_ana") == 0) mu2e_ana(config.name_, mode, max_entries, 0, n_threads, out_tag);
        else if(strcmp(function, "rmc_ana")  == 0) rmc_ana (config.name_, mode, max_entries, 0, n_threads, out_tag);
        else if(strcmp(function, "cnv_ana")  == 0) cnv_ana (config.name_, mode, max_entries, 0, n_threads, out_tag);
        else if(strcmp(function, "run1b_ana") == 0) run1b_ana(config.name_, mode, max_entries, 0, n_threads, out_tag);
      } else {
        gInterpreter->ProcessLine(Form("%s(\"%s\", %i, %lld, 0, 1, \"%s\");",
                                       function, config.name_.Data(), mode, max_entries, out_tag.Data()));
      }
    }
  }

  if(processes > 1) {
    while(CountAnalyzerProcesses() > 0) sleep(10);
  }
  printf("Finished histogramming!\n");

  return 0;
}
