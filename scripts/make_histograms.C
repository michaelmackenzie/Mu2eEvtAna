// Process the input datasets and write out histogram files
#include "Mu2eEvtAna/scripts/datasets.C"
#include "Mu2eEvtAna/scripts/utils.C"
#include "Mu2eEvtAna/scripts/functions.C"

/**
   processes: Number of parallel processes to submit
   dataset  : Specific dataset to process, if empty it will process all enabled datasets
   mode     : Histogramming mode
   function : histogramming processing function, defined in ana/scripts/
   n_threads: Number of threads per process to split file lists (default 1)
 **/
int make_histograms(int processes = 1, TString dataset = "", const int mode = 1,
                    const char* function = "mu2e_ana", int n_threads = 1, Long64_t max_entries = 1e5) {

  if(processes > 3) {
    cout << "Requested " << processes << " parallel processes, but this exceeds the interactive maximum of about 2-3!\n";
    return 1;
  }
  if(processes > 1 || n_threads > 1) {
    gSystem->Exec("[ ! -d log ] && mkdir log");
    gSystem->Exec("[ ! -d temp ] && mkdir temp");
    gSystem->Exec("[ ! -d output ] && mkdir output");
  }

  auto datasets = DATA::datasets();

  vector<TString> logs;
  for(auto config : datasets) {
    if(dataset == "" && !config.process_) continue;
    if(dataset != "" && config.name_ != dataset) continue;
    if(processes > 1) {
      while(CountAnalyzerProcesses() >= processes) sleep(10);
      TString command = Form("root.exe -q -b \"${MUSE_WORK_DIR}/Mu2eEvtAna/scripts/make_histograms.C(0, \\\"%s\\\", %i, \\\"%s\\\", %i)\" >| log/out_%s.log 2>&1 &",
                             config.name_.Data(), mode, function, n_threads, config.name_.Data());
      printf(" Submitting %-20s histogramming (%i threads)...\n", config.name_.Data(), n_threads);
      logs.push_back(Form("log/out_%s.log", config.name_.Data()));
      gSystem->Exec(command.Data());
    } else {
      if(n_threads > 1) {
        if(     strcmp(function, "mu2e_ana") == 0) mu2e_ana(dataset, mode, max_entries, 0, n_threads);
        else if(strcmp(function, "rmc_ana")  == 0) rmc_ana (dataset, mode, max_entries, 0, n_threads);
        else if(strcmp(function, "conv_ana") == 0) conv_ana(dataset, mode, max_entries, 0, n_threads);
      } else {
        gInterpreter->ProcessLine(Form("%s(\"%s\", %i);",
                                      function, config.name_.Data(), mode));
      }
    }
  }

  if(processes > 1) {
    while(CountAnalyzerProcesses() > 0) sleep(10);
  }
  printf("Finished histogramming!\n");

  return 0;
}
