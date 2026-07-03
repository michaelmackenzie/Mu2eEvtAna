#ifndef __MU2E_EVT_ANA_FUNCTIONS_C__
#define __MU2E_EVT_ANA_FUNCTIONS_C__
#include "Mu2eEvtAna/scripts/datasets.C"
#include "Mu2eEvtAna/scripts/utils.C"
#include <fstream>
#include <sstream>
#include <algorithm>

// Global analyzer pointers
Mu2eEvtAna::Mu2eEvtAna* gMu2eAna = nullptr;
Mu2eEvtAna::RMCAna* gRMCAna = nullptr;
Mu2eEvtAna::ConvAna* gConvAna = nullptr;

// Analyzer type identifiers
enum AnalyzerType {
  kMu2eAna = 0,
  kRMCAna = 1,
  kConvAna = 2
};

// Get analyzer name as string
inline TString GetAnalyzerName(AnalyzerType ana_type) {
  switch(ana_type) {
    case kMu2eAna: return "mu2e_ana";
    case kRMCAna: return "rmc_ana";
    case kConvAna: return "cnv_ana";
    default: return "unknown";
  }
}

// Split a file list into parts - called via functions.C(ana_type, dataset, mode, max_entries, first_entry, thread_id)
void SplitFileList(TString file_list, int n_parts, int part, TString output_file) {
  ifstream infile(file_list);
  if(!infile.is_open()) {
    cout << "Error: cannot open file list " << file_list << endl;
    return;
  }

  vector<string> lines;
  string line;
  while(getline(infile, line)) {
    if(!line.empty()) lines.push_back(line);
  }
  infile.close();

  int n_total = lines.size();
  if(n_total == 0) {
    ofstream outfile(output_file);
    outfile.close();
    return;
  }

  if(n_parts < 1) n_parts = n_total;
  if(part < 0 || part >= n_parts) {
    cout << "Error: part must be in range [0, " << n_parts-1 << "]" << endl;
    return;
  }

  int per_part = (n_total + n_parts - 1) / n_parts;
  int start = part * per_part;
  int end = min(start + per_part, n_total);

  ofstream outfile(output_file);
  for(int i = start; i < end; ++i) {
    outfile << lines[i] << endl;
  }
  outfile.close();
}

// Forward declaration for ProcessThreaded
int ProcessThreaded(AnalyzerType ana_type, TString dataset, int Mode, Long64_t max_entries, Long64_t first_entry, int thread_id);

// Thread worker function - called via functions.C(ana_type, dataset, mode, max_entries, first_entry, thread_id)
int functions(int ana_type, TString dataset, int Mode, Long64_t max_entries, Long64_t first_entry, int thread_id);

// Thread worker function implementation
int functions(int ana_type, TString dataset, int Mode, Long64_t max_entries, Long64_t first_entry, int thread_id) {
  return ProcessThreaded((AnalyzerType)ana_type, dataset, Mode, max_entries, first_entry, thread_id);
}

int ProcessThreaded(AnalyzerType ana_type, TString dataset, int Mode, Long64_t max_entries, Long64_t first_entry, int thread_id) {
  TString file_list = GetDatasetFileList(dataset);
  if(file_list == "") {
    cout << "Dataset " << dataset << " not found!" << endl;
    return -1;
  }

  TString input_file = Form("temp/%s_thread_%i.files", dataset.Data(), thread_id);

  // File list was already split by ProcessWithThreads, just verify it exists
  if(gSystem->AccessPathName(input_file)) {
    cout << "ERROR: Thread file not found: " << input_file << endl;
    return -1;
  }

  TString analyzer_name = GetAnalyzerName(ana_type);
  TString ana_name = Form("%s.%s.m%i.thread_%i", analyzer_name.Data(), dataset.Data(), Mode, thread_id);

  Mu2eEvtAna::Mu2eEvtAna* ana = nullptr;
  switch(ana_type) {
    case kMu2eAna:
      if(gMu2eAna) delete gMu2eAna;
      gMu2eAna = new Mu2eEvtAna::Mu2eEvtAna(0);
      ana = gMu2eAna;
      break;
    case kRMCAna:
      if(gRMCAna) delete gRMCAna;
      gRMCAna = new Mu2eEvtAna::RMCAna(0);
      ana = gRMCAna;
      break;
    case kConvAna:
      if(gConvAna) delete gConvAna;
      gConvAna = new Mu2eEvtAna::ConvAna(0);
      ana = gConvAna;
      break;
  }

  ana->AddFile(input_file, max_entries, first_entry);
  ana->SetName(ana_name);
  ana->cache_size_ = 200000000U;
  ana->load_baskets_ = false;
  ana->report_rate_ = 5000;

  const int status = ana->Process(max_entries);
  cout << "Thread " << thread_id << " status = " << status << endl;

  return status;
}

// Generic multi-threaded processing function
int ProcessWithThreads(AnalyzerType ana_type, TString dataset, int Mode,
                       Long64_t max_entries, Long64_t first_entry, int n_threads) {
  TString file_list = GetDatasetFileList(dataset);
  if(file_list == "") {
    cout << "Dataset " << dataset << " not found!" << endl;
    return -1;
  }

  TString analyzer_name = GetAnalyzerName(ana_type);

  if(n_threads <= 1) {
    Mu2eEvtAna::Mu2eEvtAna* ana = nullptr;
    switch(ana_type) {
      case kMu2eAna:
        if(gMu2eAna) delete gMu2eAna;
        gMu2eAna = new Mu2eEvtAna::Mu2eEvtAna(0);
        ana = gMu2eAna;
        break;
      case kRMCAna:
        if(gRMCAna) delete gRMCAna;
        gRMCAna = new Mu2eEvtAna::RMCAna(0);
        ana = gRMCAna;
        break;
      case kConvAna:
        if(gConvAna) delete gConvAna;
        gConvAna = new Mu2eEvtAna::ConvAna(0);
        ana = gConvAna;
        break;
    }

    ana->AddFile(file_list, max_entries, first_entry);
    ana->SetName(Form("%s.%s.m%i", analyzer_name.Data(), dataset.Data(), Mode));
    ana->cache_size_ = 200000000U;
    ana->load_baskets_ = false;
    ana->report_rate_ = 5000;

    const int status = ana->Process(max_entries);
    cout << "Status code = " << status << endl;
    return status;
  }

  // Count number of files in the input list and split before submitting threads
  ifstream infile(file_list);
  int n_input_files = 0;
  string line;
  vector<string> all_files;
  while(getline(infile, line)) {
    if(!line.empty()) {
      all_files.push_back(line);
      n_input_files++;
    }
  }
  infile.close();

  // Adjust n_threads if fewer input files than requested threads
  int actual_n_threads = n_threads;
  if(n_input_files < n_threads) {
    actual_n_threads = n_input_files;
    cout << "Note: Reducing threads from " << n_threads << " to " << actual_n_threads
         << " (only " << n_input_files << " input files)" << endl;
  }

  gSystem->Exec("[ ! -d log ] && mkdir log");
  gSystem->Exec("[ ! -d temp ] && mkdir temp");

  // Split file list before submitting threads
  vector<TString> thread_files;
  for(int t = 0; t < actual_n_threads; ++t) {
    TString thread_file = Form("temp/%s_thread_%i.files", dataset.Data(), t);
    thread_files.push_back(thread_file);

    int per_part = (n_input_files + actual_n_threads - 1) / actual_n_threads;
    int start = t * per_part;
    int end = min(start + per_part, n_input_files);

    ofstream outfile(thread_file);
    for(int i = start; i < end; ++i) {
      outfile << all_files[i] << endl;
    }
    outfile.close();
  }

  cout << "Processing " << dataset << " with " << actual_n_threads << " threads" << endl;

  TString header = "EvtAna";
  switch(ana_type) {
  case kMu2eAna: header = "EvtAna"; break;
  case kRMCAna: header = "RMCAna"; break;
  case kConvAna: header = "ConvAna"; break;
  }
  TString merged_output = Form("%s.%s.%s.m%i.root", header.Data(), analyzer_name.Data(), dataset.Data(), Mode);

  TString threaded_func;
  switch(ana_type) {
    case kMu2eAna: threaded_func = "mu2e_ana_threaded"; break;
    case kRMCAna: threaded_func = "rmc_ana_threaded"; break;
    case kConvAna: threaded_func = "cnv_ana_threaded"; break;
  }

  // XRootD environment variables to prevent permanent hangs
  gSystem->Setenv("XRD_CONNECTIONRETRY", "32");
  gSystem->Setenv("XRD_REQUESTTIMEOUT", "3600");
  gSystem->Setenv("XRD_REDIRECTLIMIT", "255");
  gSystem->Setenv("XRD_STREAMTIMEOUT", "1800");

  std::vector<int> pids;

  for(int t = 0; t < actual_n_threads; ++t) {
    TString cmd = Form("(root.exe -q -b \"${MUSE_WORK_DIR}/Mu2eEvtAna/scripts/functions.C(%i, \\\"%s\\\", %i, %lli, %lli, %i)\" >| log/out_%s_thread_%i.log 2>&1) & echo $!",
                       ana_type, dataset.Data(), Mode, max_entries, first_entry, t, dataset.Data(), t);
    cout << "Submitting thread " << t << ": " << cmd << endl;
    TString pid_str = gSystem->GetFromPipe(cmd.Data());
    int pid = pid_str.Atoi();
    if (pid > 0) pids.push_back(pid);
    // Stagger launches so the dCache server doesn't get flooded all at once
    gSystem->Sleep(2000);
  }

  // Wait for the jobs to finish
  // WaitJobs();
  printf("\n");
  while(!pids.empty()) {
    printf("\033[32mWaiting for analyzer processes to complete, %2lu remaining\033[0m\r", pids.size());
    fflush(stdout);
    gSystem->Sleep(2000); // Check every 2 seconds

    for (auto it = pids.begin(); it != pids.end(); ) {
      TString check_cmd = Form("[ -d /proc/%i ] && echo 1 || echo 0", *it);
      TString active = gSystem->GetFromPipe(check_cmd.Data());
      if (active.Atoi() == 0) {
        it = pids.erase(it); // Remove finished PID from our watch list
      } else {
        ++it;
      }
    }
  }
  printf("\nAll jobs finished successfully.\n");

  TString merge_list = Form("temp/%s_merge_list.txt", dataset.Data());
  ofstream ml(merge_list);
  int nfinished = 0;
  for(int t = 0; t < actual_n_threads; ++t) {
    const int status = gSystem->Exec(Form("grep -q 'Thread %i status = 0' log/out_%s_thread_%i.log", t, dataset.Data(), t));
    if(status == 0) {
      ml << Form("%s.%s.%s.m%i.thread_%i.root\n", header.Data(), analyzer_name.Data(), dataset.Data(), Mode, t);
      ++nfinished;
    } else {
      cout << ">>> Error! " << dataset << "_thread_" << t << " did not finish properly!\n";
    }
  }
  ml.close();

  if(nfinished > 0) {
    MergeOutputFiles(merged_output, merge_list);
    cout << "Multi-thread processing complete. Output: " << merged_output << endl;
  } else {
    cout << "Multi-thread processing failed! Did not produce output: " << merged_output << endl;
  }
  return 0;
}

// Convenience wrapper functions for CINT
int mu2e_ana(TString dataset, int Mode = 0, Long64_t max_entries = 1e6, Long64_t first_entry = 0, int n_threads = 1) {
  return ProcessWithThreads(kMu2eAna, dataset, Mode, max_entries, first_entry, n_threads);
}

int rmc_ana(TString dataset, int Mode = 0, Long64_t max_entries = -1, Long64_t first_entry = 0, int n_threads = 1) {
  return ProcessWithThreads(kRMCAna, dataset, Mode, max_entries, first_entry, n_threads);
}

int cnv_ana(TString dataset, int Mode = 0, Long64_t max_entries = -1, Long64_t first_entry = 0, int n_threads = 1) {
  return ProcessWithThreads(kConvAna, dataset, Mode, max_entries, first_entry, n_threads);
}

// C-compatible wrapper functions for ROOT
extern "C" {
  int mu2e_ana_threaded(TString dataset, int Mode, Long64_t max_entries, Long64_t first_entry, int thread_id) {
    return ProcessThreaded(kMu2eAna, dataset, Mode, max_entries, first_entry, thread_id);
  }

  int rmc_ana_threaded(TString dataset, int Mode, Long64_t max_entries, Long64_t first_entry, int thread_id) {
    return ProcessThreaded(kRMCAna, dataset, Mode, max_entries, first_entry, thread_id);
  }

  int cnv_ana_threaded(TString dataset, int Mode, Long64_t max_entries, Long64_t first_entry, int thread_id) {
    return ProcessThreaded(kConvAna, dataset, Mode, max_entries, first_entry, thread_id);
  }
}

#endif
