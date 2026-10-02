#ifndef __MU2E_EVT_ANA_FUNCTIONS_C__
#define __MU2E_EVT_ANA_FUNCTIONS_C__
#include "Mu2eEvtAna/scripts/datasets.C"
#include "Mu2eEvtAna/scripts/utils.C"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>

// Global debug level
int debug_level_ = 0;

// Global flags
bool use_xrootd_ = true;

//------------------------------------------------------------------------------------
// Analyzer registry: analyzer function name (e.g. "cnv_ana") --> analyzer class and the libraries it needs.
// The built-in Mu2eEvtAna analyzers are registered by default; analyzers in other packages register with
//   RegisterAnalyzer("bnv_ana", "Mu2eEvtAna::BNVAna", "$MUSE_BUILD_DIR/BNVAna/lib/libbnvana.so");
// The class must derive from Mu2eEvtAna::Mu2eEvtAna and have a constructor taking the verbosity.
struct AnalyzerInfo_t {
  TString class_name_; // analyzer class, including the namespace
  TString libraries_ ; // ':'-separated libraries to load before creating the analyzer (may use environment variables)
};

std::map<TString, AnalyzerInfo_t>& AnalyzerRegistry() {
  static std::map<TString, AnalyzerInfo_t> registry = {
    {"mu2e_ana" , {"Mu2eEvtAna::Mu2eEvtAna", ""}},
    {"rmc_ana"  , {"Mu2eEvtAna::RMCAna"    , ""}},
    {"cnv_ana"  , {"Mu2eEvtAna::ConvAna"   , ""}},
    {"run1b_ana", {"Mu2eEvtAna::Run1BAna", ""}},
  };
  return registry;
}

void RegisterAnalyzer(TString ana_func, TString class_name, TString libraries = "") {
  AnalyzerRegistry()[ana_func] = AnalyzerInfo_t{class_name, libraries};
}

bool IsRegisteredAnalyzer(TString ana_func) {
  return AnalyzerRegistry().count(ana_func) > 0;
}

// Load the libraries an analyzer needs
bool LoadAnalyzerLibraries(TString libraries) {
  if(libraries == "") return true;
  TObjArray* libs = libraries.Tokenize(":");
  bool status = true;
  for(int i = 0; i < libs->GetEntries(); ++i) {
    TString lib = ((TObjString*) libs->At(i))->GetString();
    gSystem->ExpandPathName(lib);
    if(gSystem->Load(lib) < 0) {
      cout << "ERROR: Failed to load library " << lib << endl;
      status = false;
    }
  }
  delete libs;
  return status;
}

// Create an analyzer by class name (only one analyzer instance is kept per process)
Mu2eEvtAna::Mu2eEvtAna* gAnalyzer = nullptr;
Mu2eEvtAna::Mu2eEvtAna* CreateAnalyzer(TString class_name, TString libraries) {
  if(!LoadAnalyzerLibraries(libraries)) return nullptr;
  if(gAnalyzer) {
    delete gAnalyzer;
    gAnalyzer = nullptr;
  }
  gAnalyzer = (Mu2eEvtAna::Mu2eEvtAna*) gROOT->ProcessLine(Form("new %s(0);", class_name.Data()));
  if(!gAnalyzer) cout << "ERROR: Failed to create analyzer " << class_name << endl;
  return gAnalyzer;
}

// Output file prefix of an analyzer (e.g. "ConvAna" for ConvAna.<name>.root)
TString AnalyzerOutputPrefix(Mu2eEvtAna::Mu2eEvtAna* ana) {
  const TString probe = "__name__";
  const TString name = ana->name_;
  ana->SetName(probe);
  TString prefix = ana->OutputFileName();
  ana->SetName(name);
  const int index = prefix.Index("." + probe);
  return (index > 0) ? TString(prefix(0, index)) : TString("EvtAna");
}

// Configure and run an analyzer on an input file list
int RunAnalyzer(Mu2eEvtAna::Mu2eEvtAna* ana, TString input, TString ana_name, Long64_t max_entries, Long64_t first_entry) {
  ana->AddFile(input, max_entries, first_entry);
  ana->SetName(ana_name);
  ana->cache_size_ = 200000000U;
  ana->load_baskets_ = false;
  ana->report_rate_ = 5000;
  ana->use_xrootd_ = use_xrootd_;
  ana->verbose_ = debug_level_;
  return ana->Process(max_entries);
}

// Split a file list into parts
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

// Process one thread's share of the input, using the file list written by ProcessWithThreads.
// Called in a child process via functions.C(ana_func, class_name, libraries, name_tag, mode, max_entries, first_entry, thread_id, xrootd)
int functions(TString ana_func, TString class_name, TString libraries, TString name_tag, int Mode,
              Long64_t max_entries, Long64_t first_entry, int thread_id, bool xrootd = true) {
  use_xrootd_ = xrootd;
  TString input_file = Form("temp/%s_thread_%i.files", name_tag.Data(), thread_id);

  // File list was already split by ProcessWithThreads, just verify it exists
  if(gSystem->AccessPathName(input_file)) {
    cout << "ERROR: Thread file not found: " << input_file << endl;
    return -1;
  }

  auto ana = CreateAnalyzer(class_name, libraries);
  if(!ana) return -1;
  const int status = RunAnalyzer(ana, input_file, Form("%s.%s.m%i.thread_%i", ana_func.Data(), name_tag.Data(), Mode, thread_id),
                                 max_entries, first_entry);
  cout << "Thread " << thread_id << " status = " << status << endl;
  return status;
}

// Generic multi-threaded processing function
//   input   : a known dataset name, a single ntuple file (*.root), or a file list of ntuples
//   name_tag: tag used to name the output files, defaults to the dataset name or the input file name
int ProcessWithThreads(TString ana_func, TString input, int Mode,
                       Long64_t max_entries, Long64_t first_entry, int n_threads,
                       TString name_tag = "") {
  if(!IsRegisteredAnalyzer(ana_func)) {
    cout << "Analyzer " << ana_func << " is not registered (see RegisterAnalyzer)!" << endl;
    return -1;
  }
  const AnalyzerInfo_t info = AnalyzerRegistry()[ana_func];

  TString file_list = ResolveInput(input);
  if(file_list == "") {
    cout << "Input " << input << " not found!" << endl;
    return -1;
  }

  TString dataset = (name_tag != "") ? name_tag : DefaultNameTag(input);
  const bool single_file = file_list.EndsWith(".root");
  if(file_list != input) cout << "Processing dataset " << input << " (file list " << file_list << ")" << endl;
  else                   cout << "Processing input file " << file_list << endl;
  cout << "Output name tag: " << dataset << endl;

  const TString analyzer_name = ana_func;

  if(n_threads <= 1) {
    auto ana = CreateAnalyzer(info.class_name_, info.libraries_);
    if(!ana) return -1;
    const int status = RunAnalyzer(ana, file_list, Form("%s.%s.m%i", analyzer_name.Data(), dataset.Data(), Mode), max_entries, first_entry);
    cout << "Status code = " << status << endl;
    return status;
  }

  // Count number of files in the input and split before submitting threads
  int n_input_files = 0;
  vector<string> all_files;
  if(single_file) { // a single ntuple file, nothing to split
    all_files.push_back(file_list.Data());
    n_input_files = 1;
  } else {
    ifstream infile(file_list);
    string line;
    while(getline(infile, line)) {
      if(!line.empty()) {
        all_files.push_back(line);
        n_input_files++;
      }
    }
    infile.close();
  }

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
    // be very sure we're not using an old file
    gSystem->Exec(Form("[ -f temp/%s ] && rm %s", thread_file.Data(), thread_file.Data()));

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

  // Output file prefix from the analyzer (e.g. ConvAna.<name>.root)
  auto prefix_ana = CreateAnalyzer(info.class_name_, info.libraries_);
  if(!prefix_ana) return -1;
  const TString header = AnalyzerOutputPrefix(prefix_ana);
  delete gAnalyzer;
  gAnalyzer = nullptr;
  TString merged_output = Form("%s.%s.%s.m%i.root", header.Data(), analyzer_name.Data(), dataset.Data(), Mode);


  // XRootD environment variables to prevent permanent hangs
  gSystem->Setenv("XRD_CONNECTIONRETRY", "32");
  gSystem->Setenv("XRD_REQUESTTIMEOUT", "3600");
  gSystem->Setenv("XRD_REDIRECTLIMIT", "255");
  gSystem->Setenv("XRD_STREAMTIMEOUT", "1800");

  std::vector<int> pids;

  for(int t = 0; t < actual_n_threads; ++t) {
    TString cmd = Form("(root.exe -q -b \"${MUSE_WORK_DIR}/Mu2eEvtAna/scripts/functions.C(\\\"%s\\\", \\\"%s\\\", \\\"%s\\\", \\\"%s\\\", %i, %lli, %lli, %i, %i)\" >| log/out_%s_thread_%i.log 2>&1) & echo $!",
                       ana_func.Data(), info.class_name_.Data(), info.libraries_.Data(), dataset.Data(), Mode, max_entries, first_entry, t,
                       (int) use_xrootd_, dataset.Data(), t);
    cout << "Submitting thread " << t << ": " << cmd << endl;
    TString pid_str = gSystem->GetFromPipe(cmd.Data());
    int pid = pid_str.Atoi();
    if (pid > 0) pids.push_back(pid);
    // Stagger launches so the dCache server doesn't get flooded all at once
    gSystem->Sleep(1000);
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
// "input" can be a known dataset name, a single ntuple file, or a file list of ntuples,
// and "name_tag" overrides the tag used to name the output files
int mu2e_ana(TString input, int Mode = 0, Long64_t max_entries = 1e6, Long64_t first_entry = 0, int n_threads = 1, TString name_tag = "") {
  return ProcessWithThreads("mu2e_ana", input, Mode, max_entries, first_entry, n_threads, name_tag);
}

int rmc_ana(TString input, int Mode = 0, Long64_t max_entries = -1, Long64_t first_entry = 0, int n_threads = 1, TString name_tag = "") {
  return ProcessWithThreads("rmc_ana", input, Mode, max_entries, first_entry, n_threads, name_tag);
}

int cnv_ana(TString input, int Mode = 0, Long64_t max_entries = -1, Long64_t first_entry = 0, int n_threads = 1, TString name_tag = "") {
  return ProcessWithThreads("cnv_ana", input, Mode, max_entries, first_entry, n_threads, name_tag);
}

int run1b_ana(TString input, int Mode = 0, Long64_t max_entries = -1, Long64_t first_entry = 0, int n_threads = 1, TString name_tag = "") {
  return ProcessWithThreads("run1b_ana", input, Mode, max_entries, first_entry, n_threads, name_tag);
}

#endif
