#ifndef __MU2E_EVT_ANA_UTILS_C__
#define __MU2E_EVT_ANA_UTILS_C__
#include "Mu2eEvtAna/scripts/datasets.C"
#include <fstream>
#include <sstream>
#include <algorithm>

/**
 * Get the input file list path for a given dataset
 */
TString GetDatasetFileList(TString dataset_name) {
  for(auto config : DATA::datasets()) {
    if(config.name_ == dataset_name) {
      TString dir = ".";
      if(gSystem) {
        TString muse_dir = gSystem->Getenv("MUSE_WORK_DIR");
        if(muse_dir != "") dir = muse_dir;
      }
      return Form("%s/Mu2eEvtAna/file_lists/%s.files", dir.Data(), config.full_name_.Data());
    }
  }
  return "";
}

/**
 * Check whether the given name corresponds to a known dataset
 */
bool IsKnownDataset(TString name) {
  for(auto config : DATA::datasets()) {
    if(config.name_ == name) return true;
  }
  return false;
}

/**
 * Resolve an input specification into something Mu2eEvtAna::AddFile can read.
 * The input can be:
 *   - a known dataset name  --> the corresponding file list
 *   - a single ntuple file  (*.root, local path or xrootd URL)
 *   - a file list of ntuples (any other existing file)
 * Returns "" if the input cannot be resolved.
 */
TString ResolveInput(TString input) {
  if(input == "") return "";
  TString file_list = GetDatasetFileList(input);
  if(file_list != "") return file_list;
  // not a dataset name, so treat it as a file path
  if(input.BeginsWith("root://") || input.BeginsWith("http://") || input.BeginsWith("https://")) return input;
  if(!gSystem->AccessPathName(input)) return input; // file exists locally
  cout << "Input " << input << " is neither a known dataset nor an existing file!" << endl;
  return "";
}

/**
 * Default output name tag for a given input specification:
 *   - the dataset name for a known dataset
 *   - the file name without its directory or extension otherwise
 */
TString DefaultNameTag(TString input) {
  if(IsKnownDataset(input)) return input;
  TString tag = gSystem->BaseName(input.Data());
  if(tag.EndsWith(".root" )) tag.Remove(tag.Length()-5);
  if(tag.EndsWith(".files")) tag.Remove(tag.Length()-6);
  tag.ReplaceAll("/", "_"); // just in case
  tag.ReplaceAll(" ", "_");
  return tag;
}

/**
 * Split a file list into N parts, writing each part to a separate file.
 * Returns a vector of output file paths.
 */
vector<TString> SplitFileList(TString file_list, int n_parts, TString output_dir, TString prefix) {
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", output_dir.Data(), output_dir.Data()));

  ifstream infile(file_list);
  if(!infile.is_open()) {
    cout << "Error: cannot open file list " << file_list << endl;
    return {};
  }

  vector<string> lines;
  string line;
  while(getline(infile, line)) {
    if(!line.empty()) lines.push_back(line);
  }
  infile.close();

  int n_total = lines.size();
  if(n_total == 0) {
    cout << "Warning: file list " << file_list << " is empty" << endl;
    return {};
  }

  if(n_parts < 1) n_parts = 1;
  if(n_parts > n_total) n_parts = n_total;

  int per_part = (n_total + n_parts - 1) / n_parts;

  vector<TString> output_files;
  for(int part = 0; part < n_parts; ++part) {
    int start = part * per_part;
    int end = min(start + per_part, n_total);

    TString output_file = Form("%s/%s_thread_%i.files", output_dir.Data(), prefix.Data(), part);
    ofstream outfile(output_file);
    for(int i = start; i < end; ++i) {
      outfile << lines[i] << endl;
    }
    outfile.close();

    output_files.push_back(output_file);
  }

  return output_files;
}

/**
 * Submit a ROOT batch job for a specific analyzer function with given parameters
 */
void SubmitJob(TString job_cmd, TString log_file, TString dataset_name, int thread_id = -1) {
  gSystem->Exec(Form("[ ! -d log ] && mkdir log"));
  TString full_cmd = Form("%s >| %s 2>&1 &", job_cmd.Data(), log_file.Data());
  cout << "Submitting job" << (thread_id >= 0 ? Form(" (thread %d)", thread_id) : "")
       << ": " << full_cmd << endl;
  gSystem->Exec(full_cmd.Data());
}

/**
 * Wait for all current analyzer processes to complete
 */
int CountAnalyzerProcesses() {
  TString cmd = "ps -eo user,args | grep -E 'functions\\.C' | grep -v grep | wc -l";
  TString res = gSystem->GetFromPipe(cmd);
  int count = res.Atoi() > 0 ? res.Atoi() : 0;
  return count > 0 ? count : 0;
}

void WaitJobs() {
  int procs = CountAnalyzerProcesses();
  printf("\n");
  while(procs > 0) {
    printf("\033[32mWaiting for analyzer processes to complete, %2i remaining\033[0m\r", procs);
    sleep(10);
    procs = CountAnalyzerProcesses();
  }
  printf("\n");
}

/**
 * Merge output files using hadd
 */
int MergeOutputFiles(TString output_file, TString merge_list_file) {
  TString cmd = Form("bash ${MUSE_WORK_DIR}/Mu2eEvtAna/scripts/merge_hist_files.sh %s %s -f -c",
                     output_file.Data(), merge_list_file.Data());
  return gSystem->Exec(cmd.Data());
}

#endif
