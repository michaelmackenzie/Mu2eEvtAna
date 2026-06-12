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
      return Form("Mu2eEvtAna/file_lists/%s.files", config.full_name_.Data());
    }
  }
  return "";
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
  TString cmd = "ps -eo user,args | grep -E '(make_histograms\\.C|functions\\.C)' | grep -v grep | wc -l";
  TString res = gSystem->GetFromPipe(cmd);
  int count = res.Atoi() > 0 ? res.Atoi() : 0;
  return count > 0 ? count - 1 : 0;
}

void WaitJobs() {
  bool first = true;
  while(CountAnalyzerProcesses() > 0) {
    if(first) { cout << "Waiting for analyzer processes to complete..." << endl; first = false; }
    sleep(10);
  }
}

/**
 * Merge output files using hadd
 */
int MergeOutputFiles(TString output_file, TString merge_list_file) {
  TString cmd = Form("bash ${MUSE_WORK_DIR}/Mu2eEvtAna/scripts/merge_hist_files.sh %s %s -f",
                     output_file.Data(), merge_list_file.Data());
  return gSystem->Exec(cmd.Data());
}

#endif
