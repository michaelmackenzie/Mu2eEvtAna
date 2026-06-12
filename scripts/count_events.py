
import sys
from pathlib import Path
import ROOT

def main():
    # Verify that a file path argument was provided
    if len(sys.argv) < 2:
        print("Error: Missing file list argument.", file=sys.stderr)
        print("Usage: python count_events.py <file_list.txt>", file=sys.stderr)
        sys.exit(1)

    # Get the file path from the first positional argument
    list_file_path = sys.argv[1]
    count = 0
    nproc = 0
    nfiles = 0
    ngoodfiles = 0

    try:
        # Open and read the list file line by line
        with open(list_file_path, 'r', encoding='utf-8') as f:
            for line in f:
                file_name = line.strip()
                if not file_name: continue
                nfiles += 1
                rf = ROOT.TFile.Open(file_name, 'READ')
                if not rf:
                    print(f'>>> Could not open {file_name}')
                    continue
                t = rf.Get('EventNtuple/ntuple')
                if not t:
                    print(f'>>> Could not retreive an ntuple from {file_name}')
                    continue
                h = rf.Get('EventNtuple/n_proc_events')
                if not t:
                    print(f'>>> Could not retreive histogram from {file_name}')
                    continue
                count += t.GetEntries()
                nproc += int(h.Integral()+0.5)
                ngoodfiles += 1
                rf.Close()
            print(f'From {nfiles} input files, found {ngoodfiles} ntuples with {count} events and {nproc} processed events: {count/ngoodfiles:.1f} events per file')
                
    except FileNotFoundError:
        print(f"Error: The file '{list_file_path}' was not found.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
