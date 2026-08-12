#!/usr/bin/env python3
"""
Cut-Flow Comparison Script

This script reads multiple cut-flow files with different syntax styles,
extracts the cut names and event counts row-by-row, aligns them, and flags
any discrepancies in counts to help quickly catch selection or processing errors.

Supported Styles:
Style 1: "Cut  X:  name  count: abs rel"
Style 2: "name  N/A/NaN  count  abs  rel  description"
"""

import sys
import re
import os
from typing import List, Dict, Tuple, Optional, Any

# --- REGEX PATTERNS FOR PARSING ---
# Style 1: Matches "Cut  1:      a_track    3830956:   92.5      92.51"
STYLE1_REGEX = re.compile(r"^\s*Cut\s+\d+:\s*([a-zA-Z0-9_\-\s]+?)\s+(\d+):")

# Style 2: Matches "has_a_track   NaN         3830956     92.510030..."
# Uses the N/A or NaN column indicator as a robust anchor point
STYLE2_REGEX = re.compile(r"^\s*(.*?)\s+(N/A|NaN)\s+(\d+)\s+([\d\.]+)")


def parse_cut_flow_file(filepath: str) -> List[Dict[str, Any]]:
    """
    Reads a cut-flow log file, auto-detects its syntax style, and extracts
    the cuts as a list of dictionaries containing 'name' and 'count'.
    """
    parsed_cuts = []
    
    if not os.path.exists(filepath):
        print(f"Error: File not found at '{filepath}'", file=sys.stderr)
        return []

    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Determine style based on content markers
    detected_style = None
    for line in lines:
        if "Cut " in line and ":" in line:
            detected_style = 1
            break
        if "Total Cut Flow" in line or ("Cut Group" in line and "Events Passing" in line):
            detected_style = 2
            break

    # If header keywords are missing, fallback to attempting regex matching
    if detected_style is None:
        for line in lines:
            if STYLE1_REGEX.search(line):
                detected_style = 1
                break
            if STYLE2_REGEX.search(line):
                detected_style = 2
                break

    # Parse file lines according to detected style
    if detected_style == 1:
        for line in lines:
            match = STYLE1_REGEX.search(line)
            if match:
                name = match.group(1).strip()
                count = int(match.group(2))
                parsed_cuts.append({"name": name, "count": count, "raw_line": line.strip()})
                
    elif detected_style == 2:
        for line in lines:
            match = STYLE2_REGEX.search(line)
            if match:
                name = match.group(1).strip()
                count = int(match.group(3))
                parsed_cuts.append({"name": name, "count": count, "raw_line": line.strip()})
                
    else:
        print(f"Warning: Could not automatically determine cut-flow style for '{filepath}'", file=sys.stderr)

    return parsed_cuts


def compare_cut_flows(files: List[str]) -> None:
    """
    Parses a list of cut-flow logs, aligns them row-by-row, and prints
    a comparative overview highlighting any errors or discrepancies.
    """
    if len(files) < 2:
        print("Error: Please provide at least two cut-flow files to compare.", file=sys.stderr)
        return

    # Parse all provided log data
    all_file_data = {}
    max_rows = 0
    for fp in files:
        base_name = os.path.basename(fp)
        cuts = parse_cut_flow_file(fp)
        all_file_data[base_name] = cuts
        if len(cuts) > max_rows:
            max_rows = len(cuts)

    if not all_file_data or all(len(c) == 0 for c in all_file_data.values()):
        print("Error: No valid cut data could be extracted from the specified files.", file=sys.stderr)
        return

    # Dynamic layout calculation based on file names
    file_headers = list(all_file_data.keys())
    
    print("\n" + "=" * 100)
    print(" CUT-FLOW DISCREPANCY COMPARISON REPORT")
    print("=" * 100)
    
    # Header format configuration
    header_str = f"{'Row':<5} | "
    for fh in file_headers:
        header_str += f"{fh[:22]:<22} | {'Count':<12} | "
    header_str += "Status"
    print(header_str)
    print("-" * len(header_str))

    mismatch_count = 0

    # Process alignment row-by-row
    for idx in range(max_rows):
        row_str = f"{idx:<5} | "
        row_counts = []
        
        for fh in file_headers:
            file_cuts = all_file_data[fh]
            if idx < len(file_cuts):
                cut_item = file_cuts[idx]
                name = cut_item["name"]
                count = cut_item["count"]
                row_counts.append(count)
                
                # Truncate long cut names gracefully for presentation display
                display_name = name[:20] + ".." if len(name) > 22 else name
                row_str += f"{display_name:<22} | {count:<12} | "
            else:
                row_str += f"{' [NO MORE CUTS]':<22} | {'-':<12} | "
                row_counts.append(None)

        # Status validation checks
        # Verify if all parsed values on this row are identical
        if len(set(row_counts)) == 1 and None not in row_counts:
            status = "✅ OK"
        else:
            status = "❌ MISMATCH"
            mismatch_count += 1

        row_str += status
        print(row_str)

    print("-" * len(header_str))
    print(f"Total Cuts Evaluated: {max_rows}")
    if mismatch_count == 0:
        print("🎉 SUCCESS: All compared cut-flows are perfectly identical!")
    else:
        print(f"⚠️  ALERT: Found {mismatch_count} rows with errors or discrepancies.")
    print("=" * len(header_str) + "\n")


if __name__ == "__main__":
    # Example execution configuration block
    # Usage: python compare_cuts.py cutflow_style1.txt cutflow_style2.txt
    if len(sys.argv) > 1:
        compare_cut_flows(sys.argv[1:])
    else:
        print("--- CUT FLOW COMPARISON DEMO MODE ---")
        print("No command-line files provided. Generating inline samples to demonstrate accuracy...")
        
        # Creating demonstration sample files matching user descriptions
        sample1_content = """
        Cut  0:               All    4141125:    100        100
        Cut  1:           a_track    3830956:   92.5      92.51
        Cut  2: a_converged_track    3830825:    100      92.51
        Cut  3:     has_front_seg    3830816:    100      92.51
        """
        
        sample2_content = """
        ================== Total Cut Flow =======================
                                 Cut Group  Events Passing  Absolute [%]  Relative [%]                                                              Description
                             No cuts   N/A         4141125    100.000000    100.000000                                                     No selection applied
                         has_a_track   NaN         3830956     92.510030     92.510030                 Event has at least one reconstructed track (ntracks > 0)
                       is_good_track   NaN         3830825     92.506867     99.996580                              Track quality: status >= 0 and goodfit != 0
                   has_trk_front_seg   NaN         3830810     92.506650     99.999765                              Track has >=1 segment intersecting TT_Front
        """
        
        fn1, fn2 = "demo_style1.txt", "demo_style2.txt"
        with open(fn1, "w") as f: f.write(sample1_content)
        # Note: Intentionally modified row 3 count to 3830810 in sample 2 to show validation error catching
        with open(fn2, "w") as f: f.write(sample2_content)
        
        compare_cut_flows([fn1, fn2])
        
        # Clean up demo files safely
        for fn in [fn1, fn2]:
            if os.path.exists(fn):
                os.remove(fn)
