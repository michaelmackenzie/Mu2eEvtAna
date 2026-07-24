#!/usr/bin/env bash

set -euo pipefail

DATASET=${1:?"Error: Dataset name must be provided as the first argument."}

if ! find log/ -name "out_${DATASET}_thread_*.log" -print -quit | grep -q .; then
    echo "No ${DATASET} log files found in log/" >&2
    exit 1
fi

awk '
# --- PHASE 1: PROCESS FILE BY FILE ---
# We use separate tracking variables to map out each file before saving globally.
# This prevents different log files from mangling each others structural signatures.

/Cut *[0-9]+/ {
    match($0, /Cut *[0-9]+/)
    cut_label = substr($0, RSTART, RLENGTH)
    match(cut_label, /[0-9]+/)
    cut_id = substr(cut_label, RSTART, RLENGTH) + 0

    split($0, parts, ":")
    sub(/^[ \t]+/, "", parts[2])
    split(parts[2], tokens, /[ \t]+/)
    cut_name = tokens[1]
    cut_count = tokens[2] + 0

    # If Cut 0 is found, it means a cut-flow block is starting in the file
    if (cut_id == 0) {
        file_block_idx++
    }

    # Store temporary file-level structure
    file_cut_name[file_block_idx, cut_id] = cut_name
    file_cut_count[file_block_idx, cut_id] = cut_count
    if (cut_id > file_max_cut[file_block_idx]) {
        file_max_cut[file_block_idx] = cut_id
    }
}

# Every time we finish reading a single thread log file, parse its blocks
ENDFILE {
    for (fb = 1; fb <= file_block_idx; fb++) {
        # Construct a unique structural signature string for this block
        # e.g., "0:All->1:de_track_count->2:is_reco_electron"
        sig = ""
        max_c = file_max_cut[fb]
        for (i = 0; i <= max_c; i++) {
            sig = sig i ":" file_cut_name[fb, i] "->"
        }

        # Check if we have encountered this exact cut-flow sequence before
        if (!(sig in sig_to_global_id)) {
            global_block_count++
            sig_to_global_id[sig] = global_block_count
            global_id_to_max_cut[global_block_count] = max_c

            # Save the names for this new structural path
            for (i = 0; i <= max_c; i++) {
                global_cut_name[global_block_count, i] = file_cut_name[fb, i]
            }
        }

        # Map the file block back to the consolidated global block ID
        g_id = sig_to_global_id[sig]

        # Accumulate metrics across files into the global pool
        for (i = 0; i <= max_c; i++) {
            global_cut_count[g_id, i] += file_cut_count[fb, i]

            # Max width tracking for printing layout alignment
            if (length(file_cut_name[fb, i]) > max_name_len) {
                max_name_len = length(file_cut_name[fb, i])
            }
        }
        if (global_id_to_max_cut[g_id] > global_max_cut) {
            global_max_cut = global_id_to_max_cut[g_id]
        }
    }

    # Reset temporary file-level arrays for the next thread log file
    delete file_cut_name
    delete file_cut_count
    delete file_max_cut
    file_block_idx = 0
}

# --- PHASE 2: PRINT THE RESULTS ---
END {
    id_width = length(global_max_cut)

    for (g_id = 1; g_id <= global_block_count; g_id++) {
        # Visual divider between different cut-flows
        if (g_id > 1) print ""

        base_count = global_cut_count[g_id, 0]
        if (base_count == 0) base_count = 1

        max_c = global_id_to_max_cut[g_id]
        prev_c = -1
        for (i = 0; i <= max_c; i++) {
            if ((g_id, i) in global_cut_count) {
                c_name = global_cut_name[g_id, i]
                c_val  = global_cut_count[g_id, i]
                pct    = (c_val / base_count) * 100
                if (prev_c < 0.) {
                  prev_c = c_val
                }
                rel_pct = (c_val / prev_c) * 100
                prev_c = c_val

                fmt = sprintf("  Cut %" id_width "d: %" max_name_len "s %10d: %6.3g %10.4g\n", i, c_name, c_val, rel_pct, pct)
                printf(fmt)
            }
        }
    }
}' log/out_${DATASET}_thread_*.log
