#!/usr/bin/env bash

set -euo pipefail

DATASET=${1:?"Error: Dataset name must be provided as the first argument."}

if ! find log/ -name "out_${DATASET}_thread_*.log" -print -quit | grep -q .; then
    echo "No ${DATASET} log files found in log/" >&2
    exit 1
fi

awk '
# --- PHASE 1: PROCESS FILE BY FILE ---
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

    # If Cut 0 is found, it means a new cut-flow block sequence is starting in this file
    if (cut_id == 0) {
        file_block_idx++
    }

    # Store file-level metrics using the sequential block index and cut index
    file_cut_name[file_block_idx, cut_id] = cut_name
    file_cut_count[file_block_idx, cut_id] = cut_count

    # Track the highest cut ID seen for this block across the file
    if (cut_id > file_max_cut[file_block_idx]) {
        file_max_cut[file_block_idx] = cut_id
    }
}

# Every time we finish reading a single thread log file, consolidate by block index
ENDFILE {
    # Keep track of the highest global block sequence index across all files
    if (file_block_idx > global_block_count) {
        global_block_count = file_block_idx
    }

    for (fb = 1; fb <= file_block_idx; fb++) {
        max_c = file_max_cut[fb]

        # Track the absolute maximum cut index seen globally for this block position
        if (max_c > global_id_to_max_cut[fb]) {
            global_id_to_max_cut[fb] = max_c
        }

        # Consolidate metrics directly into the matching global block index position
        for (i = 0; i <= max_c; i++) {
            if ((fb, i) in file_cut_count) {
                global_cut_count[fb, i] += file_cut_count[fb, i]

                # Dynamically fill/assert the master name for this slot if not set
                if (file_cut_name[fb, i] != "") {
                    global_cut_name[fb, i] = file_cut_name[fb, i]
                }

                # Max width tracking for pretty printing layout alignment
                if (length(global_cut_name[fb, i]) > max_name_len) {
                    max_name_len = length(global_cut_name[fb, i])
                }
            }
        }

        if (global_id_to_max_cut[fb] > global_max_cut) {
            global_max_cut = global_id_to_max_cut[fb]
        }
    }

    # Reset temporary file-level arrays for the next thread log file
    delete file_cut_name
    delete file_cut_count
    delete file_max_cut
    file_block_idx = 0
}

# --- PHASE 2: PRINT THE CONSOLIDATED RESULTS ---
END {
    id_width = length(global_max_cut)

    for (g_id = 1; g_id <= global_block_count; g_id++) {
        # Visual divider between sequential cut-flows
        if (g_id > 1) print ""

        base_count = global_cut_count[g_id, 0]
        if (base_count == 0) base_count = 1

        max_c = global_id_to_max_cut[g_id]
        prev_c = -1

        for (i = 0; i <= max_c; i++) {
            # Check if this index slot was populated by at least one job file
            if ((g_id, i) in global_cut_count) {
                c_name = global_cut_name[g_id, i]
                c_val  = global_cut_count[g_id, i]
                pct    = (c_val / base_count) * 100

                if (prev_c < 0) {
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
