#! /bin/bash

DATASET=$1
if [ ! -f log/out_${DATASET}_thread_0.log ]; then
    echo "No ${DATASET} log files found!"
    exit 1
fi

awk '
/Cut *[0-9]+/ {
    # Extract the cut label (e.g., "Cut  0" or "Cut 100")
    match($0, /Cut *[0-9]+/)
    cut_label = substr($0, RSTART, RLENGTH)

    # Extract the specific numeric ID of the cut
    match(cut_label, /[0-9]+/)
    cut_id = substr(cut_label, RSTART, RLENGTH) + 0

    if (cut_id > max_cut) max_cut = cut_id

    # Split by the colon to separate "  Cut  X" from everything else
    split($0, parts, ":")

    # Strip leading/trailing spaces to isolate the name and the trailing metrics
    # Rest of string format: "     is_reco_electron       6391 (     63.91%)"
    match(parts[2], /[^ ].*( )/)
    content = substr(parts[2], RSTART, RLENGTH)

    # Extract the numerical values at the end of the line
    split(parts[2], end_nums, " ")

    # The count is the first number after the name, and the percentage is inside the parenthesis
    # end_nums[1] is the text name, end_nums[2] is the raw count
    name[cut_id] = end_nums[1]
    count[cut_id] += end_nums[2]

    # Track maximum name length for dynamic column layout alignment
    if (length(end_nums[1]) > max_name_len) max_name_len = length(end_nums[1])
}
END {
    id_width = length(max_cut)

    for (i = 0; i <= max_cut; i++) {
        if (i in count) {
            pct = (count[i] / count[0]) * 100

            # Constructs a format string using the dynamic id_width and max_name_len
            fmt = sprintf("  Cut %" id_width "d: %" max_name_len "s %10d ( %10.4g%%)\n", i, name[i], count[i], pct)
            printf(fmt)
        }
    }
}' log/out_${DATASET}_thread_*.log
