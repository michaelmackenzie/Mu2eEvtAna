#!/bin/bash
# Merge Mu2eEvtAna output files from parallel processes
# Usage: merge_hist_files.sh <output_file> <input_list_file> [--force]

if [ $# -lt 2 ]; then
    echo "Usage: $0 <output_file> <input_list_file> [--force]"
    echo "  input_list_file: text file containing one .root file path per line"
    exit 1
fi

output="$1"
input_list="$2"
force=false

shift 2 || true
while [ $# -gt 0 ]; do
    case "$1" in
        --force|-f) force=true ;;
        *) echo "Unknown option: $1" ;;
    esac
    shift || true
done

if [ ! -f "$input_list" ]; then
    echo "Error: Input list file not found: $input_list"
    exit 1
fi

files=()
while IFS= read -r file || [ -n "$file" ]; do
    if [ -n "$file" ] && [ -f "$file" ]; then
        files+=("$file")
    fi
done < "$input_list"

if [ ${#files[@]} -eq 0 ]; then
    echo "Error: No valid input files found in $input_list"
    exit 1
fi

echo "Merging ${#files[@]} files into $output"
if [ "$force" = true ]; then
    hadd -f "$output" "${files[@]}"
else
    hadd "$output" "${files[@]}"
fi
