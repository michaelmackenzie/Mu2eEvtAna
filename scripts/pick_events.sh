#! /bin/bash

shopt -s extglob

EVTIDS=$1
DATASET=$2
if [ ! -f ${EVTIDS} ]; then
    echo "No file ${EVTIDS}!"
    exit 1
fi
if [[ "${DATASET}" == "" ]]; then
    echo "No dataset given!"
    exit 1
fi

FILE_LIST="${DATASET}_files.txt"
> "${FILE_LIST}" # Clear it if it already exists
while IFS= read -r line; do
    echo "Processing: $line";
    cleaned_line="${line//:+([[:space:]])//}"
    echo ${cleaned_line}
    pickEvent -e ${DATASET} "${cleaned_line}"

    # Parse the expected file name out of your line structure (e.g. 1430/42011/206664 -> 1430_42011_206664)
    # Using bash expansion to swap the forward slashes to underscores
    file_suffix="${cleaned_line//\//_}"
    expected_file="${DATASET}_${file_suffix}.art"

    # Check if the file was successfully produced
    if [ -f "${expected_file}" ]; then
        echo "--> Successfully produced: ${expected_file}"
        # Append to our merge tracking list
        echo "${expected_file}" >> "${FILE_LIST}"
    else
        echo "--> [WARNING] Expected file not found: ${expected_file}"
    fi
done < ${EVTIDS}

# Check if we actually recorded any files before running the merge
if [ -s "${FILE_LIST}" ]; then
    echo "--------------------------------------------------"
    echo "Merging $(wc -l < "${FILE_LIST}") files..."

    # Combine the files using the -o/--output flag for the merged name
    art -c /dev/null -S "${FILE_LIST}" -o "${DATASET}.evtids.art"

    echo "Merge complete! Created: ${DATASET}.evtids.art"
else
    echo "Error: No .art files were generated, skipping merge step."
fi

# Optional cleanup: remove the temporary file list text file
# rm "${FILE_LIST}"
