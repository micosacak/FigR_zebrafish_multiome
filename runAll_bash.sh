#!/bin/bash

# Input template file
TEMPLATE_FILE="01_jobXX_idxYY.sh"

# Check if start and end range and index range(s) are provided
if [ $# -ne 3 ]; then
    echo "Error: Please provide start and end job numbers and index range(s) (e.g., ./runAll_bash.sh 1 3 1:5 or ./runAll_bash.sh 1 3 '1:1 3:4 8:12')"
    exit 1
fi

START_JOB=$1
END_JOB=$2
INDEX_RANGES=$3

# Validate job number inputs
if ! [[ "$START_JOB" =~ ^[0-9]+$ ]] || ! [[ "$END_JOB" =~ ^[0-9]+$ ]]; then
    echo "Error: Start and end job numbers must be positive integers"
    exit 1
fi

if [ "$START_JOB" -gt "$END_JOB" ]; then
    echo "Error: Start job number must be less than or equal to end job number"
    exit 1
fi

# Parse index ranges (e.g., "1:5" or "1:1 3:4 8:12")
# Split by spaces to handle multiple ranges
IFS=' ' read -r -a ranges <<< "$INDEX_RANGES"
the_indices=()

for range in "${ranges[@]}"; do
    if ! [[ "$range" =~ ^[0-9]+:[0-9]+$ ]]; then
        echo "Error: Each range must be in format 'start:end' (e.g., 1:5 or 1:1 3:4 8:12)"
        exit 1
    fi
    START_IDX=$(echo "$range" | cut -d':' -f1)
    END_IDX=$(echo "$range" | cut -d':' -f2)
    if ! [[ "$START_IDX" =~ ^[0-9]+$ ]] || ! [[ "$END_IDX" =~ ^[0-9]+$ ]]; then
        echo "Error: Each range must contain positive integers"
        exit 1
    fi
    if [ "$START_IDX" -gt "$END_IDX" ]; then
        echo "Error: Start index must be less than or equal to end index in range $range"
        exit 1
    fi
    # Add indices from this range to the_indices
    the_indices+=($(seq "$START_IDX" "$END_IDX"))
done

# Remove duplicates and sort indices
the_indices=($(echo "${the_indices[@]}" | tr ' ' '\n' | sort -nu))

# Validate that indices array is not empty
if [ ${#the_indices[@]} -eq 0 ]; then
    echo "Error: No valid indices provided in the third argument"
    exit 1
fi

# Loop through the range of job numbers
for (( JOB_NUM=$START_JOB; JOB_NUM<=$END_JOB; JOB_NUM++ )); do
    for IDX in "${the_indices[@]}"; do
        # Zero-pad the index for consistent naming
        PADDED_IDX=$(printf "%02d" "$IDX")
        # New file name
        NEW_SCRIPT="01_job${JOB_NUM}_idx${PADDED_IDX}.sh"

        # Copy template and modify
        cp "$TEMPLATE_FILE" "$NEW_SCRIPT"
        sed -i "s/--job-name=01jobXX/--job-name=01job${JOB_NUM}_${PADDED_IDX}/" "$NEW_SCRIPT"
        sed -i "s/--output=01_jobXX_idxYY.txt/--output=01_job${JOB_NUM}_idx${PADDED_IDX}.txt/" "$NEW_SCRIPT"
        sed -i "s/--args value=YY/--args value=${JOB_NUM}/" "$NEW_SCRIPT"
        sed -i "s/nidxs=XX:XX/nidxs=${IDX}:${IDX}/" "$NEW_SCRIPT"
        sed -i "s/01_jobXX_idxYY.R/01_job${JOB_NUM}_idx${PADDED_IDX}.R/" "$NEW_SCRIPT"

        # Submit the job
        sbatch "$NEW_SCRIPT"
        echo "Generated and submitted $NEW_SCRIPT with job number $JOB_NUM and nidxs=${IDX}:${IDX}"
    done
done

echo "All jobs from $START_JOB to $END_JOB with indices $INDEX_RANGES have been submitted"
