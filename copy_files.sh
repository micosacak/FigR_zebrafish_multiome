#!/bin/bash
# Check if there are 3 arguments
# arg 1 is an integer to start
# arg 2 is an integer to end
# arg 3 is one or more ranges of indices (e.g., 1:5 or 1:1 3:4 8:12)
if [ $# -ne 3 ]; then
    echo "Error: Please provide start and end job numbers and index range(s) (e.g., cp_file.sh 1 3 1:5 or cp_file.sh 1 3 '1:1 3:4 8:12')"
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

echo "Indices: ${the_indices[@]}"
original_file="01_jobXX_idxYY.R"

for (( JOB_NUM=$START_JOB; JOB_NUM<=$END_JOB; JOB_NUM++ )); do
    for IDX in "${the_indices[@]}"; do
        # Zero-pad the index to ensure consistent naming (e.g., 01, 02, etc.)
        PADDED_IDX=$(printf "%02d" "$IDX")
        new_file="01_job${JOB_NUM}_idx${PADDED_IDX}.R"
        cp "$original_file" "$new_file"
        echo "Copied $original_file to $new_file"
    done
done
