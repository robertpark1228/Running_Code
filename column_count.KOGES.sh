#!/bin/bash

# Check input
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

input_file="$1"
header=($(head -n 1 "$input_file"))
num_cols=${#header[@]}

# Known categorical values
known_cats=("0" "1" "2" "66666" "77777")

# Output header row
echo -e "Column\tType\tNA_Count\tKnown_Categories\tOther_Values\tAverage"

# Loop over columns (excluding ID)
for ((i=2; i<=num_cols; i++)); do
    col_name=${header[$((i-1))]}

    # Extract column values, skipping header
    col_data=$(awk -v col=$i 'NR > 1 { print $col }' "$input_file")

    # Count NA values
    na_count=$(echo "$col_data" | grep -c -E '^\s*$|NA')

    # Clean non-NA values
    clean_values=$(echo "$col_data" | grep -v -E '^\s*$|NA')

    # Count unique values
    unique_count=$(echo "$clean_values" | sort -u | wc -l)

    # If <= 10 unique values → categorical
    if [ "$unique_count" -le 10 ]; then
        type="Categorical"

        # Count known values
        known_counts=()
        for val in "${known_cats[@]}"; do
            count=$(echo "$clean_values" | awk -v v="$val" '$1 == v' | wc -l)
            [ "$count" -gt 0 ] && known_counts+=("${val}(${count})")
        done

        known_str=$(IFS=','; echo "${known_counts[*]}")

        # Detect other values
        other_vals=$(echo "$clean_values" | grep -v -x -e "${known_cats[0]}" -e "${known_cats[1]}" -e "${known_cats[2]}" -e "${known_cats[3]}" -e "${known_cats[4]}" | sort -u | tr '\n' ',' | sed 's/,$//')

        echo -e "${col_name}\t${type}\t${na_count}\t${known_str:--}\t${other_vals:--}\t-"

    else
        type="Continuous"

        avg=$(echo "$clean_values" | awk '{sum += $1; n++} END { if (n > 0) printf "%.2f", sum/n; else print "-"}')

        echo -e "${col_name}\t${type}\t${na_count}\t-\t-\t${avg}"
    fi
done
