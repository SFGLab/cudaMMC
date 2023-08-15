#!/bin/bash

# Default output directory
OUTPUT_DIR="chr"

# Function to display help message
display_help() {
    echo "Usage: $0 [-c chr11,chr13 or -c chr10-chr12] [-d output_directory]"
    echo "   -c: Define chromosomes either as comma-separated list or as a range."
    echo "   -d: Define the output directory (default is 'chr')."
    echo "   -h: Display this help message."
    exit 0
}

# Function to process the data file for given chromosomes
process_data_file() {
    local data_file="data.txt"
    local output_dir="$1"
    shift
    local chromosomes=("$@")

    mkdir -p "$output_dir"

    for chr in "${chromosomes[@]}"; do
        # Check if range is provided like chr10-chr12
        if [[ $chr == *-* ]]; then
            IFS="-" read -ra RANGE <<< "$chr"
            start_chr=$(echo "${RANGE[0]}" | awk -F "chr" '{print $2}')
            end_chr=$(echo "${RANGE[1]}" | awk -F "chr" '{print $2}')
            for num in $(seq "$start_chr" "$end_chr"); do
                awk -v chr="chr$num" '{if ($1 == chr) print $0}' "$data_file" > "$output_dir/data.txt.chr$num"
            done
        else
            awk -v chr="$chr" '{if ($1 == chr) print $0}' "$data_file" > "$output_dir/data.txt.$chr"
        fi
    done
}

# Parse input arguments
while getopts "c:d:h" option; do
    case "$option" in
        c) IFS="," read -ra CHROMOSOMES <<< "$OPTARG";;
        d) OUTPUT_DIR="$OPTARG";;
        h) display_help;;
    esac
done

# Check for provided chromosomes
if [ -z "${CHROMOSOMES}" ]; then
    echo "Please provide chromosomes. Use -h for help."
    exit 1
fi

# Main
process_data_file "$OUTPUT_DIR" "${CHROMOSOMES[@]}"