#!/bin/bash

# Usage: ./script.sh .per-base.bed.gz per-base-no-cov.bed.gz

input=$1
output=$2
threshold=${3:-2}

echo "Input file: $input"
echo "Output file: $output"
# Check if the input file exists
if [[ ! -f "$input" ]]; then
    echo "Error: Input file $input not found!"
    exit 1
fi

zcat "$input" | awk '$4 + 0 < '$threshold'' | gzip > "$output"