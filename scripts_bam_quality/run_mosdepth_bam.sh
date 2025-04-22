#!/bin/bash

source /home/b/buit/miniconda3/etc/profile.d/conda.sh
conda_env_name="HiWi"
conda activate "$conda_env_name"

if [[ $? -ne 0 ]]; then
  echo "Error: Failed to activate conda environment '$conda_env_name'"
  exit 1
fi

output_prefix="$1"
bam_file="$2"

echo "Running mosdepth with:"
echo "  Output prefix: $output_prefix"
echo "  BAM file:      $bam_file"

# Actually run the command
mosdepth -x -t 8 "$output_prefix" "$bam_file"
#touch "${output_prefix}.dummy.test"