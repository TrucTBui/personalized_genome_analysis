#!/bin/bash

per_base_list="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/QC_BAM/list_per_base_files_all.txt"

# Check if the BAM list file exists
if [[ ! -f "$per_base_list" ]]; then
    echo "Error: Input file $per_base_list not found!"
    exit 1
fi

mkdir -p /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/logs/
rm -f /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/logs/run_filter_low_cov_*.out
rm -f /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/logs/run_filter_low_cov_*.err

# Loop through each BAM file in the list
while IFS= read -r per_base_file; do
    echo "submitted $per_base_file"

    if [[ ! -f "$per_base_file" ]]; then
        echo "Error: BAM file $per_base_file not found!"
        continue
    fi

    base_name=$(basename $(dirname "$per_base_file"))

    output_dir="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/QC_BAM/Low_coverage_results"
    out_prefix="${output_dir}/${base_name}/"
    out_file="${out_prefix}${base_name}_per-base-low-cov.bed.gz"

    mkdir -p "$out_prefix"

    qsub -N low-cov-$base_name -b y -q dell.q \
        -o /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/logs/run_filter_low_cov_$base_name.out \
        -e /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/logs/run_filter_low_cov_$base_name.err \
        -l vf=0.5G \
        /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/run_filter_low_cov_individual.sh "$per_base_file" "$out_file"

done < "$per_base_list"
