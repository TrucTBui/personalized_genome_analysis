#!/bin/bash

bam_list="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genomes_all.txt"
#bam_list="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/test_input_genome.txt"

# Check if the BAM list file exists
if [[ ! -f "$bam_list" ]]; then
    echo "Error: Input file $bam_list not found!"
    exit 1
fi

mkdir -p /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/logs/
rm -f /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/logs/run_*.out
rm -f /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/logs/run_*.err

# Loop through each BAM file in the list
while IFS= read -r bam_file; do
    echo "submitted $bam_file"

    if [[ ! -f "$bam_file" ]]; then
        echo "Error: BAM file $bam_file not found!"
        continue
    fi

    base_name=$(basename "$bam_file" .bam)

    output_dir="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/QC_BAM"
    out_prefix="${output_dir}/${base_name}/"

    mkdir -p "$out_prefix"

    qsub -N $base_name -b y -q dell.q \
        -o /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/logs/run_$base_name.out \
        -e /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/logs/run_$base_name.err \
        -l vf=0.5G \
        /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/run_mosdepth_bam.sh "$out_prefix" "$bam_file"

done < "$bam_list"
