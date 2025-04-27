source /etc/sge.sh

bam_list="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genomes_all.txt"

# Check if the file exists
if [ ! -f "$bam_list" ]; then
    echo "Error: File $bam_list does not exist."
    exit 1
fi

mkdir -p /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/logs/
rm -f /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/logs/run_mapq_*.out
rm -f /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/logs/run_mapq_*.err

job_ids=()

# Loop through each BAM file in the list to calculate MAPQ
while IFS= read -r bam_file; do
    echo "submitted $bam_file"

    if [[ ! -f "$bam_file" ]]; then
        echo "Error: BAM file $bam_file not found!"
        continue
    fi

    base_name=$(basename "$bam_file" .bam)

    output_dir="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/QC_BAM/MapQ"

    mkdir -p "$output_dir"

    qsub_command="qsub -N mapq-$base_name -b y -q dell.q \
        -o /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/logs/run_mapq_$base_name.out \
        -e /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/logs/run_mapq_$base_name.err \
        -l vf=0.5G \
        /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/run_mapq_bam.sh "$bam_file" "$output_dir/$base_name.mapq.txt""  

    job_id=$(eval "$qsub_command" | awk '{print $3}') 
    job_ids+=("$job_id") 
 

done < "$bam_list"

for job_id in "${job_ids[@]}"; do
    # Use qstat to check job status.  Loop until job is not found.
    while qstat -j "$job_id" > /dev/null 2>&1; do
        sleep 20 
    done
    echo "Job $job_id finished"
done

# Combine all MAPQ files into one
combined_file="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/QC_BAM/MapQ/all_mapq.txt"
rm -f "$combined_file"
for mapq_file in /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/QC_BAM/MapQ/*.mapq.txt; do
    if [[ -f "$mapq_file" ]]; then
        cat "$mapq_file" >> "$combined_file"
    fi
done

/home/b/buit/miniconda3/envs/HiWi/bin/python /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/plot_mapq.py "$combined_file" "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/QC_BAM/Plots/mapq_plot.png"