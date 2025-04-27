threshold=$1
if [ -z "$threshold" ]; then
    echo "Usage: $0 <threshold>"
    exit 1
fi
#!/bin/bash
source /etc/sge.sh

# Filter out low coverage regions (running on Grid Engine)
/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/run_filter_low_cov_all.sh "$threshold"

# Produce stattistics and plots
output_dir="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/QC_BAM/Low_coverage_results/threshold_$threshold"
/home/b/buit/miniconda3/envs/HiWi/bin/python /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/run_low_cov_stats_all.py "$output_dir"

# Combine all stats file into one csv file
cd "$output_dir" 

first_file=true
echo -n > low_cov_stats_plot_combined.tsv

for file in $(find . -name "*_low_cov_stats_plot.tsv" | sort); do
    if $first_file; then
        head -n 1 "$file" >> low_cov_stats_plot_combined.tsv  # keep the header
        tail -n +2 "$file" >> low_cov_stats_plot_combined.tsv
        first_file=false
    else
        tail -n +2 "$file" >> low_cov_stats_plot_combined.tsv
    fi
done

plot_dir="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/QC_BAM/Plots"

# Make heatmap
/home/b/buit/miniconda3/envs/HiWi/bin/python /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/plot_heat_map_low_cov.py \
-i "$output_dir"/low_cov_stats_plot_combined.tsv \
-o "$plot_dir"/low_cov_heatmap_threshold"$threshold".png \
--title "Low Coverage Statistics per Chromosome (<"$threshold"x Coverage)"

