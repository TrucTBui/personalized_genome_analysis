# Run all low coverage stats that in /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/QC_BAM/Low_coverage_results
import os

# Find all low coverage stats files in the specified directory
directory = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/QC_BAM/Low_coverage_results"
low_cov_files = []
for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith("per-base-low-cov.bed.gz"):
                low_cov_files.append(os.path.join(root, file))
if not low_cov_files:
    print("No low coverage stats files found in the specified directory.")
else:
    for low_cov_file in low_cov_files:
        output_stat_file = low_cov_file.replace("per-base-low-cov.bed.gz", "low_cov_stats.txt")
        print(f"Processing {low_cov_file} and saving results to {output_stat_file}")
        
        # Call the low coverage stats script with the current file
        os.system(f"/home/b/buit/miniconda3/envs/HiWi/bin/python /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_bam_quality/low_cov_stats.py -l {low_cov_file} -o {output_stat_file}")