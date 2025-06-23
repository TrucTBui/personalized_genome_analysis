import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Plot heatmap of low coverage statistics.")
parser.add_argument("-i", "--input", type=str, required=True, help="Input TSV file with low coverage statistics.")
parser.add_argument("-o", "--output", type=str, required=True, help="Output PNG file for the heatmap.")
parser.add_argument("-t", "--title", type=str, default="Low Coverage Heatmap", help="Title for the heatmap.")
args = parser.parse_args()

# Load data
df = pd.read_csv(args.input, sep="\t")
# For testing, use the hardcoded path
#df = pd.read_csv("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/QC_BAM/Low_coverage_results/low_cov_stats_plot_combined.tsv", sep="\t")

# Pivot for heatmap
pivot = df.pivot(index="Chromosome", columns="Sample", values="Percent_Low_Coverage_Filtered")

# Sort chromosomes correctly
chrom_order = [str(i) for i in range(1, 23)] + ["X", "Y"]
pivot = pivot.reindex(chrom_order)

# Plot heatmap
plt.figure(figsize=(12, 8))
sns.heatmap(pivot, cmap="YlOrRd", annot=True, vmax = 3 , fmt="g", linewidths=0.5, cbar_kws={'label': '% Low Coverage'})
plt.title(args.title)
plt.xticks(rotation=45, ha="right")
plt.xlabel("Sample")
plt.ylabel("Chromosome")
plt.tight_layout()
plt.show()
plt.savefig(args.output, dpi=300)
plt.close()
