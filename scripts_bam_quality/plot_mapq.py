import pandas as pd
import matplotlib.pyplot as plt
import sys

# Read input arguments
mapq_file = sys.argv[1]
output_plot = sys.argv[2]

if len(sys.argv) != 3:
    print("Usage: python plot_mapq.py <mapq_file> <output_plot>")
    sys.exit(1)

# Load data
df = pd.read_csv(mapq_file, sep="\t", names=["Sample", "MAPQ", "Count"])

# Assign MAPQ group
def categorize_mapq(q):
    if q >= 1:
        return "MAPQ >= 1"
    else:
        return "MAPQ = 0"

df["MAPQ_Group"] = df["MAPQ"].apply(categorize_mapq)

# Aggregate by sample and MAPQ group
grouped = df.groupby(["Sample", "MAPQ_Group"], as_index=False)["Count"].sum()

# Calculate percentage per sample
grouped["Percentage"] = grouped["Count"] / grouped.groupby("Sample")["Count"].transform("sum") * 100

# Plot: stacked bar chart
pivot = grouped.pivot(index="Sample", columns="MAPQ_Group", values="Percentage").fillna(0)
pivot = pivot[["MAPQ = 0","MAPQ >= 1"]]  # consistent order

pivot.plot(kind="bar", stacked=True, figsize=(14, 7), colormap="Set2", edgecolor="black")

plt.ylabel("Percentage of Reads")
plt.title("MAPQ Quality Distribution per Sample")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.legend(title="MAPQ Group")
plt.savefig(output_plot)
