import pybedtools
import pandas as pd
from collections import defaultdict
import argparse
import os
import csv

argparser = argparse.ArgumentParser(description="Calculate statistics for low coverage regions in a BED file.")
argparser.add_argument("-l", "--low_cov_file", type=str, required=True, help="Path to the low coverage BED file.")
argparser.add_argument("-o", "--output_stat_file", type=str, help="Path to the output statistics file. If not provided, will be saved in the same directory as the input file.")

args = argparser.parse_args()
low_cov_file = args.low_cov_file
output_stat_file = args.output_stat_file
#low_cov_dir = os.path.dirname(low_cov_file)
#output_stat_file = os.path.join(low_cov_dir, "low_cov_stats.txt")
#low_cov_file = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/QC_BAM/Low_coverage_results/56001808050343/per-base-no-cov.bed.gz"  # Your low coverage BED file

special_regions_file = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/QC_BAM/Low_coverage_results/heterochromoatin_annotation.txt"  # BED file for centromeres, telomeres, heterochromatin

total_length_chromosome = {"1":249250621, "2":243199373, "3":198022430, "4":191154276, "5":180915260,
                           "6":171115067, "7":159138663, "8":146364022, "9":141213431, "10":135534747,
                           "11":135006516, "12":133851895, "13":115169878, "14":107349540, "15":102531392,
                           "16":90354753, "17":81195210, "18":78077248, "19":59128983, "20":63025520,
                           "21":48129895, "22":51304566, "X":155270560, "Y":59373566}

df = pd.read_csv(special_regions_file, sep="\t", comment='#')
df.columns = ['bin', 'chrom', 'chromStart', 'chromEnd', 'ix', 'n', 'size', 'type', 'bridge']
df_filtered = df[['chrom', 'chromStart', 'chromEnd', 'type']]
df_filtered['chrom'] = df_filtered['chrom'].str.replace("^chr", "", regex=True)


# === Load BED files using pybedtools ===
low_cov_bed = pybedtools.BedTool(low_cov_file)
special_bed = pybedtools.BedTool.from_dataframe(df_filtered)

low_cov_bed = low_cov_bed.sort().merge()
special_bed = special_bed.sort().merge()

# === Total stats for original low coverage regions ===
def total_bases(bedtool_obj):
    return sum([(int(r.end) - int(r.start)) for r in bedtool_obj])

# === Filter: remove any overlapping region ===
low_cov_filtered = low_cov_bed.intersect(special_bed, v=True)

# Filter out anything not in allowed chromosomes
allowed_chroms = [str(c) for c in range(1, 23)] + ["X", "Y"]
low_cov_filtered = low_cov_filtered.filter(lambda x: x.chrom in allowed_chroms).saveas()
low_cov_bed = low_cov_bed.filter(lambda x: x.chrom in allowed_chroms).saveas()

stats_original = defaultdict(lambda: {"regions": 0, "bases": 0})
stats_filtered = defaultdict(lambda: {"regions": 0, "bases": 0})

for r in low_cov_bed:
    stats_original[r.chrom]["regions"] += 1
    stats_original[r.chrom]["bases"] += int(r.end) - int(r.start) 

for r in low_cov_filtered:
    stats_filtered[r.chrom]["regions"] += 1
    stats_filtered[r.chrom]["bases"] += int(r.end) - int(r.start) 

print("Output file: ", output_stat_file)

with open(output_stat_file, "w") as f:
    f.write("Input bed file: " + low_cov_file + "\n")
    f.write("Statistics for low coverage regions\n")
    f.write("== Overall Statistics ==\n")
    f.write(f"Total regions: {len(low_cov_bed)}\n")
    f.write(f"Total base pairs: {total_bases(low_cov_bed):,}\n")
    f.write(f"Remaining regions: {len(low_cov_filtered)}\n")
    f.write(f"Remaining base pairs: {total_bases(low_cov_filtered):,}\n")
    f.write("\n== Per-Chromosome Statistics ==\n")
    f.write(f"{'Chr':<6}{'Orig. Regions':>15}{'Orig. BP':>15}{'Filt. Regions':>15}{'Filt. BP':>15}{'%Orig BP':>12}{'%Filt BP':>12}\n")
    f.write("-" * 90 + "\n")

    order = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
    all_chroms = sorted(set(stats_original.keys()) | set(stats_filtered.keys()), key=lambda x: order.index(x) if x in order else float('inf'))

    for chrom in all_chroms:
        if chrom not in total_length_chromosome:
            continue
        o = stats_original[chrom]
        f_stats = stats_filtered[chrom]
        total_len = total_length_chromosome[chrom]
        
        heterochromatin_length = sum([int(r.end) - int(r.start) for r in special_bed if r.chrom == chrom])
        filtered_len = total_len - heterochromatin_length 
        
        percent_orig = 100 * o['bases'] / total_len
        percent_filt = 100 * f_stats['bases'] / filtered_len

        f.write(f"{chrom:<6}{o['regions']:>15,}{o['bases']:>15,}{f_stats['regions']:>15,}{f_stats['bases']:>15,}{percent_orig:12.2f}%{percent_filt:12.2f}%\n")

sample_id = os.path.basename(low_cov_file).split("_")[0]  
output_table_path = output_stat_file.replace(".txt", "_plot.tsv")
print(output_table_path)
with open(output_table_path, "w", newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    writer.writerow(["Sample", "Chromosome", "Percent_Low_Coverage_Original", "Percent_Low_Coverage_Filtered"])
    
    for chrom in all_chroms:
        if chrom not in total_length_chromosome:
            continue
        o = stats_original[chrom]
        f_stats = stats_filtered[chrom]
        total_len = total_length_chromosome[chrom]
        heterochromatin_length = sum([int(r.end) - int(r.start) for r in special_bed if r.chrom == chrom])
        filtered_len = total_len - heterochromatin_length 

        percent_orig = 100 * o['bases'] / total_len if total_len > 0 else 0
        percent_filt = 100 * f_stats['bases'] / filtered_len if filtered_len > 0 else 0

        writer.writerow([sample_id, chrom, round(percent_orig, 2), round(percent_filt, 2)])