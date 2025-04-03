"""
version 28.03.2025
"""
import subprocess
import argparse
import os
from collections import Counter, defaultdict
import time
from intervaltree import IntervalTree, Interval
import re
import pandas as pd
import tempfile

parser = argparse.ArgumentParser()

parser.add_argument("-g", "--genome", type=str, required=True, help="file containing path to genome(s) in bam format")
parser.add_argument("-r", "--reference", type=str, required=True, help="reference genome in fasta format")
parser.add_argument("-l", "--location", type=str, required=True, help="genome location to extract reads")
parser.add_argument("-o", "--output", type=str, required=False, help="output folder for the result", default="./")
parser.add_argument("-p", "--person", type=str, required=True, help="which person in the family")

args = parser.parse_args()
genome_path = args.genome
ref_path = args.reference
output_path = args.output
person = args.person
location_path = args.location

def parse_genome(path):
    genomes = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            genomes.append(line.strip())
    return genomes

def create_mpileup(sam_file, region_bed):
    """
    Efficiently generates a DataFrame from a SAM file using samtools mpileup,
    writing directly to a CSV file and loading it with Pandas.

    Args:
        sam_file (str): Path to the SAM file.
        region_bed (str): A bed file containing genomic region in the format "chr   start end".


    Returns:
        pd.DataFrame: A DataFrame with columns: Chromosome, Position, Parsed_Bases
    """
    genome_basename = os.path.splitext(os.path.basename(sam_file))[0]

    # Create a temporary file to store the results of the mpileup
    with tempfile.NamedTemporaryFile(delete=False, mode="w+", suffix=".csv") as temp_file:
        temp_path = temp_file.name  # Get temp file path

        try:
            # Open the temp file for writing results
            with open(temp_path, "w") as out_file:
                # Open the BED file and iterate through each region
                with open(region_bed, "r") as bed_file:
                    for line in bed_file:
                        line = line.strip()
                        chromosome, start, end = line.split("\t")
                        region = f"{chromosome}:{start}-{end}"

                        # Run samtools mpileup for the current region
                        samtools_command = ["samtools", "mpileup", "-aa", sam_file, "-r", region]
                        result = subprocess.check_output(samtools_command, text=True, stderr=subprocess.DEVNULL)

                        out_file.write(result)

            # Now read the results into a Pandas DataFrame
            df = pd.read_csv(temp_path, sep="\t", header=None, usecols=[0, 1, 4],
                             names=["Chromosome", "Position", "Bases"])

            # Apply parsing function (assuming `format_frequencies` is defined)
            df[f"Frequency_{genome_basename}"] = df["Bases"].apply(lambda x: format_frequencies(x))
            df.drop(columns=["Bases"], inplace=True)  # Remove raw bases after processing

        except FileNotFoundError:
            print(f"Error: SAM file '{sam_file}' or BED file '{region_bed}' not found.")
            return pd.DataFrame(columns=["Chromosome", "Position", f"Frequency_{genome_basename}"])
        except subprocess.CalledProcessError as e:
            print(f"Error: samtools mpileup failed: {e}")
            return pd.DataFrame(columns=["Chromosome", "Position", f"Frequency_{genome_basename}"])

        finally:
            # Clean up the temporary file
            if os.path.exists(temp_path):
                os.remove(temp_path)

    return df

def format_frequencies(bases):
    """
    Parses bases, computes base frequencies, and formats them as 'A:3;T:2' strings.
    """
    counter = Counter(parse_mpileup_bases(bases))
    if len(counter) > 0:
        return ";".join(f"{key}:{value}" for key, value in counter.items())
    else:
        return "NA"


def parse_mpileup_bases(bases):
    """
    Parses the base string from samtools mpileup, including indels.

    Args:
        bases (str): The base string from samtools mpileup.

    Returns:
        list: A list of observed bases and indels at the position.
    """
    parsed_bases = []
    i = 0
    while i < len(bases):
        base = bases[i]
        # Match standard bases
        if base in "ATCGatcg":
            parsed_bases.append(base.upper())
            i += 1
        # Match insertions (+n[ATCG...])
        elif base == '+':
            match = re.match(r'\+(\d+)([ATCGatcg]+)', bases[i:])
            if match:
                length = int(match.group(1))
                insertion = match.group(2)[:length].upper() # extract only the number of letters specified by the number
                #parsed_bases.append(f"+{insertion}")
                i += len(match.group(1)) + length + 1 #increment by the length of the number, the insertion, and the + sign.
            else:
                i += 1  # handle unexpected +
        # Match deletions (-n[ATCG...])
        elif base == '-':
            match = re.match(r'\-(\d+)([ATCGNatcgn]+)', bases[i:])
            if match:
                length = int(match.group(1))
                deletion = match.group(2)[:length].upper() #extract only the number of letters specified by the number
                #parsed_bases.append(f"-{deletion}")
                i += len(match.group(1)) + length + 1 #increment by the length of the number, the deletion, and the - sign.
            else:
                i += 1  # handle unexpected -
        # match ^ and $ which are read quality and end of read symbols.
        elif base == '^' or base == '$':
            i += 1
        else:
            i += 1  # handle other characters
    return parsed_bases


def merge_dataframes(dfs):
    """
    Merges multiple DataFrames on 'Chromosome' and 'Position'.

    Args:
        dfs (list): List of DataFrames to merge.

    Returns:
        pd.DataFrame: Merged DataFrame.
    """
    for df in dfs:
        df["Chromosome"] = df["Chromosome"].astype(str)

    merged_df = dfs[0]  # Start with the first DataFrame

    for df in dfs[1:]:  # Loop over remaining DataFrames
        merged_df = pd.merge(merged_df, df, on=["Chromosome", "Position"], how="outer")

    return merged_df

def determine_final_base(freq_counts_per_replicate):
    combined_counts = Counter()
    for counts in freq_counts_per_replicate:
        combined_counts.update(counts)

    if not combined_counts:
        return 'N'  # No counts available

    sorted_bases = sorted(combined_counts.items(), key=lambda x: x[1], reverse=True)
    most_frequent_base = sorted_bases[0][0]
    second_frequent_base = sorted_bases[1][0] if len(sorted_bases) > 1 else None

    total_counts = sum(combined_counts.values())
    ratio1 = combined_counts[most_frequent_base] / total_counts
    ratio12 = (combined_counts[most_frequent_base] + (
        combined_counts[second_frequent_base] if second_frequent_base else 0)) / total_counts

    if ratio1 >= 0.75 or (total_counts < 15 and ratio1 >= 0.65):
        return most_frequent_base
    elif ratio12 > 0.65 and second_frequent_base:
        return f"{most_frequent_base}/{second_frequent_base}"
    else:
        return 'N'


def evaluate_consistency(freq_counts_per_replicate):
    # Check for any empty replicate counts
    if any(len(replicate_counts) == 0 for replicate_counts in freq_counts_per_replicate):
        return "no_counts"

    major_alleles_per_replicate = []
    for counts in freq_counts_per_replicate:
        total = sum(counts.values())
        sorted_alleles = sorted(counts.items(), key=lambda x: x[1], reverse=True)

        major_alleles = []
        cumulative_count = 0
        for base, count in sorted_alleles:
            major_alleles.append(base)
            cumulative_count += count
            if cumulative_count / total >= 0.70:
                break
        major_alleles_per_replicate.append(set(major_alleles))

    if len(major_alleles_per_replicate) == 0:
        return "no_counts"

    first_major = major_alleles_per_replicate[0]

    # If the first replicate's major alleles are a superset of the others, it's considered consistent
    for major in major_alleles_per_replicate[1:]:
        if not first_major.issubset(major) and not major.issubset(first_major):
            return "inconsistent"

    return "consistent"

def calc_coverage(freq_counts_per_replicate):
    total_coverage = 0
    for rep_count in freq_counts_per_replicate:
        for _,count in rep_count.items():
            total_coverage+=count
    return total_coverage


def process_replicates(df):
    freq_columns = df.filter(regex='^Frequency').columns
    final_bases = []
    consistencies = []
    total_coverage = []

    for _, row in df.iterrows():
        freq_counts_per_replicate = []
        for col in freq_columns:
            freq = row[col]
            if pd.isna(freq) or freq == "NA" or freq == "":
                continue

            replicate_counts = {}
            for allele in freq.split(";"):
                if allele:
                    base, count = allele.split(":")
                    replicate_counts[base] = int(count)
            freq_counts_per_replicate.append(replicate_counts)

        final_bases.append(determine_final_base(freq_counts_per_replicate))
        consistencies.append(evaluate_consistency(freq_counts_per_replicate))
        total_coverage.append(calc_coverage(freq_counts_per_replicate))

    df['Coverage'] = total_coverage
    df["Final_Base"] = final_bases
    df["Consistency"] = consistencies
    return df


def transform_location_input(path):
    """
    Transform the location information from ISAR database to the suitable format for samtools
    """
    locations= []
    locations_with_type = {}

    with open(path, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                fields = line.split()
                chromosome = fields[0]
                type = fields[2]
                #if type in {"exon", "gene", "transcript", "intron"}:
                if type in {"exon", "gene", "transcript"}:
                    continue

                if chromosome.startswith("chr"):
                    chromosome = chromosome[3:]
                start = int(fields[3])
                end = int(fields[4]) - 1  # End is exclusive

                locations.append(f"{chromosome}\t{start}\t{end}")

                location_string = f"{chromosome}:{start}-{end}"
                if location_string not in locations_with_type:
                    locations_with_type[location_string] = set()
                locations_with_type[location_string].add(type)
    merged_locations = merge_regions(locations)
    with open(f"{os.path.dirname(path)}/locations_transformed.bed", 'w') as w:
        for location in merged_locations:
            w.write(f"{location}\n")

    return locations, locations_with_type


def merge_regions(regions):
    chrom_trees = defaultdict(IntervalTree)

    # Populate IntervalTree for each chromosome
    for region in regions:
        chrom, start, end = region.split("\t")
        start, end = int(start), int(end)
        chrom_trees[chrom][start:end] = None  # No additional data needed, just the interval

    # Function to iteratively merge adjacent intervals
    def merge_adjacent(intervals):
        merged = []
        for iv in intervals:
            if merged and merged[-1][1] >= iv.begin -1:  # Merge overlapping or adjacent
                merged[-1] = (merged[-1][0], max(merged[-1][1], iv.end))
            else:
                merged.append((iv.begin, iv.end))
        return merged

    # Merge overlapping and adjacent intervals
    merged_regions = []
    for chrom, tree in chrom_trees.items():
        tree.merge_overlaps()  # Merge overlapping intervals

        # Sort and merge adjacent intervals iteratively
        sorted_intervals = sorted(tree)
        prev_merged = []
        merged = merge_adjacent(sorted_intervals)
        while prev_merged != merged:  # Keep merging until no further changes
            prev_merged = merged
            merged = merge_adjacent([Interval(start, end) for start, end in merged])

        # Store final merged intervals
        for start, end in merged:
            merged_regions.append(f"{chrom}\t{start}\t{end}")

    return merged_regions

def assign_and_merge_types(locations_with_type):
    """
    Expands location ranges into individual positions and assigns types.
    Merges types for positions appearing in multiple regions.

    :param locations_with_type: Dictionary {location (range) -> set of types}
    :return: Dictionary {position -> "Type1,Type2,..."} where position is "chromosome:pos"
    """
    position_to_type = defaultdict(set)

    for location, types in locations_with_type.items():
        chrom, range_part = location.split(":")
        start, end = map(int, range_part.split("-"))  # Convert range to integers

        for pos in range(start, end + 1):  # Expand range to individual positions
            position_to_type[(chrom, pos)].update(types)


    df = pd.DataFrame(
            [(chrom, pos, ",".join(sorted(types))) for (chrom, pos), types in position_to_type.items()],
            columns=["Chromosome", "Position", "Type"]
    )

    return df

def get_reference_bases(chromosome, start_pos, end_pos, reference_genome="/mnt/raidinput/input/own/ReferenceGenomes/human_g1k_v37.fasta.gz"):
    """Fetch the reference bases for a given region (start_pos to end_pos)."""
    cmd = f"samtools faidx {reference_genome} {chromosome}:{start_pos}-{end_pos}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    lines = result.stdout.strip().split("\n")[1:]  # Split and ignore the first header line

    # Join the sequence lines together
    ref_sequence = "".join(lines)  # Combine all parts of the sequence into one continuous string

    #print(ref_sequence)
    return ref_sequence

def integrate_reference_bases(df, reference_genome="/mnt/raidinput/input/own/ReferenceGenomes/human_g1k_v37.fasta.gz"):
    # Group by chromosome and get the min and max positions
    groups = df.groupby('Chromosome')['Position']
    reference_bases = {}

    for chrom, positions in groups:
        start_pos = positions.min()
        end_pos = positions.max()
        ref_sequence = get_reference_bases(chrom, start_pos, end_pos, reference_genome)

        # Store the reference sequence for the chromosome
        reference_bases[chrom] = ref_sequence

    def get_base_for_position(row):
        chrom = row['Chromosome']
        pos = row['Position']
        ref_sequence = reference_bases[chrom]
        # Return the base corresponding to the position
        return ref_sequence[pos - start_pos]  # Adjust the position within the region

    df['Reference_Base'] = df.apply(get_base_for_position, axis=1)
    return df


def identify_special_cases(df):
    """
    Identifies special cases in the DataFrame based on the following criteria:
    - When the reference base is not the final base
    - When the final base contains a haplotype (e.g., A/G)
    - When the final base is 'N'

    Returns:
        pd.DataFrame: DataFrame with an additional column "Special_Case" indicating the special cases.
    """

    # Create a new column to flag special cases
    def check_special_case(row):
        special_case = None

        # Check if the final base contains a haplotype (e.g., A/G)
        if "/" in row["Final_Base"]:
            special_case = "Haplotype"

        # Check if the final base is 'N'
        elif row["Final_Base"] == "N":
            special_case = "Unknown_Base"

        # Check if reference base is different from final base
        elif row["Reference_Base"] != row["Final_Base"]:
            special_case = "Ref!=Final_Base"


        return special_case if special_case else "-"

    # Apply the function to the DataFrame
    df["Special_Case"] = df.apply(check_special_case, axis=1)

    return df

def compute_stats(df, output_file):
    # Total number of positions
    total_positions = len(df)

    # Number of inconsistent positions
    inconsistent_positions = len(df[df["Consistency"] == "inconsistent"])

    # Number of no_count positions
    no_count_positions = len(df[df["Consistency"] == "no_counts"])

    # Number of haplotype positions (Final_Base contains "/")
    haplotype_positions = len(df[df["Final_Base"].str.contains("/", na=False)])

    # Number of positions where Reference_Base != Final_Base
    ref_not_final_base = len(df[df["Reference_Base"] != df["Final_Base"]])

    # Calculate percentages
    inconsistent_percentage = (inconsistent_positions / total_positions) * 100
    no_count_percentage = (no_count_positions / total_positions) * 100

    stats = {
        "Total Positions": int(total_positions),
        "Inconsistent Positions": int(inconsistent_positions),
        "No Count Positions": int(no_count_positions),
        "Haplotype Positions": int(haplotype_positions),
        "Reference != Final Base": ref_not_final_base,
        "Inconsistent Percentage": inconsistent_percentage,
        "No Count Percentage": no_count_percentage
    }
    with open(output_file, "w") as o:
        for key, val in stats.items():
            o.write(f"{key}: {val}\n")
    return stats

start_time = time.perf_counter()  # runtime measurement


locations, locations_with_type = transform_location_input(location_path)
region_bed = f"{os.path.dirname(location_path)}/locations_transformed.bed"

type_df = assign_and_merge_types(locations_with_type)
type_df['Person'] = person
type_df = integrate_reference_bases(type_df, reference_genome=ref_path)

genomes = parse_genome(genome_path)
rep_df_list =[]
for genome in genomes:
    rep_df_list.append(create_mpileup(genome, region_bed))

rep_df_merged = merge_dataframes(rep_df_list)
df_merged = merge_dataframes([type_df,rep_df_merged])
result = process_replicates(df_merged)
result = identify_special_cases(result)

if not os.path.exists(output_path):
    try:
        os.makedirs(output_path)
        print(f"Created output directory at {output_path}")
    except OSError as e:
        print(f"Error creating output folder '{output_path}': {e}")

result.to_csv(f"{output_path}/results.tsv", index=False)
compute_stats(result, f"{output_path}/stats.tsv")

end_time = time.perf_counter()
runtime = end_time - start_time
print(f"Runtime: {runtime:.5f} seconds")

