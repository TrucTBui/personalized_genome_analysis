"""
version 05.05.2025
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
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import gzip

parser = argparse.ArgumentParser()

parser.add_argument("-g", "--genome", type=str, required=True, help="file containing path to genome(s) in bam format")
parser.add_argument("-r", "--reference", type=str, required=True, help="reference genome in fasta format")
parser.add_argument("-l", "--location", type=str, required=True, help="genome location to extract reads (gtf file)")
parser.add_argument("-o", "--output", type=str, required=False, help="output folder for the result")
parser.add_argument("-p", "--person", type=str, required=True, help="which person in the family to be processed")

args = parser.parse_args()
genome_path = args.genome
ref_path = args.reference
output_path = args.output
person = args.person
location_path = args.location
gene_id = os.path.basename(location_path).split(".")[0]
MALES = ["child", "father", "grandfather_father", "grandfather_mother"]

def parse_genome(path):
    genomes = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            genomes.append(line.strip())
    return genomes

def create_mpileup(sam_file, locations):
    """
    Generates a DataFrame from a SAM file using samtools mpileup,
    writing directly to a CSV file and loading it with Pandas.

    Args:
        sam_file (str): Path to the SAM file.
        locations (list): List of locations


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
                for location in locations:
                    samtools_command = ["samtools", "mpileup", "-aa", "-B" , sam_file, "-r", location, "-f", ref_path]
                    result = subprocess.check_output(samtools_command, text=True, stderr=subprocess.DEVNULL)

                    out_file.write(result)

            # Now read the results into a Pandas DataFrame
            df = pd.read_csv(temp_path, sep="\t", header=None, usecols=[0, 1, 2, 4],
                             names=["Chromosome", "Position", "ReferenceBase" , "Bases"])

            # Apply parsing function, use the ReferenceBase to get the reference base
            df[f"Frequency"] = df.apply(lambda row: parse_mpileup_bases(row["Bases"], row["ReferenceBase"]), axis=1)
            # Format the frequencies to dictionary format
            df[f"Frequency_{genome_basename}"] = df[f"Frequency"].apply(count_frequencies)
            # Annotate the deletions
            df = annotate_deletions(df, genome_basename)
            df.drop(columns=["ReferenceBase","Bases", "Frequency"], inplace=True)  # Remove raw bases after processing

        except FileNotFoundError:
            print(f"Error: SAM file '{sam_file}' or BED file '{locations}' not found.")
            return pd.DataFrame(columns=["Chromosome", "Position", f"Frequency_{genome_basename}"])
        except subprocess.CalledProcessError as e:
            print(f"Error: samtools mpileup failed: {e}")
            return pd.DataFrame(columns=["Chromosome", "Position", f"Frequency_{genome_basename}"])

        finally:
            # Clean up the temporary file
            if os.path.exists(temp_path):
                os.remove(temp_path)

    return df

def count_frequencies(bases):
    # Count the occurrences of each base and return dictionary in the format "base:count"
    if not bases:
        return {}
    allele_count={}
    for base in bases:
        if base in allele_count:
            allele_count[base] += 1
        else:
            allele_count[base] = 1
    return allele_count

def format_frequencies(frequencies):
    """
    Formats the frequencies dictionary to a string representation.
    """
    if not frequencies:
        return "NA"
    formatted = []
    for base, count in frequencies.items():
            formatted.append(f"{base}:{count}")
    return ";".join(formatted)

def annotate_deletions(df, genome_basename):
    # Go through each row and annotate deletions in the Frequency column. 
    # If there is a "-" in the Frequency column, it means there is a deletion.
    # Add "D" to the Frequency column(s) that contain deletions.
    freq_col_name = f"Frequency_{genome_basename}"
    if freq_col_name not in df.columns:
        print(f"Warning: Frequency column '{freq_col_name}' not found. Skipping deletion marking.")
        return df

    deleted_genomic_positions = {}
    anchor_positions = set()  # To track anchor positions

    # Identify all genomic sites that are deleted based on anchor base reports
    for index, row in df.iterrows():
        # anchor_pos is the reference position *before* the deleted bases
        anchor_pos = row['Position'] 
        
        frequency_data = row[freq_col_name]
        #print(frequency_data)
        
        # Skip if frequency_data is not a dictionary (e.g., already 'NA' or 'D')
        if not isinstance(frequency_data, dict):
            continue

        for event_key in frequency_data.keys():
            if event_key.startswith('-'):
                #print(event_key)
                # event_key is like "-GCA" (meaning GCA deleted *after* anchor_pos)
                deleted_sequence = event_key[1:] # e.g., "GCA"
                deletion_length = len(deleted_sequence)
                anchor_positions.add(anchor_pos)  # Add anchor position to the set

                # Determine the genomic positions of the deleted bases
                for i in range(deletion_length):
                    # The first deleted base is at anchor_pos + 1, second at anchor_pos + 2, etc.
                    deleted_ref_pos = anchor_pos + 1 + i
                    if deleted_ref_pos not in deleted_genomic_positions:
                        deleted_genomic_positions[deleted_ref_pos] = frequency_data[event_key]
                    else:
                        deleted_genomic_positions[deleted_ref_pos] += frequency_data[event_key]


    for idx in df.index: # Use df.index to iterate over indices
        current_pos = df.loc[idx, 'Position']
        if current_pos in deleted_genomic_positions:
            # Add "D" to the Frequency column dictionary with the number of occurrences of the deletion
            current_dict = df.at[idx, freq_col_name]
            if isinstance(current_dict, dict):
                current_dict["D"] = deleted_genomic_positions[current_pos]
            df.at[idx, freq_col_name] = current_dict
            #print(df.at[idx, freq_col_name])
    
    for anchor_pos in anchor_positions:
        anchor_index = df[df['Position'] == anchor_pos].index[0]
        current_dict = df.at[anchor_index, freq_col_name]
        # Delete all the entries starting with "-" from the row of anchor positions
        if isinstance(current_dict, dict):
            # Remove all keys that start with "-"
            keys_to_remove = [key for key in current_dict.keys() if key.startswith('-')]
            for key in keys_to_remove:
                del current_dict[key]

            
    return df


def parse_mpileup_bases(bases, ref):
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
        if base in ['.', ',']:  # matches to the reference
            # Append reference base
            parsed_bases.append(ref)
            i += 1
        elif base in 'ATCGatcg': 
            parsed_bases.append(base.upper())
            i += 1
        # Match insertions (+n[ATCG...])
        elif base == '+':
            match = re.match(r'\+(\d+)([ATCGatcg]+)', bases[i:])
            if match:
                length = int(match.group(1))
                insertion = match.group(2)[:length].upper() # extract only the number of letters specified by the number
                parsed_bases.append(f"{ref}{insertion}")
                # extract the base right before the insertion
                last_base = bases[i-1].upper() if i > 0 else None
                #print(last_base)
                #if last_base and last_base in ['.', ',']:
                #    parsed_bases.remove(ref)
                i += len(match.group(1)) + length + 1 #increment by the length of the number, the insertion, and the + sign.

            else:
                i += 1  # handle unexpected +
            
        # Match deletions (-n[ATCG...])
        elif base == '-':
            match = re.match(r'\-(\d+)([ATCGNatcgn]+)', bases[i:])
            if match:
                length = int(match.group(1))
                deletion = match.group(2)[:length].upper() #extract only the number of letters specified by the number
                parsed_bases.append(f"-{deletion}")
                i += len(match.group(1)) + length + 1 #increment by the length of the number, the deletion, and the - sign.
            else:
                i += 1  # handle unexpected -
        else:
            i += 1  # handle other characters
    #print(parsed_bases)    
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
    if total_counts < 2:
        return 'N'  # Not enough counts to determine a base
    ratio1 = combined_counts[most_frequent_base] / total_counts
    ratio12 = (combined_counts[most_frequent_base] + (
        combined_counts[second_frequent_base] if second_frequent_base else 0)) / total_counts

    if (total_counts > 50 and ratio1 > 0.85) or (total_counts <= 50 and total_counts > 15 and ratio1 > 0.8) or (total_counts <= 15 and ratio1 >= 0.75):
        return most_frequent_base
    elif ratio12 > 0.70 and second_frequent_base:
        return f"{most_frequent_base}/{second_frequent_base}"
    else:
        return 'N'


def evaluate_consistency(freq_counts_per_replicate):
    # Check for any empty replicate counts
    #if any(len(replicate_counts) == 0 for replicate_counts in freq_counts_per_replicate):
    #    return "no_counts"

    major_alleles_per_replicate = []
    total_counts = 0
    for counts in freq_counts_per_replicate:
        total = sum(counts.values())
        total_counts += total
        sorted_alleles = sorted(counts.items(), key=lambda x: x[1], reverse=True)

        major_alleles = []
        cumulative_count = 0
        for base, count in sorted_alleles:
            major_alleles.append(base)
            cumulative_count += count
            if cumulative_count / total >= 0.70:
                break
        major_alleles_per_replicate.append(set(major_alleles))

    if len(major_alleles_per_replicate) == 0 or total_counts < 2:
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
            for allele, count in dict(freq).items():
                replicate_counts[allele] = int(count)
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
                if type in {"gene", "transcript"}:
                    continue

                if chromosome.startswith("chr"):
                    chromosome = chromosome[3:]
                start = int(fields[3])
                end = int(fields[4])   # End is exclusive

                locations.append(f"{chromosome}\t{start}\t{end}")

                location_string = f"{chromosome}:{start}-{end}"
                if location_string not in locations_with_type:
                    locations_with_type[location_string] = set()
                locations_with_type[location_string].add(type)
    merged_locations = merge_regions(locations)

    #with open(f"{os.path.dirname(path)}/locations_transformed.bed", 'w') as w:
    #    for location in merged_locations:
    #        w.write(f"{location}\n")

    return merged_locations, locations_with_type


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
            merged_regions.append(f"{chrom}:{start}-{end}")

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
            [(gene_id, chrom, pos, ";".join(sorted(types))) for (chrom, pos), types in position_to_type.items()],
            columns=["Gene", "Chromosome", "Position", "Type"]
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

        # Adjust the position within the region
        return ref_sequence[pos - start_pos]  # Make sure this returns a single value (a base)

    # Ensure the apply method returns a single column
    df['Reference_Base'] = df.apply(lambda row: get_base_for_position(row), axis=1)
    return df


def identify_special_cases(df):
    """
    Identifies special cases in the DataFrame based on the following criteria:
    - When the reference base is not the final base
    - When the final base is heterozygous (e.g., A/G)
    - When the final base is 'N'

    Returns:
        pd.DataFrame: DataFrame with an additional column "Special_Case" indicating the special cases.
    """

    # Create a new column to flag special cases
    def check_special_case(row):
        special_case = []

        if "/" in row["Final_Base"]:
            special_case.append("Heterozygous")
        else:
            if row["Final_Base"] == "N":
                special_case.append("Unknown_Base")
            elif row["Final_Base"] != row["Reference_Base"]:
                special_case.append("Homozygous!=Ref")

        if "D" in row["Final_Base"]:
            special_case.append("Deletion")
        bases = row["Final_Base"].split("/")
        # if there is a base longer than 1, it is an insertion
        if any(len(base) > 1 for base in bases):
            special_case.append("Insertion")


        return ";".join(special_case) if special_case else "-"

    # Apply the function to the DataFrame
    df["Special_Case"] = df.apply(check_special_case, axis=1)

    return df

def identify_ambigious_positions(df):
    """
    Filters for ambiguous positions from the DataFrame.
    Ambiguous positions are those where the final base is 'N' or the base calling is not consistent.
    """
    ambigious_positions = df[df["Final_Base"].isin(["N"]) | df["Consistency"].isin(["inconsistent", "no_counts"])]
    # If the pos is in the chrom Y and the special case is heterozygous, then it is also ambiguous (because there is only one Y chromosome)
    chrom = df["Chromosome"].values[0]
    person = df["Person"].values[0]
    if person in MALES:
        if chrom in ["Y","X"]:
            heterozygous_positions = df[df["Special_Case"] == "Heterozygous"]
            ambigious_positions = pd.concat([ambigious_positions, heterozygous_positions])

    return ambigious_positions

def identify_variant_positions(df):
    """
    Filters for variants positions from the DataFrame (Heterozygous + SNP).
    Heterozygous positions are those where the final base contains a '/'.
    """
    chrom = df["Chromosome"].values[0]
    person = df["Person"].values[0]

    if person in MALES and chrom in ["Y","X"]:
        variant_cases = ["Homozygous!=Ref", "Insertion", "Deletion"]
    else:
        variant_cases = ["Heterozygous", "Homozygous!=Ref", "Insertion", "Deletion"]

    escaped_patterns = [re.escape(p) for p in variant_cases]
    regex_pattern = '|'.join(escaped_patterns)
    
    not_na_mask = df["Special_Case"].notna()
    
    combined_conditions = pd.Series(False, index=df.index)

    if not_na_mask.any():
        relevant_series = df.loc[not_na_mask, "Special_Case"].astype(str)
        combined_conditions.loc[not_na_mask] = relevant_series.str.contains(regex_pattern, regex=True)
            
    return df[combined_conditions]

def compute_stats(df, output_file):
    # Total number of positions
    total_positions_all = len(df)

    # Number of inconsistent positions
    inconsistent_positions_all = len(df[df["Consistency"] == "inconsistent"])

    # Number of no_count positions
    no_count_positions_all = len(df[df["Consistency"] == "no_counts"])

    # Number of Heterozygous positions (Final_Base contains "/")
    heterozygous_positions_all = len(df[df["Final_Base"].str.contains("/", na=False)])

    # Number of positions where Reference_Base != Final_Base
    ref_not_final_base_all = len(df[(df["Reference_Base"] != df["Final_Base"]) & (df["Final_Base"] != "N")])
    
    # Calculate percentages
    inconsistent_percentage_all = (inconsistent_positions_all / total_positions_all) * 100
    no_count_percentage_all = (no_count_positions_all / total_positions_all) * 100

    stats_all = {
        "Total Positions": int(total_positions_all),
        "Inconsistent Positions": int(inconsistent_positions_all),
        "No Count Positions": int(no_count_positions_all),
        "Heterozygous Positions": int(heterozygous_positions_all),
        "Reference != Final Base": ref_not_final_base_all,
        "Inconsistent Percentage": round(inconsistent_percentage_all,4),
        "No Count Percentage": round(no_count_percentage_all,4)
    }

    df['Type_List'] = df['Type'].apply(lambda x: [t.strip() for t in x.split(';')])

    # Collect all unique gene types
    all_types = set(t for types in df['Type_List'] for t in types)

    stats_by_type = {}

    for gene_type in sorted(all_types):
        # Filter rows that contain the current gene type
        subset = df[df['Type_List'].apply(lambda lst: gene_type in lst)]

        if subset.empty:
            continue

        total_positions = len(subset)
        inconsistent_positions = len(subset[subset["Consistency"] == "inconsistent"])
        no_count_positions = len(subset[subset["Consistency"] == "no_counts"])
        heterozygous_positions = len(subset[subset["Final_Base"].str.contains("/", na=False)])
        ref_not_final_base = len(subset[subset["Reference_Base"] != subset["Final_Base"]])

        inconsistent_percentage = (inconsistent_positions / total_positions) * 100
        no_count_percentage = (no_count_positions / total_positions) * 100

        stats = {
            "Total Positions": int(total_positions),
            "Inconsistent Positions": int(inconsistent_positions),
            "No Count Positions": int(no_count_positions),
            "Heterozygous Positions": int(heterozygous_positions),
            "Reference != Final Base": int(ref_not_final_base),
            "Inconsistent Percentage": round(inconsistent_percentage, 4),
            "No Count Percentage": round(no_count_percentage, 4)
        }

        stats_by_type[gene_type] = stats

    # Write to output file
    with open(output_file, "w") as f:
        f.write("#Base-calling Statistics\n")
        f.write("---Overall---\n")
        for key, val in stats_all.items():
            f.write(f"{key}: {val}\n")
        f.write("\n")
        for gene_type, stats in stats_by_type.items():
            f.write(f"---{gene_type}---\n")
            for key, val in stats.items():
                f.write(f"{key}: {val}\n")
            f.write("\n")

def make_line_plot_cov(df, output_line_plot_path, gene_id=gene_id):  # Added gene_id parameter
    """
    Creates a line plot of coverage for a single person across genomic positions,
    highlighting positions with coverage < 5.

    Args:
        df (pd.DataFrame): DataFrame with 'Position' and 'Coverage' columns.
        output_line_plot_path (str): Path to save the plot.
        gene_id (str, optional): Gene ID for the plot title. Defaults to gene_id.
    """

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # Suppress all warnings in this block

        plt.figure(figsize=(8, 6))
        sns.lineplot(x="Position", y="Coverage", data=df, errorbar=None, alpha=0.7)

        # Highlight low coverage positions
        low_cov_data = df[df['Coverage'] < 2]
        if not low_cov_data.empty:
            plt.scatter(low_cov_data['Position'], low_cov_data['Coverage'],
                        color='red', marker='x', s=50, zorder=5, label="Coverage < 2")

        plt.grid(False)
        plt.title(f"Coverage for Gene {gene_id}", y=1.02)
        plt.xlabel("Position")
        plt.ylabel("Coverage")
        plt.legend(loc='upper right')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(output_line_plot_path)
        plt.close()
        


start_time = time.perf_counter()  # runtime measurement


locations, locations_with_type = transform_location_input(location_path)
#region_bed = f"{os.path.dirname(location_path)}/locations_transformed.bed"

type_df = assign_and_merge_types(locations_with_type)
type_df['Person'] = person
type_df = integrate_reference_bases(type_df, reference_genome=ref_path)

genomes = parse_genome(genome_path)
rep_df_list =[]
for genome in genomes:
    rep_df_list.append(create_mpileup(genome, locations))

rep_df_merged = merge_dataframes(rep_df_list)
df_merged = merge_dataframes([type_df,rep_df_merged])
result = process_replicates(df_merged)
result = identify_special_cases(result)
# Format all Frequency columns to string representation
freq_columns = result.filter(regex='^Frequency').columns
for col in freq_columns:
    result[col] = result[col].apply(format_frequencies)

ambiguous_positions = identify_ambigious_positions(result)
if not ambiguous_positions.empty:
    ambiguous_positions.to_csv(f"{output_path}/ambiguous_positions.tsv", sep="\t", index=False)
else:
    # If the ambiguous positions file exists, delete it
    ambiguous_file_path = f"{output_path}/ambiguous_positions.tsv"
    if os.path.exists(ambiguous_file_path):
        os.remove(ambiguous_file_path)


variants_positions = identify_variant_positions(result)
if not variants_positions.empty:
    variants_positions.to_csv(f"{output_path}/variants.tsv", sep="\t", index=False)
else:
    with open(f"{output_path}/variants.tsv", "w") as f:
        f.write("#No variant positions found.\n")

if not os.path.exists(output_path):
    try:
        os.makedirs(output_path)
        print(f"Created output directory at {output_path}")
    except OSError as e:
        print(f"Error creating output folder '{output_path}': {e}")

result.to_csv(f"{output_path}/results.tsv.gz", index=False, compression='gzip')
compute_stats(result, f"{output_path}/stats.tsv")
make_line_plot_cov(result, f"{output_path}/coverage_plot.png")

end_time = time.perf_counter()
runtime = end_time - start_time
print(f"Runtime {person}: {runtime:.5f} seconds")

