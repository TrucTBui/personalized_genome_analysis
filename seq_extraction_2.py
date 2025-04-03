"""
version 25.03.2025
Diff from the 1st version (seq_extraction.py):
- the genomic base counting process is not done naively, but with samtools
- use new annotation version from ISAR
"""
import subprocess
import argparse
import os
from collections import Counter, defaultdict
import time
from intervaltree import IntervalTree
import collections
import re
import pandas as pd


def variant_dict_from_samtools(sam_file, region=None):
    """
    Generates a variant dictionary from a SAM file using samtools mpileup,
    handling insertions, deletions, and sequencing errors.

    Args:
        sam_file (str): Path to the SAM file.
        region (str, optional): Genomic region in the format "chr:start-end".
                                If None, processes the entire SAM file.

    Returns:
        dict: A dictionary where keys are genomic positions and values are lists of bases.
    """
    variant_dict = collections.defaultdict(list)
    samtools_command = ["samtools", "mpileup", "-aa", sam_file]

    if region:
        samtools_command.extend(["-r", region])

    try:
        mpileup_output = subprocess.check_output(samtools_command, text=True)
        lines = mpileup_output.strip().split('\n')

        for line in lines:
            parts = line.split('\t')
            if len(parts) < 5:
                continue  # Skip lines that don't have the required number of fields

            chromosome = int(parts[0])
            position = int(parts[1])  # 1-based genomic position
            bases = parts[4]

            parsed_bases = parse_mpileup_bases(bases)
            variant_dict[position] = parsed_bases

    except FileNotFoundError:
        print(f"Error: SAM file '{sam_file}' not found.")
        return {}
    except subprocess.CalledProcessError as e:
        print(f"Error: samtools mpileup failed: {e}")
        return {}

    return variant_dict

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

def find_consensus_base_new (sams, reference, location):
    list_variants_dict = []
    list_id = []
    for id,sam in sams.items():
        variant_dict = variant_dict_from_samtools(sam, location)
        list_variants_dict.append(variant_dict)
        list_id.append(id)

    seq = ""
    final_bases_dict = {}
    reference_dict = {}
    coverage_dict = {}
    frequency_dict ={}

    chromosome, _ = location.split(":")
    counter = 0  # for extracting bases at specific positions of the reference
    # go through every position of the variant_dict
    # all replicate have the same variant dict format so just take the first one of the list
    for pos in sorted(list_variants_dict[0].keys()):
        pos_name = f"{chromosome}:{pos}"
        r = reference[counter]
        reference_dict[pos] = r
        counter+=1
        list_frequecy = []
        list_processed_bases = []
        coverage_dict[pos] = 0
        frequency_dict[pos] = {} # Initialize the inner dictionary for the current position

        for i in range(len(list_variants_dict)):
            variant_dict = list_variants_dict[i]
            id = list_id[i]
            bases = variant_dict[pos]
            coverage_dict[pos] += len(bases)
            frequency = Counter(bases)  # Count how frequent each variant is
            list_frequecy.append(frequency)

            # Add the base frequency of each rep to the dict
            frequency_dict[pos][id] = frequency

            unique_bases = set(bases)

            if len(unique_bases) == 0:  # Empty counts
                list_processed_bases.append("N")

            elif len(unique_bases) > 1:  # More than one unique nucleotide means a variant/ sequencing error
                sorted_bases = sorted(frequency.items(), key=lambda item: item[1], reverse=True)
                ordered_bases = [base for base, count in sorted_bases]

                most_frequent_base = ordered_bases[0]  # Extract the most significant base
                second_frequent_base = ordered_bases[1]
                ratio = (frequency[most_frequent_base] / sum(frequency.values()))
                ratio12 = ((frequency[most_frequent_base] + frequency[second_frequent_base]) / sum(frequency.values()))

                # Process the reads results: The most frequent base is supposed to be the correct one
                if (ratio >= 0.75 or (sum(frequency.values()) < 15 and ratio >= 0.6) or (
                        len(unique_bases) == 3 and (ratio >= 0.65 or frequency[most_frequent_base] / frequency[second_frequent_base] >= 2)) or
                        (len(unique_bases) >= 4 and frequency[most_frequent_base] / frequency[second_frequent_base] >= 2)):
                    list_processed_bases.append(most_frequent_base)

                else:  # The base is not clear, can be haplotype
                    if ratio12 > 0.65:
                        list_processed_bases.append(f"{most_frequent_base}/{second_frequent_base}")
                    else:
                        list_processed_bases.append("/".join(ordered_bases))
            else:  # clear base
                list_processed_bases.append(bases[0])


        # Compare all replicates to find consensus bases at the position
        if "N" in list_processed_bases:
            final_bases_dict[pos] = 'N'  # Unknown bases
            seq += 'N'
            all_special_cases[pos_name] = {}  # Empty dict

        # Absolutely clear situation, no conflicts at all
        elif all('/' not in base for base in list_processed_bases) and len(set(list_processed_bases)) == 1:
            final_bases_dict[pos] = list_processed_bases[0]
            seq += final_bases_dict[pos]
            if list_processed_bases[0] != r:
                combined_frequencies = Counter()
                # NOTE: Sum up all dictionaries (KEY STEP)
                for freq_dict in list_frequecy:
                    combined_frequencies.update(freq_dict)
                all_special_cases[pos_name] = dict(combined_frequencies)

        else:  # Sum up the frequency and do the ratio process again
            combined_frequencies = Counter()
            # NOTE: Sum up all dictionaries (KEY STEP)
            for freq_dict in list_frequecy:
                combined_frequencies.update(freq_dict)

            # No counts at the position (error handling)
            if len(combined_frequencies) == 0:
                final_bases_dict[pos] = 'N'  # Unknown bases
                seq += 'N'
                all_special_cases[pos_name] = {}  # Empty dict

            else:
                sorted_bases = sorted(combined_frequencies.items(), key=lambda item: item[1], reverse=True)

                ordered_bases = [base for base, count in sorted_bases]

                most_frequent_base = ordered_bases[0]  # Extract the most significant base
                second_frequent_base = ordered_bases[1]
                ratio = (combined_frequencies[most_frequent_base] / sum(combined_frequencies.values()))
                ratio12 = ((combined_frequencies[most_frequent_base] + combined_frequencies[second_frequent_base]) / sum(combined_frequencies.values()))

                # Process the reads results: The most frequent base is supposed to be the correct one
                if (ratio >= 0.75 or (sum(combined_frequencies.values()) < 15 and ratio >= 0.6) or (
                        len(combined_frequencies) == 3 and (ratio >= 0.65 or combined_frequencies[most_frequent_base]  / combined_frequencies[second_frequent_base] >= 2)) or
                        (len(combined_frequencies) >= 4 and combined_frequencies[most_frequent_base]  / combined_frequencies[second_frequent_base] >= 2)):
                    final_bases_dict[pos] = most_frequent_base
                    seq += most_frequent_base

                    if most_frequent_base is not r:
                        all_special_cases[pos_name] = sorted_bases

                else:  # The base is not clear, can be haplotype
                    if ratio12 > 0.65:
                        final_bases_dict[pos] = (f"{most_frequent_base}/{second_frequent_base}")
                        seq += (f"[{most_frequent_base}/{second_frequent_base}]")
                        all_special_cases[pos_name] = sorted_bases
                    else:
                        final_bases_dict[pos] = 'N'  # Unknown bases
                        seq += 'N'
                        all_special_cases[pos_name] = sorted_bases
    return final_bases_dict, reference_dict, seq, coverage_dict, frequency_dict

def extract_ref_genome(fa, locations):
    """
    :param fa: path to fasta file containing reference genome
    :param locations: a file containing loctions to extract
    :return: a dictionary, where key=location, value=sequence at the location
    """
    cmd = f"samtools faidx {fa} -r {locations}"
    result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)
    output = result.stdout.decode("utf-8")  # Decode the `stdout` attribute

    seqs= {}
    seq = ""
    location = None

    for line in output.splitlines():
        if line.startswith(">"):  # Header line
            if location:
                seqs[location] = seq
            location = line[1:]  # Extract location name (without '>')
            seq = ""
        else:
            seq += line.strip()

    # Add the last sequence to the dictionary
    if location:
        seqs[location] = seq

    return seqs


def parse_genome(path):
    genomes = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            genomes.append(line.strip())

    return genomes

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

                location_string = f"{chromosome}:{start}-{end}"
                locations.append(location_string)

                if location_string not in locations_with_type:
                    locations_with_type[location_string] = set()
                locations_with_type[location_string].add(type)

    with open(f"{os.path.dirname(path)}/locations_transformed_v2.txt", 'w') as w:
         for location in locations:
             w.write(f"{location}\n")
    return locations, locations_with_type

def merge_regions(regions):
    chrom_trees = defaultdict(IntervalTree)
    for region in regions:
        chrom, coords = region.split(":")
        start, end = map(int, coords.split("-"))
        chrom_trees[chrom][start:end] = None  # No need to store extra data

    # Merge intervals for each chromosome
    merged_regions = []
    for chrom, tree in chrom_trees.items():
        tree.merge_overlaps()
        for iv in sorted(tree):
            merged_regions.append(f"{chrom}:{iv.begin}-{iv.end}")

    #print(merged_regions)
    return merged_regions

def count_reads(bam, region=None):
    """Counts reads in a BAM file, optionally within a specified region."""
    if region:
        all_count_cmd = f"samtools view -c {bam} {region}"
        filtered_out_count_cmd = f"samtools view {bam} {region} | awk -F'\t' '$6 ~ /[IDNPX*]/' | wc -l"
    else:
        all_count_cmd = f"samtools idxstats {bam} | awk 'BEGIN {{sum=0}} {{sum+=$3+$4}} END {{print sum}}'"
        filtered_out_count_cmd = f"samtools idxstats {bam} | awk 'BEGIN {{sum=0}} {{sum+=$4}} END {{print sum}}'"

    try:
        all_result = subprocess.run(all_count_cmd, shell=True, capture_output=True, text=True, check=True)
        filtered_out_result = subprocess.run(filtered_out_count_cmd, shell=True, capture_output=True, text=True, check=True)
        return int(all_result.stdout.strip()), int(filtered_out_result.stdout.strip())
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        return None

def print_output(output_path, results, seqs, all_special_cases):
    os.makedirs(output_path, exist_ok=True)  # Create the directory if it doesn't exist

    with open(f"{output_path}/results_merged.tsv", "w") as o:
        o.write("\n".join(results))
    with open(f"{output_path}/sequence_merged.tsv", "w") as o:
        o.write("\n".join(seqs))
    freq_file_path = f"{output_path}/frequency.tsv"
    with open(freq_file_path, "w") as o:
        o.write("\n".join(freqs))
    process_replicate(freq_file_path)
    with open(f"{output_path}/special_cases.tsv", "w") as o:
        o.write("Person\tPosition\tA\tT\tG\tC\n")
        for pos in sorted(all_special_cases.keys()):
            case = all_special_cases[pos]
            case = dict(case)
            A = case["A"] if "A" in case else 0
            T = case["T"] if "T" in case else 0
            G = case["G"] if "G" in case else 0
            C = case["C"] if "C" in case else 0
            o.write(f"{person}\t{pos}\t{A}\t{T}\t{G}\t{C}\n")

def print_read_stats(output_path, reads_count):
    with open(f"{output_path}/reads_stats.tsv", "w") as o:
        o.write("".join(reads_count))

def process_replicate(freq_file):
    df = pd.read_csv(freq_file, sep="\t")
    df = df.drop_duplicates()
    consistency = []
    freq_columns = df.filter(regex='^Frequency').columns

    for _, row in df.iterrows():
        freq_values = row[freq_columns].values

        if any(pd.isna(freq) or freq == "NA" for freq in freq_values):
            consistency.append("no_counts")
            continue  # Skip further processing for this row, since no counts

        freq_counts_per_replicate = []

        # Parse frequency values for each replicate
        for freq in freq_values:
            replicate_counts = {}
            alleles = freq.split(";")
            for allele in alleles:
                if allele:
                    base, count = allele.split(":")
                    count = int(count)
                    replicate_counts[base] = count
            freq_counts_per_replicate.append(replicate_counts)

        # Determine major alleles in each replicate
        major_alleles_per_replicate = []
        for counts in freq_counts_per_replicate:
            total = sum(counts.values())
            sorted_alleles = sorted(counts.items(), key=lambda x: x[1], reverse=True)

            # Identify major alleles (top 1 or 2) covering at least 90% of total counts
            major_alleles = []
            cumulative_count = 0
            for base, count in sorted_alleles:
                major_alleles.append(base)
                cumulative_count += count
                if cumulative_count / total >= 0.80:
                    break
            major_alleles_per_replicate.append(set(major_alleles))

        # Check if major alleles are consistent across all replicates
        first_replicate_major = major_alleles_per_replicate[0]
        if all(major == first_replicate_major for major in major_alleles_per_replicate):
            consistency.append("normal")
        else:
            consistency.append("inconsistent")

    df["Consistency"] = consistency
    df.to_csv(freq_file, sep="\t", index=False)

    #TODO: Stats (bash). Make this faster

    return df

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


start_time = time.perf_counter()  # runtime measurement


# Parse the input path to a list of bam files (index 0) and vcf files(index1)
genomes= parse_genome(genome_path)

if os.path.isfile(location_path):  # Input is the annotation file from ISAR
    locations, locations_with_type = transform_location_input(location_path)
    non_overlapping_regions = merge_regions(locations)
    transformed_path = f"{os.path.dirname(location_path)}/locations_transformed_v2.txt"
    ref_genome = extract_ref_genome(ref_path, transformed_path)

else:
    raise Exception("Location file not found or not correct")


results = [f"Person\tChromosome\tPosition\tType\tCoverage\tReference\tAlternative"]
seqs = [f"Person\tLocation\tSequence"]
reads_count = [f"Person\t#All_Reads\t#All_Removed_Reads\tRemoval_Percentage\tWhole_Genome\n"]
freqs=[]
freqs_str = f"Person\tChromosome\tPosition\tType\tReference\t"

all_reads_of_replicates = {}
all_special_cases = {}
all_reads_count = 0
all_removed_reads_count = 0
all_reads_count_genome = 0
all_removed_reads_count_genome = 0

sam_dict = {}
fred_dict_all ={}
genome_basename_list = []
combined_locations = ' '.join(locations)
combined_non_overlapping_regions = " ".join(non_overlapping_regions)

for genome in genomes:

    genome_basename = os.path.splitext(os.path.basename(genome))[0]
    genome_basename_list.append(genome_basename)

    freqs_str += f"Frequency_{genome_basename}\t"

    # Extract all reads for the specified locations
    temp_sam_path = f"/home/b/buit/tmp/{genome_basename}.bam"
    command1 = f"(cd /mnt/proj/software && samtools view -b -h {genome} {combined_non_overlapping_regions} > {temp_sam_path})"
    subprocess.run(command1, shell=True, check=True)
    command_index = f"samtools sort -@ 8 {temp_sam_path} -o {temp_sam_path} && samtools index {temp_sam_path}"
    subprocess.run(command_index, shell=True, check=True)

    sam_dict[genome_basename] = temp_sam_path

    # count reads:
    read_count, filtered_out_count = count_reads(genome, combined_non_overlapping_regions)
    all_reads_count += read_count
    all_removed_reads_count += filtered_out_count
    read_count_genome, filtered_out_count_genome = count_reads(genome)
    all_reads_count_genome += read_count_genome
    all_removed_reads_count_genome += filtered_out_count_genome

    # For each person
    reads_count.append(f"{person}:{genome_basename}\t{read_count}\t{filtered_out_count}\t{round(all_removed_reads_count/all_reads_count, 4)}\t-\n")
    reads_count.append(f"{person}:{genome_basename}\t{read_count_genome}\t{filtered_out_count_genome}\t{round(filtered_out_count_genome/read_count_genome, 4)}\t+\n")

freqs.append(freqs_str)

for location in locations:
    chrom, pos_range = location.split(":")
    start_pos, end_pos = map(int, pos_range.split("-"))

    # Extract ref genome at the current range
    ref = ref_genome.get(location)

    # Process reads and find variants for the current location
    type = locations_with_type[f"{chrom}:{start_pos}-{end_pos}"]
    #type_string_region = ";".join(sorted(type))

    final_bases_dict, reference_dict, seq, coverage_dict, fred_dict = find_consensus_base_new(sam_dict, ref, location)

    for pos in sorted(final_bases_dict.keys()):

        types_for_pos = []
        for location_str, types in locations_with_type.items():
            chr_loc, pos_range_loc = location_str.split(":")
            start_loc, end_loc = map(int, pos_range_loc.split("-"))

            if chr_loc == chrom and start_loc <= pos <= end_loc:
                types_for_pos.extend(types)

        type_string = "/".join(sorted(set(types_for_pos)))  # Merge and sort types

        r = reference_dict[pos]
        p = final_bases_dict[pos]  # Final base
        c = coverage_dict[pos]  # coverage/ sequencing depth at that pos
        f = fred_dict[pos]
        results.append(f"{person}\t{chrom}\t{pos}\t{type_string}\t{c}\t{r}\t{p}")
        freq_string = f"{person}\t{chrom}\t{pos}\t{type_string}\t{r}\t"
        for id in genome_basename_list:
            counter = f[id]
            if len(counter) >0:
                for idx, (key, value) in enumerate(counter.items()):
                    if idx == len(counter) - 1:
                        freq_string += f"{key}:{value}"
                    else:
                        freq_string += f"{key}:{value};"
            else:
                freq_string+="NA"
            freq_string+="\t"
        freqs.append(freq_string)
    seqs.append(f"{person}\t{location}\t{seq}")

print_output(output_path, results, seqs, all_special_cases)
print_read_stats(output_path, reads_count)

end_time = time.perf_counter()
runtime = end_time - start_time
print(f"Runtime: {runtime:.5f} seconds")