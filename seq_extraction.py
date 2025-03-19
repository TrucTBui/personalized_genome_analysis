"""
version 19.03.2025
"""
import subprocess
import argparse
import os
from collections import Counter
import time
import pandas as pd

def extract_reads_from_sam(sam_file):
    reads = []

    with open(sam_file, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            chromosome = columns[2]
            start_position = int(columns[3])  # Start position (1-based)
            sequence = columns[9]          # Sequence of the read
            cigar = columns[5]

            """
            matched = True
            
            for letter in ["I", "D", "N", "S", "H", "P", "X"]:
                if letter in cigar:
                    matched = False
            if matched:
                reads.append([chromosome, start_position, sequence])
            """
            reads.append([chromosome, start_position, sequence])
    return reads

def reads_processing(all_reads, location):
    chromosome, pos_range = location.split(":")
    start_range, end_range = map(int, pos_range.split("-"))

    variant_dict = {}

    for read in all_reads:
        read_chromosome = read[0]
        read_start = read[1]  # Start pos of read
        read_end = read_start + len(read[2]) - 1

        # check for suitable reads
        if chromosome != read_chromosome:
            continue
        if read_end < start_range or read_start > end_range:
            continue

        # determine the actual start and end positions of the trimmed read
        actual_start = max(read_start, start_range)
        actual_end = min(read_end, end_range)

        # trim the sequence
        trimmed_sequence = read[2][actual_start - read_start:actual_end - read_start + 1]

        # find variants
        for i, base in enumerate(trimmed_sequence):
            # pos = f"{chromosome}:{actual_start + i}"
            pos = actual_start + i
            if pos not in variant_dict:
                variant_dict[pos] = []
            variant_dict[pos].append(base)

    return variant_dict

def find_consensus_base (replicates_dict, reference, location):
    seq = ""
    list_variants_dict = []
    for rep, reads in replicates_dict.items():
        variant_dict = reads_processing(reads, location)
        list_variants_dict.append(variant_dict)

    final_bases_dict = {}
    reference_dict = {}
    coverage_dict = {}
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

        for variant_dict in list_variants_dict:
            bases = variant_dict[pos]
            coverage_dict[pos] += len(bases)
            frequency = Counter(bases)  # Count how frequent each variant is
            list_frequecy.append(frequency)
            unique_bases = set(bases)

            if len(unique_bases) > 1:  # More than one unique nucleotide means a variant/ sequencing error
                sorted_bases = sorted(frequency.items(), key=lambda item: item[1], reverse=True)
                ordered_bases = [base for base, count in sorted_bases]

                most_frequent_base = ordered_bases[0]  # Extract the most significant base
                second_frequent_base = ordered_bases[1]
                ratio = (frequency[most_frequent_base] / sum(frequency.values()))
                ratio2 = (frequency[second_frequent_base] / sum(frequency.values()))
                ratio12 = ((frequency[most_frequent_base] + frequency[second_frequent_base]) / sum(frequency.values()))

                # Process the reads results: The most frequent base is supposed to be the correct one
                if (ratio >= 0.7 or (sum(frequency.values()) < 15 and ratio >= 0.6) or (
                        len(unique_bases) == 3 and (ratio >= 0.6 or ratio / ratio2 >= 1.5)) or
                        (len(unique_bases) >= 4 and ratio / ratio2 >= 1.5)):
                    list_processed_bases.append(most_frequent_base)

                else:  # The base is not clear, can be haplotype
                    if ratio12 > 0.65:
                        list_processed_bases.append(f"{most_frequent_base}/{second_frequent_base}")
                    else:
                        list_processed_bases.append("/".join(ordered_bases))
            else:  # clear base
                list_processed_bases.append(bases[0])

        # Compare all replicates to find consensus bases at the position
        # Absolutely clear situation, no conflicts at all
        if all('/' not in base for base in list_processed_bases) and len(set(list_processed_bases)) == 1:
            final_bases_dict[pos] = list_processed_bases[0]
            seq += final_bases_dict[pos]
        else:  # Sum up the frequency and do the ratio process again
            combined_frequencies = Counter()
            # NOTE: Sum up all dictionaries (KEY STEP)
            for freq_dict in list_frequecy:
                combined_frequencies.update(freq_dict)

            sorted_bases = sorted(combined_frequencies.items(), key=lambda item: item[1], reverse=True)

            ordered_bases = [base for base, count in sorted_bases]

            most_frequent_base = ordered_bases[0]  # Extract the most significant base
            second_frequent_base = ordered_bases[1]
            ratio = (combined_frequencies[most_frequent_base] / sum(combined_frequencies.values()))
            ratio2 = (combined_frequencies[second_frequent_base] / sum(combined_frequencies.values()))
            ratio12 = ((combined_frequencies[most_frequent_base] + combined_frequencies[second_frequent_base]) / sum(combined_frequencies.values()))

            # Process the reads results: The most frequent base is supposed to be the correct one
            if (ratio >= 0.7 or (sum(combined_frequencies.values()) < 15 and ratio >= 0.6) or (
                    len(combined_frequencies) == 3 and (ratio >= 0.6 or combined_frequencies[most_frequent_base]  / combined_frequencies[second_frequent_base] >= 1.5)) or
                    (len(combined_frequencies) >= 4 and combined_frequencies[most_frequent_base]  / combined_frequencies[second_frequent_base] >= 1.5)):
                final_bases_dict[pos] = most_frequent_base
                seq += most_frequent_base

            else:  # The base is not clear, can be haplotype
                if ratio12 > 0.65:
                    final_bases_dict[pos] = (f"{most_frequent_base}/{second_frequent_base}")
                    seq += (f"[{most_frequent_base}/{second_frequent_base}]")
                    all_special_cases[pos_name] = sorted_bases
                else:
                    final_bases_dict[pos] = 'N'  # Unknown bases
                    seq += 'N'
                    all_special_cases[pos_name] = sorted_bases
    return final_bases_dict, reference_dict, seq, coverage_dict

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
                type = fields[2]
                if type == "exon":  # Ignore, since UTR and CDS are already in exons
                    continue
                chromosome = fields[6]
                if chromosome.startswith("chr"):
                    chromosome = chromosome[3:]
                start = int(fields[3])
                end = int(fields[4]) - 1  # End is exclusive
                locations.append(f"{chromosome}:{start}-{end}")
                locations_with_type[(chromosome, start, end)] = type
    with open(f"{os.path.dirname(path)}/locations_transformed.txt", 'w') as w:
         for location in locations:
             w.write(f"{location}\n")
    return locations, locations_with_type

def print_output(output_path, results, seqs, all_special_cases):
    os.makedirs(output_path, exist_ok=True)  # Create the directory if it doesn't exist

    with open(f"{output_path}/results_merged.tsv", "w") as o:
        o.write("\n".join(results))
    with open(f"{output_path}/sequence_merged.tsv", "w") as o:
        o.write("\n".join(seqs))
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

"""
# The input file is already a csv file, each line containing the position
if os.path.isfile(location_path):  # Input is a file
    with open(location_path, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                locations.append(line.strip())
else:
    raise Exception("Location file not found or not correct")

ref_genome = extract_ref_genome(ref_path, location_path)
"""

if os.path.isfile(location_path):  # Input is the annotation file from ISAR
    locations, locations_with_type = transform_location_input(location_path)
    transformed_path = f"{os.path.dirname(location_path)}/locations_transformed.txt"
    ref_genome = extract_ref_genome(ref_path, transformed_path)

else:
    raise Exception("Location file not found or not correct")


results = [f"Person\tChromosome\tPosition\tType\tCoverage\tReference\tAlternative"]
seqs = [f"Person\tLocation\tType\tSequence"]

all_reads_of_replicates = {}
all_special_cases = {}
for genome in genomes:

    genome_basename = os.path.splitext(os.path.basename(genome))[0]

    # Extract all reads for the specified locations
    combined_locations = ' '.join(locations)
    temp_sam_path = f"/tmp/{genome_basename}.sam"
    command1 = f"(cd /mnt/proj/software && samtools view {genome} {combined_locations} | sort | uniq > {temp_sam_path})"
    subprocess.run(command1, shell=True, check=True)

    all_reads = extract_reads_from_sam(temp_sam_path)
    all_reads_of_replicates[genome_basename] = all_reads
    os.remove(f"/tmp/{genome_basename}.sam")

for location in locations:
    chrom, pos_range = location.split(":")
    start_pos, end_pos = map(int, pos_range.split("-"))

    # Extract ref genome at the current range
    ref = ref_genome.get(location)

    # Process reads and find variants for the current location
    type = locations_with_type[(chrom, start_pos, end_pos)]

    final_bases_dict, reference_dict, seq, coverage_dict = find_consensus_base(all_reads_of_replicates, ref, location)

    for pos in sorted(final_bases_dict.keys()):
        r = reference_dict[pos]
        p = final_bases_dict[pos]  # Final base
        c = coverage_dict[pos]  # coverage/ sequencing depth at that pos
        results.append(f"{person}\t{chrom}\t{pos}\t{type}\t{c}\t{r}\t{p}")
    seqs.append(f"{person}\t{location}\t{type}\t{seq}")

end_time = time.perf_counter()
runtime = end_time - start_time
print(f"Runtime: {runtime:.5f} seconds")

print_output(output_path, results, seqs, all_special_cases)