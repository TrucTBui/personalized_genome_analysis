"""
Version: 14.04.2025
"""
import subprocess
from typing import Counter
import pandas as pd
from itertools import combinations
import re
import argparse
import os

def get_bam_paths_from_input(input_path):
    with open(input_path, 'r') as f:
        bam_paths = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    return bam_paths

def filter_heterozygous_positions(variants_df):
    hetero_df = variants_df[variants_df["Special_Case"] == "Heterozygous"]
    return hetero_df

def parse_cigar(cigar_string):
    
    if cigar_string == "*":  # Unmapped read
        return []
    return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHPX=])', cigar_string)]  # d = digit, M=match, I=insertion, D=deletion, N=skipped region, S=soft clip, H=hard clip, P=padding, X=mismatch, ==sequence match

def get_base_from_read_at_ref_pos_1based(read_seq, read_ref_start_1based, cigar_ops, target_ref_pos_1based):
    # Finds the base in the read sequence corresponding to a target reference position, considering CIGAR operations.
    
    current_ref_pos_1based = read_ref_start_1based
    current_read_idx_0based = 0  

    for length, op in cigar_ops:
        if op == 'M' or op == '=' or op == 'X':  # Match, Sequence Match, Sequence Mismatch
            # Reference segment for this op: [current_ref_pos_1based, current_ref_pos_1based + length - 1]
            if current_ref_pos_1based <= target_ref_pos_1based < current_ref_pos_1based + length:
                offset_in_op = target_ref_pos_1based - current_ref_pos_1based
                read_char_idx = current_read_idx_0based + offset_in_op
                if 0 <= read_char_idx < len(read_seq):
                    return read_seq[read_char_idx]
                else:
                    return None
            current_ref_pos_1based += length
            current_read_idx_0based += length
        elif op == 'I':  # Insertion in read (consumes read bases, not reference bases)
            current_read_idx_0based += length
        elif op == 'D':  # Deletion in read (consumes reference bases, not read bases)
            # Reference segment for this op: [current_ref_pos_1based, current_ref_pos_1based + length - 1]
            if current_ref_pos_1based <= target_ref_pos_1based < current_ref_pos_1based + length:
                return None  # Target reference position is deleted in this read
            current_ref_pos_1based += length
        elif op == 'S':  # Soft clip (consumes read bases, not reference bases; bases are in SEQ)
            current_read_idx_0based += length
        elif op == 'N':  # Skipped region from reference; consumes reference, not read
            if current_ref_pos_1based <= target_ref_pos_1based < current_ref_pos_1based + length:
                return None # Target reference position is in a skipped region
            current_ref_pos_1based += length
        elif op == 'H':  # Hard clip (consumes neither ref nor read bases; bases are NOT in SEQ)
            pass # No change 

        # Optimization: if we've passed the target position in the reference alignment
        if current_ref_pos_1based > target_ref_pos_1based and op not in ('I', 'S', 'H'):
            return None # Target position not covered

    return None # Target position not covered by CIGAR operations

def resolve_haplotype_sequences(haplotype_df, bam_paths, distance_threshold=100):
   
    if haplotype_df.empty:
        return pd.DataFrame(columns=["Gene", "Chromosome", "Positions", "Haplotype", "ReadCount"])

    haplotype_df = haplotype_df.sort_values("Position").reset_index(drop=True)
    positions = haplotype_df["Position"].tolist() 
    if not positions: 
        return pd.DataFrame(columns=["Gene", "Chromosome", "Positions", "Haplotype", "ReadCount"])
        
    chromosome = haplotype_df["Chromosome"].iloc[0]
    gene = haplotype_df["Gene"].iloc[0]
    person = haplotype_df["Person"].iloc[0] 
    if chromosome == "Y":
        return pd.DataFrame(columns=["Gene", "Chromosome", "Positions", "Haplotype", "ReadCount"])
    if chromosome == "X" and person in ["father", "grandfather_father", "grandfather_mother", "child"]:
        return pd.DataFrame(columns=["Gene", "Chromosome", "Positions", "Haplotype", "ReadCount"])

    # Get blocks of positions that are close enough
    haplotype_blocks = []
        
    current_block = [positions[0]]
    for i in range(1, len(positions)):
        if positions[i] - positions[i - 1] <= distance_threshold:
            current_block.append(positions[i])
        else:
            if len(current_block) > 1: # Only consider blocks with at least 2 positions
                haplotype_blocks.append(current_block)
            current_block = [positions[i]]
    if len(current_block) > 1: # Add the last block if it has more than one position
        haplotype_blocks.append(current_block)

    if not haplotype_blocks:
        return pd.DataFrame(columns=["Gene", "Chromosome", "Positions", "Haplotype", "ReadCount"])

    haplotype_list = []
    for block in haplotype_blocks: 
        start_pos, end_pos = block[0], block[-1]

        all_reads_output = []

        for bam_path in bam_paths:
        
            cmd_list = ["samtools", "view", bam_path, f"{chromosome}:{start_pos}-{end_pos}"]
            try:
                result = subprocess.run(cmd_list, capture_output=True, text=True, check=False)
                if result.returncode != 0:
                    print(f"Warning: samtools command failed for region {chromosome}:{start_pos}-{end_pos}: {result.stderr}")
                    continue
            except FileNotFoundError:
                raise FileNotFoundError(f"samtools not found")

            reads_output = result.stdout.strip().split('\n')
            all_reads_output.extend(reads_output)

        for read_line in all_reads_output:
            if not read_line:
                continue
            
            fields = read_line.split('\t')
            if len(fields) < 10: 
                continue

            read_align_start_ref_1based = int(fields[3])
            cigar_str = fields[5]
            seq = fields[9]

            if cigar_str == "*" or seq == "*" or seq is None: # Skip unmapped reads or reads with no sequence
                continue
            
            cigar_ops = parse_cigar(cigar_str)
            if not cigar_ops: 
                continue

            read_haplotype = {} # Stores {position: base} for the current read

            for block_pos_1based in block: # These are 1-based positions
                base_in_read = get_base_from_read_at_ref_pos_1based(
                    seq,
                    read_align_start_ref_1based,
                    cigar_ops,
                    block_pos_1based
                )

                if base_in_read is not None:
                    # Check if this base is one of the expected heterozygous alleles
                    valid_alleles_series = haplotype_df.loc[haplotype_df["Position"] == block_pos_1based, "Final_Base"]
                    if not valid_alleles_series.empty:
                        allowed_alleles_str = valid_alleles_series.iloc[0]
                        allowed_alleles = allowed_alleles_str.split("/")
                        if base_in_read in allowed_alleles:
                            read_haplotype[block_pos_1based] = base_in_read
                            
                    else:
                        raise ValueError(f"Position {block_pos_1based} not found in haplotype DataFrame")
            
            if len(read_haplotype) >= 2:  # Need at least 2 positions covered to be a valid haplotype
                covered_positions = sorted(list(read_haplotype.keys()))
                                
                hap_str = ''.join(read_haplotype[p] for p in covered_positions)
                haplotype_list.append((tuple(covered_positions), hap_str))

    if not haplotype_list or len(haplotype_list) < 1:
        return pd.DataFrame(columns=["Gene", "Chromosome", "Positions", "Haplotype", "ReadCount"])

    # Count occurrences of each haplotype pair in BAM reads
    #print(haplotype_list)
    pair_counts = Counter(haplotype_list)

    processed_haplotypes_by_pair = {} # Key: (pos1, pos2), Value: list of (hap_str, count)

    for (pos_tuple, hap_str), count in pair_counts.items():
        if pos_tuple not in processed_haplotypes_by_pair:
            processed_haplotypes_by_pair[pos_tuple] = []
        processed_haplotypes_by_pair[pos_tuple].append((hap_str, count))

    output_data = []
    for pos_tuple, hap_list_with_counts in processed_haplotypes_by_pair.items():
        # Sort the haplotypes for this position_pair by count in descending order
        sorted_haps_for_pair = sorted(hap_list_with_counts, key=lambda x: x[1], reverse=True)
        
        # Select the top 2 (or fewer if less than 2 exist), because there can only be 2 haplotypes per pair 
        #if len(pos_tuple) == 2:
        sorted_haps_for_pair = sorted_haps_for_pair[:2]  # Only keep the top 2 haplotypes for pairs
        
        for hap_str, count in sorted_haps_for_pair:
            output_data.append((
                pos_tuple,
                hap_str,
                count
            ))
    # Merging inter-reads haplotypes
    def merge_haps(hap1, hap2):
        pos_tuple1, seq1, _ = hap1
        pos_tuple2, seq2, _ = hap2

        # Find overlap
        overlap = set(pos_tuple1) & set(pos_tuple2)
        if not overlap:
            return None

        for o in sorted(overlap):
            i1 = pos_tuple1.index(o)
            i2 = pos_tuple2.index(o)
            if seq1[i1] != seq2[i2]:
                return None  # Conflict

        # Merge positions and sequence
        new_pos = []
        new_seq = []
        for p, b in zip(pos_tuple1, seq1):
            new_pos.append(p)
            new_seq.append(b)
        for p, b in zip(pos_tuple2, seq2):
            if p not in new_pos:  # Add the rest of the second haplotype seq
                new_pos.append(p)
                new_seq.append(b)

        # Sort positions and reorder sequence
        zipped = sorted(zip(new_pos, new_seq), key=lambda x: x[0])
        merged_pos, merged_seq = zip(*zipped)
        return (list(merged_pos), ''.join(merged_seq), -1)  # -1 as placeholder for count, since we merged
    
    # Iteratively merge overlaping haplotypes until no more merges are possible
    merged = output_data.copy()
    #print(merged)
    while True:
        next_merged = []
        for hap in merged:
            merged_any = False
            for i, existing in enumerate(next_merged):
                #print(existing)
                result = merge_haps(hap, existing)
                if result:
                    next_merged[i] = result
                    merged_any = True
                    output_data.append(result)
                    break
            if not merged_any:
                next_merged.append(hap)
        if len(next_merged) == len(merged):
            break  # no more merging possible
        merged = next_merged
   
    output_df = pd.DataFrame(
        [(gene, chromosome, ','.join(map(str, pos)), hap, count) for pos, hap,count in output_data],
        columns=["Gene", "Chromosome", "Positions", "Haplotype", "ReadCount"]
    )
    #output_df = output_df.drop_duplicates()
    output_df = output_df.sort_values(by=["Positions", "ReadCount"], ascending=[True, False]).reset_index(drop=True)
    output_df.drop_duplicates(subset=["Chromosome", "Positions", "Haplotype"], inplace=True, keep='first')

    # If there is a Haplotype with length greater or equal 8: print the row
    output_df_filter = output_df[output_df["Haplotype"].str.len() >= 8]
    if not output_df_filter.empty:
        first_row_series = output_df_filter.iloc[0]
        string_values = [str(value) for value in first_row_series]
        tab_separated_row_string = '\t'.join(string_values)
        
        print(f"{tab_separated_row_string}\t{person}")
     
    return pd.DataFrame(output_df)
    #return pd.DataFrame([(chromosome, ','.join(map(str, pos)), hap, count) for pos, hap,count in output_data], columns=["Chromosome", "Positions", "Haplotype", "ReadCount"]).sort_values(by=["Positions", "ReadCount"], ascending=[True, False]).reset_index(drop=True)

"""
def resolve_haplotype_sequences2(haplotype_df, bam_path, distance_threshold=100):
    
    if haplotype_df.empty:
        return pd.DataFrame()

    haplotype_df = haplotype_df.sort_values("Position").reset_index(drop=True)
    positions = haplotype_df["Position"].tolist()
    chromosome = haplotype_df["Chromosome"].iloc[0]

    # Get blocks of positions that are close enough to be considered part of the same haplotype
    haplotype_blocks = []
    current_block = [positions[0]]
    for i in range(1, len(positions)):
        if positions[i] - positions[i - 1] <= distance_threshold:
            current_block.append(positions[i])
        else:
            if len(current_block) > 1:
                haplotype_blocks.append(current_block)
            current_block = [positions[i]]
    if len(current_block) > 1:
        haplotype_blocks.append(current_block)

    # Collect pairwise haplotypes
    haplotype_pairs = []
    for block in haplotype_blocks:
        start, end = block[0], block[-1]
        cmd = f"samtools view {bam_path} {chromosome}:{start}-{end}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        reads = result.stdout.strip().split('\n')

        for read in reads:  
            if not read:
                continue
            fields = read.split('\t')
            pos = int(fields[3])
            seq = fields[9]
            cigar = fields[5]
            read_haplotype = {}


            for block_pos in block:
                offset = block_pos - pos
                # TODO: Handle deletion and insertion cases

                if 0 <= offset < len(seq):
                    valid_allels = haplotype_df.loc[haplotype_df["Position"] == block_pos, "Final_Base"]
                    if not valid_allels.empty:
                        allowed_allels = valid_allels.iloc[0].split("/")
                        if seq[offset] in allowed_allels:
                            read_haplotype[block_pos] = seq[offset]
                        else:
                            print(f"Read {read} at position {block_pos} has invalid allele {seq[offset]} not in {allowed_allels}")

            if len(read_haplotype) >= 2:  # At least 2 positions covered in the read
                covered_positions = list(read_haplotype.keys())
                for sub_block in combinations(covered_positions, 2):
                    hap_str = ''.join(read_haplotype[p] for p in sub_block)
                    haplotype_pairs.append((tuple(sub_block), hap_str))


    def merge_haps(hap1, hap2):
        pos1, seq1 = hap1
        pos2, seq2 = hap2
        #print(pos2)

        # Find overlap
        overlap = set(pos1) & set(pos2)
        if not overlap:
            return None

        for o in sorted(overlap):
            i1 = pos1.index(o)
            i2 = pos2.index(o)
            if seq1[i1] != seq2[i2]:
                return None  # Conflict

        # Merge positions and sequence
        new_pos = []
        new_seq = []
        for p, b in zip(pos1, seq1):
            new_pos.append(p)
            new_seq.append(b)
        for p, b in zip(pos2, seq2):
            if p not in new_pos:  # Add the rest of the second haplotype seq
                new_pos.append(p)
                new_seq.append(b)

        # Sort positions and reorder sequence
        zipped = sorted(zip(new_pos, new_seq), key=lambda x: x[0])
        merged_pos, merged_seq = zip(*zipped)
        return (list(merged_pos), ''.join(merged_seq))
    
    # Iteratively merge overlaping haplotypes until no more merges are possible
    merged = haplotype_pairs.copy()
    while True:
        next_merged = []
        for hap in merged:
            merged_any = False
            for i, existing in enumerate(next_merged):
                result = merge_haps(hap, existing)
                if result:
                    next_merged[i] = result
                    merged_any = True
                    break
            if not merged_any:
                next_merged.append(hap)
        if len(next_merged) == len(merged):
            break  # no more merging possible
        merged = next_merged
   
    output_df = pd.DataFrame(
        [(chromosome, ','.join(map(str, pos)), hap) for pos, hap in merged],
        columns=["Chromosome", "Positions", "Haplotype"]
    )
    #output_df = output_df.drop_duplicates()
    output_df = output_df.sort_values(by="Positions").reset_index(drop=True)
    print(output_df)
    return output_df
"""

if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description="Resolve haplotypes from BAM files based on heterozygous positions.")
    parser.add_argument("-v", "--variant", type=str, help="Path to the input file containing the variant table.")
    parser.add_argument("-b", "--bam", type=str, help="Paths to the txt file containing BAM file paths, one per line.")
    args = parser.parse_args()
    variant_path = args.variant
    bams = args.bam
    
    with open(variant_path, "r") as f:
        lines = f.readlines()
    
    if lines[0].startswith("#No variant positions found"):
        hetero_df = pd.DataFrame(columns=["Gene", "Chromosome", "Position", "Type", "Person", "Reference_Base", "Final_Base", "Special_Case"])
    else:
        variants_df = pd.read_csv(variant_path, sep="\t")
        hetero_df = filter_heterozygous_positions(variants_df)

    bam_paths = get_bam_paths_from_input(bams)

    output_path = os.path.join(os.path.dirname(variant_path), "resolved_haplotypes.tsv")
    resolved_haplotypes_df = resolve_haplotype_sequences(hetero_df, bam_paths)
    resolved_haplotypes_df.to_csv(output_path, sep="\t", index=False)

