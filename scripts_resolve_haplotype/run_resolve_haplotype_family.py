import pandas as pd
import subprocess
import argparse
import os
import time

FAMILY = [
    "grandfather_father", "grandmother_father",  # paternal grandparents
    "father", "child", "mother", "aunt",         # nuclear family
    "grandmother_mother", "grandfather_mother"   # maternal grandparents
]

FAMILY_MALE = ["grandfather_father", "father", "child", "grandfather_mother"]
argparser = argparse.ArgumentParser(description="Resolve haplotypes from BAM files based on heterozygous positions for all samples")
argparser.add_argument("-c","--chunk", type=str, help="Chunk to process")


args = argparser.parse_args()
chunk = args.chunk

chrom = chunk.split("_")[0]
analysis_folder = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results/{chrom}/"
input_genome_bam_folder = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genomes/BAM"
gene_id_with_name_file = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genes/gene_id_with_name/gene_id_with_name_{chunk}.txt"

def get_gene_id_and_name_dict(gene_id_with_name_file):
    """
    :param gene_id_with_name_file: a tsv file containing IDs (1st column) and gene names(2nd column)
    :return: a dictionary mapping gene IDs to gene names
    """
    gene_dict = {}
    try:
        with open(gene_id_with_name_file, 'r') as tsv_file:
            for line in tsv_file:
                line = line.strip()
                if line and not line.startswith("#"):  # Skip empty lines and comments
                    fields = line.split()
                    if len(fields) == 2:
                        gene_dict[fields[0]] = fields[1]
                    else:
                        print(f"Warning: Skipping malformed line in gene list: {line}")
    except FileNotFoundError:
        print(f"Error: Gene list file not found: {gene_id_with_name_file}")
    return gene_dict

def run_resolve_haplotype(gene, person):
    input_bam = os.path.join(input_genome_bam_folder, f"{person}.txt")
    person_analysis_folder = os.path.join(analysis_folder, gene, person)
    
    cmd = ["/home/b/buit/miniconda3/envs/HiWi/bin/python", "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_resolve_haplotype/resolve_haplotype_individual.py",
            "-b", input_bam,
            "-v", os.path.join(person_analysis_folder, "variants.tsv")
            ]
    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    gene_dict = get_gene_id_and_name_dict(gene_id_with_name_file)
    start_time = time.time()
    for ID, _ in gene_dict.items():
        print(f"Resolve haplotypes for gene: {ID}")
        if chrom != "chrY":
            for person in FAMILY:
                run_resolve_haplotype(ID, person)
        else:
            for person in FAMILY_MALE:
                run_resolve_haplotype(ID, person)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Runtime for chunk {chunk}: {elapsed_time:.2f} seconds.")
