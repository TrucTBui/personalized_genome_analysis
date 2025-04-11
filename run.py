"""
Version 03.04.2025
"""
import subprocess
import os
import argparse
import time
parser = argparse.ArgumentParser()

parser.add_argument("-c", "--chromosome", type=str, required=True, help="which chromosome to analyse")

args = parser.parse_args()
chromosome = args.chromosome

def check_base_calling_outputs_exist(gene_name):
    """
    Checks if any output folders exist for the given gene from the base calling step.
    """
    gene_output_folder = os.path.join(base_output_folder, gene_name)
    if os.path.exists(gene_output_folder) and os.listdir(gene_output_folder):
        return True
    return False

def check_analysis_outputs_exist(gene_name):
    """
    Checks if the main analysis output file exists for the given gene.
    """
    analysis_output_file = os.path.join(base_output_folder, gene_name, "results_all.tsv.gz")
    return os.path.exists(analysis_output_file)

def run_base_calling(gene_name, ID, chromosome):
    for filename in os.listdir(input_genome_folder):
        person = os.path.splitext(os.path.basename(filename))[0]
        file_path = os.path.join(input_genome_folder, filename)
        output_folder = os.path.join(base_output_folder,gene_name,person)
        cmd = f"/home/b/buit/miniconda3/envs/HiWi/bin/python  /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/seq_extraction_3.py --genome {file_path} --reference /mnt/raidinput/input/own/ReferenceGenomes/human_g1k_v37.fasta.gz --location /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genes/{chromosome}/{ID}.tsv --output {output_folder}/ --person {person}"
        subprocess.run(cmd, shell=True, check=True)
    print(f"Finished base calling for gene: {gene_name}")

def run_analysis(gene_name):
    cmd = f"/home/b/buit/miniconda3/envs/HiWi/bin/python  /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/analyze.py --gene {gene_name} --chromosome {chromosome}"
    subprocess.run(cmd, shell=True, check=True)
    print(f"Finished analysis for gene: {gene_name}")


def get_gene_list(gene_list_file):
    """
    :param gene_list_file: a tsv file containing IDs (1st column) and gene names(2nd column)
    :return: a dictionary mapping gene IDs to gene names
    """
    gene_list = {}
    try:
        with open(gene_list_file, 'r') as tsv_file:
            for line in tsv_file:
                line = line.strip()
                if line and not line.startswith("#"):  # Skip empty lines and comments
                    fields = line.split()
                    if len(fields) == 2:
                        gene_list[fields[0]] = fields[1]
                    else:
                        print(f"Warning: Skipping malformed line in gene list: {line}")
    except FileNotFoundError:
        print(f"Error: Gene list file not found: {gene_list_file}")
    return gene_list


input_genome_folder = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/genomes"
base_output_folder = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/{chromosome}"
gene_list = get_gene_list(f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genes/chunks/gene_list_{chromosome}.txt")

"""
chromosome_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                   'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                   'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']
"""

if __name__ == "__main__":
    start_time = time.perf_counter()
    if not os.path.exists(base_output_folder):
        os.mkdir(base_output_folder)
        print(f"Created base output folder: {base_output_folder}")

    if gene_list:
        for id, gene_name in gene_list.items():
            run_base_calling(gene_name, id, chromosome)
            if not check_base_calling_outputs_exist(gene_name):
                run_base_calling(gene_name, id, chromosome)
            else:
                print(f"Base calling outputs already exist for gene: {gene_name}. Skipping base calling.")

            if not check_analysis_outputs_exist(gene_name):
                run_analysis(gene_name)
            else:
                print(f"Analysis outputs already exist for gene: {gene_name}. Skipping analysis.")
    else:
        print("No genes found in the gene list file.")

    end_time = time.perf_counter()
    runtime = end_time - start_time
    print(f"Runtime: {runtime:.5f} seconds")

