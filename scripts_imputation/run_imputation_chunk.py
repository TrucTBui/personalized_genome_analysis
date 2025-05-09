import os
import argparse
import subprocess
import time

parser = argparse.ArgumentParser(description="Run imputation for a specific chunk.")
parser.add_argument("-c", "--chunk", type=str, required=True, help="which chunk to analyse")

args = parser.parse_args()
chunk = args.chunk
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

if __name__ == "__main__":
    start_time = time.perf_counter()

    gene_dict = get_gene_id_and_name_dict(gene_id_with_name_file)
    for ID, gene in gene_dict.items():
        print(f"Imputation for gene: {ID}")
        cmd = f"/home/b/buit/miniconda3/envs/HiWi/bin/python  /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_imputation/impute_family.py -c {chunk} -g {ID}"
        subprocess.run(cmd, shell=True, check=True)

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Total time taken: {elapsed_time:.2f} seconds")

