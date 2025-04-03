"""
Version 03.04.2025
"""
import subprocess
import os

def run_base_calling(gene_name, ID):
    for filename in os.listdir(input_genome_folder):
        person = os.path.splitext(os.path.basename(filename))[0]
        file_path = os.path.join(input_genome_folder, filename)
        output_folder = os.path.join(base_output_folder,gene_name,person)
        cmd = f"python3 /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/seq_extraction_3.py --genome {file_path} --reference /mnt/raidinput/input/own/ReferenceGenomes/human_g1k_v37.fasta.gz --location /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genes/{gene_name}/annot_{ID}_v2.tsv --output {output_folder}/ --person {person}"
        subprocess.run(cmd, shell=True, check=True)

def run_analysis(gene_name):
    cmd = f"python3 /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/analyze.py --gene {gene_name}"
    subprocess.run(cmd, shell=True, check=True)


def get_gene_list(gene_list_file):
    """
    :param gene_list_file: a tsv file containing IDs (1st column) and gene names(2nd column)
    :return: a dictionary mapping gene IDs to gene names
    """
    with open(gene_list_file, 'r') as tsv_file:
        gene_list = {}
        for line in tsv_file:
            line = line.strip()
            fields = line.split()
            gene_list[fields[0]] = fields[1]

    return gene_list


input_genome_folder = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/genomes"
base_output_folder = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family"
gene_list = get_gene_list("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/gene_list")


for id, gene_name in gene_list.items():
    run_base_calling(gene_name, id)
    run_analysis(gene_name)

