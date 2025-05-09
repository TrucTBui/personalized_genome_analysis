"""
Version 03.04.2025
"""
import subprocess
import os
import argparse
import time
parser = argparse.ArgumentParser()

parser.add_argument("-c", "--chunk", type=str, required=True, help="which chunk to analyse")

args = parser.parse_args()
chunk = args.chunk
chrom = chunk.split("_")[0]
MALES = ["child", "father", "grandfather_father", "grandfather_mother"]

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

def check_base_calling_outputs_exist(ID):
    """
    Check if the coverage_plot.png file exists for the given gene for every person.
    """
    FAMILY = ["father", "mother", "child", "grandmother_father", "grandmother_mother", "grandfather_father", "grandfather_mother", "aunt"]
    if not chrom == "chrY":
        for person in FAMILY:
            coverage_plot_file = os.path.join(base_output_folder, ID, person, "coverage_plot.png")
            if not os.path.exists(coverage_plot_file):
                return False
    else:
        for person in MALES:
            coverage_plot_file = os.path.join(base_output_folder, ID, person, "coverage_plot.png")
            if not os.path.exists(coverage_plot_file):
                return False
    
    return True

def check_analysis_outputs_exist(ID):
    """
    Checks if the main analysis output file exists for the given gene.
    """
    analysis_output_file = os.path.join(base_output_folder, ID, "results_all.tsv.gz")
    filtered_analysis_output_file = os.path.join(base_output_folder, ID, "merged_results_filtered.tsv")

    #check if the filtered_analysis_output_file really contains the filtered results or only blank lines
    """
    if os.path.exists(filtered_analysis_output_file):
        with open(filtered_analysis_output_file, 'r') as f:
            lines = f.readlines()
            if len(lines) < 2:
                return False
    """
    return os.path.exists(analysis_output_file) and os.path.exists(filtered_analysis_output_file)

def run_base_calling(ID, chunk, chrom):
    for filename in os.listdir(input_genome_folder):
        person = os.path.splitext(os.path.basename(filename))[0]
        file_path = os.path.join(input_genome_folder, filename)
        output_folder = os.path.join(base_output_folder,ID,person)
        os.makedirs(output_folder, exist_ok=True)
        cmd = f"/home/b/buit/miniconda3/envs/HiWi/bin/python  /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_base_calling_rep_assessment/seq_extraction_person.py --genome {file_path} --reference /mnt/raidinput/input/own/ReferenceGenomes/human_g1k_v37.fasta.gz --location /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genes/gtf_all_genes/{chunk}/{ID}.tsv --output {output_folder}/ --person {person}"
        if not chrom == "chrY":
            subprocess.run(cmd, shell=True, check=True)
        else:
            if person in MALES: # Only run for males in the Y chromosome
                subprocess.run(cmd, shell=True, check=True)
    print(f"Finished base calling for gene: {ID}")

def run_analysis(ID):
    cmd = f"/home/b/buit/miniconda3/envs/HiWi/bin/python  /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_base_calling_rep_assessment/analyze.py --gene {ID} --chunk {chunk}"
    subprocess.run(cmd, shell=True, check=True)
    print(f"Finished analysis for gene: {ID}")


input_genome_folder = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genomes"
base_output_folder = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results/{chrom}"

gene_id_with_name_file = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genes/gene_id_with_name/gene_id_with_name_{chunk}.txt"
#gtf_folder = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genes/gtf_all_genes/{chunk}"

gene_list = get_gene_id_and_name_dict(gene_id_with_name_file)

if __name__ == "__main__":
    start_time = time.perf_counter()
    os.makedirs(base_output_folder, exist_ok=True)

    if gene_list:
        for id, gene_name in gene_list.items():
            print(f"Processing gene: {id}")
            run_base_calling(id, chunk, chrom)
            run_analysis(id)
            
            if not check_base_calling_outputs_exist(id):
                run_base_calling(id, chunk, chrom)
            #else:
            #    print(f"Base calling outputs already exist for gene: {gene_name}. Skipping base calling.")
            """
            for person in os.listdir(os.path.join(base_output_folder, id)):
                person_folder = os.path.join(base_output_folder, id, person)
                unzipped_file = os.path.join(person_folder, "results.tsv")
                gzipped_file = unzipped_file + ".gz"
                if os.path.isfile(unzipped_file) and not os.path.exists(gzipped_file):
                    import gzip
                    import shutil
                    with open(unzipped_file, 'rb') as f_in, gzip.open(gzipped_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                    os.remove(unzipped_file)
            """
            if not check_analysis_outputs_exist(id):
                run_analysis(id)
            #else:
            #    print(f"Analysis outputs already exist for gene: {gene_name}. Skipping analysis.")
            
    else:
        print("No genes found in the gene list file.")

    end_time = time.perf_counter()
    runtime = end_time - start_time
    print(f"Runtime: {runtime:.5f} seconds")

