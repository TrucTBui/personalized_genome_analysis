"""
Version 31.01.2025
"""
import subprocess
import argparse
import os
import glob
import pandas as pd



input_genome_folder = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/genomes"
base_output_folder = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family"

def run(gene_name, ID):
    for filename in os.listdir(input_genome_folder):
        person = os.path.splitext(os.path.basename(filename))[0]
        file_path = os.path.join(input_genome_folder, filename)
        output_folder = os.path.join(base_output_folder,gene_name,person)
        #print(output_folder)
        cmd = f"python3 /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/seq_extraction.py --genome {file_path} --reference /mnt/raidinput/input/own/ReferenceGenomes/human_g1k_v37.fasta.gz --location /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/{gene_name}/annot_{ID}.tsv --output {output_folder}/ --person {person}"
        #print(cmd)
        subprocess.run(cmd, shell=True, check=True)

def find_diff_in_results(result_df, output):
    differences=[]
    final_base_cols = [col for col in result_df.columns if col not in ["Chromosome", "Position", "Gene", "Type", "Reference"]]

    for index, row in result_df.iterrows():
        ref = row.Reference  # Reference base

        alts = [row[col] for col in final_base_cols]  # Collect all bases from the family


        if len(set(alts)) > 1 or any('/' in alt for alt in alts) or alts[0] is not ref:  # Diff in read
           differences.append(row)


    diff_df = pd.DataFrame(differences)

    diff_df.to_csv(output, sep="\t", index=False)

gene_list = {"ENSG00000160791":"CCR5"}


for id, gene_name in gene_list.items():
    run(gene_name, id)


df_results_list = []
df_special_cases_list = []
df_seqs_list =[]

# Iterate over gene subfolders
for gene_name in gene_list.values():  # Based on your `gene_list`
    gene_folder = os.path.join(base_output_folder, gene_name)

    # Find all "*.results_merged.tsv" files inside the gene folder (recursively)
    result_files = glob.glob(os.path.join(gene_folder, "**", "results_merged.tsv"), recursive=True)
    special_cases_files = glob.glob(os.path.join(gene_folder, "**", "special_cases.tsv"), recursive=True)
    seq_files = glob.glob(os.path.join(gene_folder, "**", "sequence_merged.tsv"), recursive=True)

    for file in result_files:
        # Read the file into a DataFrame
        df = pd.read_csv(file, sep="\t")  # Ensure tab-separated format

        # Add a "Gene" column to keep track of the associated gene
        df["Gene"] = gene_name

        df_results_list.append(df)

    for file in special_cases_files:
        # Read the file into a DataFrame
        df = pd.read_csv(file, sep="\t")  # Ensure tab-separated format

        df_special_cases_list.append(df)

    for file in seq_files:
        # Read the file into a DataFrame
        df = pd.read_csv(file, sep="\t")  # Ensure tab-separated format

        df_seqs_list.append(df)

# Merge all DataFrames into one
merged_results_df = pd.concat(df_results_list, ignore_index=True)

merged_special_cases_df = pd.concat(df_special_cases_list, ignore_index=True)
merged_special_cases_df = merged_special_cases_df.sort_values(by="Position")
merged_special_cases_df.reset_index(drop=True, inplace=True)
merged_special_cases_df.to_csv(os.path.join(base_output_folder, "CCR5", "merged_special_cases.tsv"), sep="\t", index=False)

merged_seqs_df = pd.concat(df_seqs_list, ignore_index=True)
merged_seqs_df = merged_seqs_df.sort_values(by="Location")
merged_seqs_df.reset_index(drop=True, inplace=True)
merged_seqs_df.to_csv(os.path.join(base_output_folder, "CCR5", "merged_seqs.tsv"), sep="\t", index=False)

pivoted_df = merged_results_df.pivot_table(
    index=["Chromosome", "Position", "Gene", "Type","Reference"],
    columns="Person",
    values="Alternative",
    aggfunc="first"  # Choose the first occurrence if there are duplicates
)

# Reset index to make it a flat DataFrame
pivoted_df.reset_index(inplace=True)

output_file = os.path.join(base_output_folder, "CCR5", "merged_results.tsv")  # TODO A bit harcoding
pivoted_df.to_csv(output_file, sep="\t", index=False)

# Filter for special cases (where the bases are not similar across all samples)
output_filtered_file = os.path.join(base_output_folder, "CCR5", "merged_results_filtered.tsv")  # TODO A bit harcoding
find_diff_in_results(pivoted_df, output_filtered_file)

"""
sc_pivoted_df = merged_special_cases_df.pivot_table(
    index=["Position"],
    columns="Person",
    values="Bases",
    aggfunc="first"  # Choose the first occurrence if there are duplicates
)

# Reset index to make it a flat DataFrame
sc_pivoted_df.reset_index(inplace=True)

sc_output_file = os.path.join(base_output_folder, "CCR5", "merged_special_cases.tsv")  # TODO A bit harcoding
#sc_pivoted_df.to_csv(sc_output_file, sep="\t", index=False)
"""

