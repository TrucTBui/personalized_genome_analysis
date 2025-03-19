"""
Script used to add clinvar annotation and integrate into analysis data.
The annotation is downloaded on ncbi clinvar database by searching for a specific gene.
version: 12.03.2025
"""

import pandas as pd
import os
import re
import gzip
import subprocess


def clinvar_to_df(clinvar_folder, output_tsv=None, skip_first_line=True):
    """
    Extracts information from all .tsv files in a folder and returns a concatenated pandas DataFrame.

    Args:
        clinvar_folder (str): The path to the folder containing .tsv files.
        output_tsv (str, optional): The path to save the combined DataFrame as a TSV file. Defaults to None.
        skip_first_line (bool): Whether to skip the first line of each file. Defaults to True.

    Returns:
        pd.DataFrame: A DataFrame containing the extracted information from all files.
    """
    all_data = []
    for filename in os.listdir(clinvar_folder):
        if filename.endswith(".tsv") and filename.startswith("clinvar"):
            file_path = os.path.join(clinvar_folder, filename)
            data = []
            first_line_skipped = False
            with open(file_path, "r") as f:
                for line in f:
                    if skip_first_line and not first_line_skipped:
                        first_line_skipped = True
                        continue

                    fields = line.strip().split("\t") #strip to remove trailing newlines, avoid empty fields.
                    if len(fields) >= 16:
                        name = fields[0]
                        condition = fields[3]
                        chrom = fields[5]
                        pos = fields[6]
                        dbSNP_ID = fields[11]
                        variant_type = fields[13]
                        molecular_consequence = fields[14]
                        classification = fields[15]
                        if variant_type == "single nucleotide variant":
                            data.append({
                                "name": name,
                                "condition": condition,
                                "chromosome": chrom,
                                "position": pos,
                                "dbSNP_ID": dbSNP_ID,
                                "molecular_consequence": molecular_consequence,
                                "classification": classification,
                            })
                    else:
                        print(f"Warning: Skipping line with insufficient fields in {filename}: {line.strip()}")

            all_data.extend(data) #add all data from current file to the all_data list.

    df = pd.DataFrame(all_data)
    if output_tsv:
        df.to_csv(output_tsv, sep='\t', index=False)
    return df

def filter_vcf_by_genes_grep(clinvar_vcf_file, output_vcf_file, gene_list_file):
    """
           Filters a VCF file using bcftools and a list of genes from a file.

           Args:
               clinvar_vcf_file (str): Path to the input VCF file (.vcf.gz).
               output_vcf_file (str): Path to the output VCF file (.vcf). That output will then
               bill compressed into a .vcf.gz file.
               gene_list_file (str): Path to the file containing gene IDs (one per line).
    """

    try:
        # Read gene names from file
        with open(gene_list_file, "r") as f:
            gene_names = [line.strip() for line in f if line.strip()]

        if not gene_names:
            raise ValueError("Gene list file is empty.")

        # Run grep
        gene_pattern = "|".join(gene_names)
        grep_command = f" zgrep -E 'GENEINFO=([^:|]+):({gene_pattern})' {clinvar_vcf_file} > filtered_variants.txt"

        print(f"Running command: {grep_command}")
        subprocess.run(grep_command, shell=True, check=True)

        # Extract header & filtered variants, then compress
        with open(output_vcf_file, "wb") as out_f:
            subprocess.run(f"zgrep '^#' {clinvar_vcf_file}", shell=True, stdout=out_f, check=True)
            subprocess.run(f"cat filtered_variants.txt", shell=True, stdout=out_f, check=True)
            subprocess.run(f"bgzip -f {output_vcf_file}", shell=True, check=True)

        os.remove("filtered_variants.txt")

    except subprocess.CalledProcessError as e:
        print(f"Error filtering VCF file: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def transform_clinvar_vcf(clinvar_vcf_file, output_tsv):
    data = []
    if clinvar_vcf_file.endswith(".gz"):
        open_func = gzip.open
    else:
        open_func = open

    with open_func(clinvar_vcf_file, "rt") as infile, open(output_tsv, "w") as outfile:
        # Write TSV header
        header = [
            "name", "condition", "chromosome",
            "position", "dbSNP_ID", "molecular_consequence", "classification"
        ]

        outfile.write("\t".join(header) + "\n")
        # Process each line
        for line in infile:
            if line.startswith("#"):
                continue
            fields = line.strip().split()
            if len(fields) < 8:
                continue  # Skip malformed lines

            chrom = fields[0]  # Chromosome
            pos = fields[1]  # Position
            #ref = fields[3]
            #alt = fields[4]
            info = fields[7]  # INFO field

            # Extract relevant fields from INFO
            name_match = re.search(r"CLNHGVS=([^;]+)", info)
            name = name_match.group(1) if name_match else "NA"

            condition_match = re.search(r"CLNDN=([^;]+)", info)
            condition = condition_match.group(1).replace("_", " ") if condition_match else "not specified"

            rs_match = re.search(r"RS=(\d+)", info)
            rs_id = f"rs{rs_match.group(1)}" if rs_match else ""

            mc_match = re.search(r"MC=([^;]+)", info)
            mc = mc_match.group(1).replace("|", ", ") if mc_match else ""
            if "," in mc:
                mc1 = mc.split(",")
                mc = mc1[1]

            classification_match = re.search(r"CLNSIG=([^;]+)", info)
            classification = classification_match.group(1).replace("_", " ") if classification_match else "NA"

            variant_match = re.search(r"CLNVC=([^;]+)", info)
            variant_type = variant_match.group(1) if variant_match else "NA"

            if variant_type == "single_nucleotide_variant":
                row = [name, condition, chrom, pos, rs_id, mc, classification]
                data.append(row)
                # Write to TSV file
                outfile.write(
                    f"{name}\t{condition}\t{chrom}\t{pos}\t{rs_id}\t{mc}\t{classification}\n")

    df = pd.DataFrame(data, columns=header)
    return df

def search_variance_in_clinvar(merged_results_df, cinvar_df):
    """
    Adds a 'Clinvar' column to df1 based on chromosome and position matches with df2.

    Args:
        merged_results_df (pd.DataFrame): The first DataFrame containing all analysed positions with heterozygosity or dissimilarities with the reference.
        cinvar_df (pd.DataFrame): The second DataFrame (ClinVar data).

    Returns:
        pd.DataFrame: df1 with the added 'Clinvar' column.
    """
    merged_results_df.columns = [col.lower() for col in merged_results_df.columns]
    cinvar_df.columns = [col.lower() for col in cinvar_df.columns]


    def find_clinvar_info(chrom, pos):
        """Helper function to find ClinVar information."""
        matches = cinvar_df[(cinvar_df['chromosome'] == str(chrom)) & (cinvar_df['position'].astype(str).str.contains(str(pos)))]
        if not matches.empty:
            classifications = matches['classification'].unique().tolist()
            molecular_consequences = matches['molecular_consequence'].unique().tolist()
            if len(classifications) == 1:
                if "Conflicting classifications of pathogenicity" in classifications or "Uncertain significance" in classifications:
                    return "?", classifications[0], "; ".join(molecular_consequences) if molecular_consequences else "NA"
                elif "Benign" not in classifications and "Likely benign" not in classifications:
                    return "+", classifications[0], "; ".join(molecular_consequences) if molecular_consequences else "NA"
                else:
                    return "-", classifications[0], "; ".join(molecular_consequences) if molecular_consequences else "NA"
            elif all("Pathogenic" in c or "Likely pathogenic" in c for c in classifications):
                return "+", "Pathogenic/Likely pathogenic", "; ".join(molecular_consequences) if molecular_consequences else "NA"
            elif all("Benign" in c or "Likely benign" in c for c in classifications):
                return "-", "Benign/Likely benign", "; ".join(molecular_consequences) if molecular_consequences else "NA"
            else:
                return "?", "; ".join(classifications), "; ".join(molecular_consequences) if molecular_consequences else "NA"  # Return all classifications as a string
        else:
            return "-", "NA", "NA"

    results = merged_results_df.apply(lambda row: find_clinvar_info(row['chromosome'], row['position']), axis=1)
    merged_results_df['clinvar'], merged_results_df['classification'], merged_results_df['molecular_consequence'] = zip(*results)
    return merged_results_df


"""
# clinvar vcf file downloaded from ncbi web
cinvar_df = clinvar_to_df("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/clinvar/", "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/clinvar_results_merged.tsv")
merged_result = pd.read_csv("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/merged_results_filtered.tsv", sep="\t")
result = search_variance_in_clinvar(merged_result, cinvar_df)
result.to_csv("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/merged_results_w_clinvar.tsv", sep='\t', index=False)  # Save as TSV
"""

merged_result = pd.read_csv("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/merged_results_filtered.tsv", sep="\t")
filter_vcf_by_genes_grep("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/clinvar/clinvar.vcf.gz", "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/clinvar/filtered_clinvar.vcf", "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/clinvar/gene_list.txt")
cinvar_df = transform_clinvar_vcf("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/clinvar/filtered_clinvar.vcf.gz", "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/clinvar/filtered_clinvar_transformed.tsv")
result = search_variance_in_clinvar(merged_result, cinvar_df)
result.to_csv("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/merged_results_w_clinvar_2.tsv", sep='\t', index=False)

