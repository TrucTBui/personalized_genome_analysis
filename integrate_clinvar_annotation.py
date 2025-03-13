"""
Script used to add clinvar annotation
"""

import pandas as pd

def extract_info_to_df(file, output_tsv, skip_first_line=True):
    """
    Extracts information from a file and returns it as a pandas DataFrame.

    Args:
        file (str): The path to the input file.

    Returns:
        pd.DataFrame: A DataFrame containing the extracted information.
    """
    data = []
    first_line_skipped = False  # Track if the first line has been skipped
    with open(file, "r") as f:
        for line in f:
            if skip_first_line and not first_line_skipped:
                first_line_skipped = True
                continue  # Skip to the next iteration

            fields = line.split("\t")
            if len(fields) >= 16:  # Check if line has enough fields
                name = fields[0]
                canon_name = fields[12]
                condition = fields[3]
                chrom = fields[5]
                pos = fields[6]
                dbSNP_ID = fields[11]
                variant_type = fields[13]
                molecular_consequence = fields[14]
                classification = fields[15]
                data.append({
                    "name": name,
                    "canonical_name": canon_name,
                    "condition": condition,
                    "chrom": chrom,
                    "pos": pos,
                    "dbSNP_ID": dbSNP_ID,
                    "variant_type": variant_type,
                    "molecular_consequence": molecular_consequence,
                    "classification": classification,
                })
            else:
                print(f"Warning: Skipping line with insufficient fields: {line.strip()}") #add a warning if a line is skipped.

    df = pd.DataFrame(data)
    df.to_csv(output_tsv, sep='\t', index=False)  # Save as TSV

extract_info_to_df("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/clinvar/clinvar_result_CCR5.tsv", "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/clinvar/clinvar_result_CCR5_filtered.tsv")