"""
Version 11.03.2025
"""

import pandas as pd
import re
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")

def compute_concordance(result_pos_wise_df):
    """Calculate the concordance rate of each family member with the reference genome."""
    members = result_pos_wise_df.columns[5:]
    total_positions = len(result_pos_wise_df)
    concordance = {}
    num_variants = {}
    for member in members:
        matches = (result_pos_wise_df[member] == result_pos_wise_df["Reference"]).sum()
        concordance[member] = round(float(matches / total_positions), 5)  # Convert to Python float
        num_variants[member] = total_positions - matches

    total_pos_row = {"Metric": "Total_Position"}
    concordance_row = {"Metric": "Similarities_to_Reference"}
    variants_row = {"Metric": "Dissimilarities_Count"}

    for member in members:
        total_pos_row[member] = total_positions
        concordance_row[member] = concordance[member]
        variants_row[member] = num_variants[member]

    return pd.DataFrame([total_pos_row, concordance_row, variants_row])

def heterozygosity_rate(df):
    """Compute the heterozygosity rate for each individual."""
    members = df.columns[5:]
    heterozygosity = {}
    for member in members:
        heter_count = df[member].apply(lambda x: bool(re.match(r"^[ATGC]/[ATGC]$", str(x)))).sum()
        heterozygosity[member] = round(float(heter_count / len(df)), 5)  # Convert to Python float
    return heterozygosity

def heterozygosity_rate_and_count(df):
    """Compute the heterozygosity rate and count for each individual."""
    members = df.columns[5:]
    results = {}
    for member in members:
        heter_count = df[member].apply(lambda x: bool(re.match(r"^[ATGC]/[ATGC]$", str(x)))).sum()
        heterozygosity = round(float(heter_count / len(df)), 5)  # Convert to Python float
        results[member] = {"rate": heterozygosity, "count": heter_count}

    rate_row = {"Metric": "Heterozygosity_Rate"}
    count_row = {"Metric": "Heterozygosity_Count"}

    for member in members:
        rate_row[member] = results[member]["rate"]
        count_row[member] = results[member]["count"]

    return pd.DataFrame([rate_row, count_row])

def heterozygosity_rate_by_type(df, type_column="Type", genotype_columns=None):
    """
    Compute the heterozygosity rate for each individual, broken down by gene type.

    Args:
        df (pd.DataFrame): The input DataFrame.
        type_column (str): The name of the column containing gene types.
        genotype_columns (list): A list of column names containing genotype data.
                                If None, columns starting from index 5 are used.

    Returns:
        dict: A dictionary where keys are gene types and values are dictionaries
              containing heterozygosity rates for each individual.
    """

    if genotype_columns is None:
        genotype_columns = df.columns[5:]

    gene_types = df[type_column].unique()
    results = []

    for gene_type in gene_types:
        subset_df = df[df[type_column] == gene_type]
        for member in genotype_columns:
            heter_count = subset_df[member].apply(lambda x: bool(re.match(r"^[ATGC]/[ATGC]$", str(x)))).sum()
            heterozygosity = round(float(heter_count / len(subset_df)), 5) if len(subset_df) > 0 else 0.0
            results.append({"type": gene_type, "individual": member, "heterozygosity_rate": heterozygosity})

    heterozygosity_df = pd.DataFrame(results)
    reshaped_df = heterozygosity_df.pivot(index="type", columns="individual", values="heterozygosity_rate")

    # Rename the index
    reshaped_df.index = "Hetero_" + reshaped_df.index.astype(str)

    # Reset the index and rename the first column
    reshaped_df = reshaped_df.reset_index()
    reshaped_df = reshaped_df.rename(columns={'type': 'Metric'})

    return reshaped_df

def heterozygosity_rate_and_count_by_type(df, type_column="Type", genotype_columns=None):
    """
    Compute the heterozygosity rate and count for each individual, broken down by gene type.

    Args:
        df (pd.DataFrame): The input DataFrame.
        type_column (str): The name of the column containing gene types.
        genotype_columns (list): A list of column names containing genotype data.
                                 If None, columns starting from index 5 are used.

    Returns:
        pd.DataFrame: A DataFrame with heterozygosity rates and counts for each individual and gene type.
    """

    if genotype_columns is None:
        genotype_columns = df.columns[5:]

    gene_types = df[type_column].unique()
    results = []

    for gene_type in gene_types:
        subset_df = df[df[type_column] == gene_type]
        for member in genotype_columns:
            heter_count = subset_df[member].apply(lambda x: bool(re.match(r"^[ATGC]/[ATGC]$", str(x)))).sum()
            heterozygosity = round(float(heter_count / len(subset_df)), 5) if len(subset_df) > 0 else 0.0
            results.append({"type": gene_type, "individual": member, "heterozygosity_rate": heterozygosity, "heterozygosity_count": heter_count})

    heterozygosity_df = pd.DataFrame(results)

    # Pivot for rates
    reshaped_rate_df = heterozygosity_df.pivot(index="type", columns="individual", values="heterozygosity_rate")
    reshaped_rate_df = reshaped_rate_df.reset_index()
    reshaped_rate_df['type'] = "Hetero_rate_" + reshaped_rate_df['type'].astype(str)
    reshaped_rate_df = reshaped_rate_df.rename(columns={'type': 'Metric'})

    # Pivot for counts
    reshaped_count_df = heterozygosity_df.pivot(index="type", columns="individual", values="heterozygosity_count")
    reshaped_count_df = reshaped_count_df.reset_index()
    reshaped_count_df['type'] = "Hetero_count_" + reshaped_count_df['type'].astype(str)
    reshaped_count_df = reshaped_count_df.rename(columns={'type': 'Metric'})

    # Merge rate and count DataFrames
    merged_df = pd.concat([reshaped_rate_df, reshaped_count_df], ignore_index=True)

    return merged_df

df = pd.read_csv("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/merged_results.tsv", sep="\t")


concordance_df = compute_concordance(df)
heterozygosity_df = heterozygosity_rate_and_count(df)
heterozygosity_by_type = heterozygosity_rate_and_count_by_type(df)

# Merge the DataFrames
merged_df = pd.concat([concordance_df, heterozygosity_df, heterozygosity_by_type], ignore_index=True)
#merged_df.to_csv("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/CCR5/quality_check.tsv", sep="\t", index=False)
merged_df.to_csv("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Family/quality_check.tsv", sep="\t", index=False)
