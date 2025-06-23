import pandas as pd
import os
import argparse

def read_protein_coding_genes_list(protein_coding_genes_file):
    # Gene stable ID  Gene name       Chromosome/scaffold name
    df = pd.read_csv(protein_coding_genes_file, sep="\t", header=0)
    # rename: Gene stable ID to Gene_ID, Gene name to Gene_Name, Chromosome/scaffold name to Chromosome
    df.rename(columns={
        "Gene stable ID": "Gene_ID",
        "Gene name": "Gene_Name",
        "Chromosome/scaffold name": "Chromosome"
    }, inplace=True)
    # Sort by Chromosome and Gene_ID
    df.sort_values(by=["Chromosome", "Gene_ID"], inplace=True)
    # Reset index
    df.reset_index(drop=True, inplace=True)
    return df

def make_statistics(protein_coding_genes_df, chromosome, person):
    
    df_chromosome = protein_coding_genes_df[protein_coding_genes_df["Chromosome"] == chromosome]

    results_list = []

    outer_folder = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results/chr{chromosome}/"

    # Iterate all gene IDs in the df_chromosome, find if there are variant files for them
    for gene_id in df_chromosome["Gene_ID"].unique():
        variant_file = os.path.join(outer_folder, gene_id, person, "all_variants.tsv")
        if not os.path.exists(variant_file):
            continue
        df_variants = pd.read_csv(variant_file, sep="\t", header=0)
        results_list.append(df_variants)

    # Concatenate all results into a single DataFrame
    if results_list:
        all_varianst_df = pd.concat(results_list, ignore_index=True)
    else:
        all_varianst_df = pd.DataFrame(columns=["Gene", "Chromosome", "Position", "Region", "Person", "Ref", "Alt", "Type", "Genotype", "AllelDepth", "Imputation"]) 

    # Remove duplicates based on Chromosome, Position, Ref, Alt, and Person
    all_varianst_df.drop_duplicates(subset=["Chromosome", "Position", "Ref", "Alt", "Person"], inplace=True)
    # Reset index
    all_varianst_df.reset_index(drop=True, inplace=True)

    # Chromosome    Person	Gene_Count	All_Variants_Count	Intron_Variants_Count	Exon_Variants_Count	CDS_Variants_Count	UTR_Variants_Count	SNPs_Count	Insertions_Count	Deletions_Count	SNPs_Exon_Count	SNPs_CDS_Count	Deletion_Exon_Count	Deletion_CDS_Count	Insertion_Exon_Count	Insertion_CDS_Count	Frameshift_CDS_Count

    gene_count = all_varianst_df["Gene"].nunique()
    all_variants_count = all_varianst_df.shape[0]
    intron_variants_count = all_varianst_df[all_varianst_df["Region"].str.contains("intron", case=False)].shape[0]
    exon_variants_count = all_varianst_df[all_varianst_df["Region"].str.contains("exon", case=False)].shape[0]
    cds_variants = all_varianst_df[all_varianst_df["Region"].str.contains("CDS", case=False)]
    cds_variants_count = cds_variants.shape[0]
    utr_variants_count = all_varianst_df[all_varianst_df["Region"].str.contains("utr", case=False)].shape[0]
    snps_count = all_varianst_df[all_varianst_df["Type"] == "SNP"].shape[0]
    insertion_count = all_varianst_df[all_varianst_df["Type"] == "Insertion"].shape[0]
    deletion_count = all_varianst_df[all_varianst_df["Type"] == "Deletion"].shape[0]
    snps_exon_count = all_varianst_df[(all_varianst_df["Type"] == "SNP") & (all_varianst_df["Region"].str.contains("exon", case=False))].shape[0]
    snps_cds_count = all_varianst_df[(all_varianst_df["Type"] == "SNP") & (all_varianst_df["Region"].str.contains("CDS", case=False))].shape[0]
    deletion_exon_count = all_varianst_df[(all_varianst_df["Type"] == "Deletion") & (all_varianst_df["Region"].str.contains("exon", case=False))].shape[0]
    deletion_CDS_count = all_varianst_df[(all_varianst_df["Type"] == "Deletion") & (all_varianst_df["Region"].str.contains("CDS", case=False))].shape[0]
    insertion_exon_count = all_varianst_df[(all_varianst_df["Type"] == "Insertion") & (all_varianst_df["Region"].str.contains("exon", case=False))].shape[0]
    insertion_CDS_count = all_varianst_df[(all_varianst_df["Type"] == "Insertion") & (all_varianst_df["Region"].str.contains("CDS", case=False))].shape[0]
    frameshift_cds_variants = all_varianst_df[
            (all_varianst_df['Region'].str.contains("CDS", na=False)) &
            (((all_varianst_df['Type'] == "Insertion") & (all_varianst_df['Alt'].apply(lambda x: (len(x) - 1) % 3 != 0))) |
              (all_varianst_df['Type'] == "Deletion") & (all_varianst_df['Ref'].apply(lambda x: (len(x) - 1) % 3 != 0)))]
    frameshift_cds_count = frameshift_cds_variants.shape[0]
    
    summary_data = {
        "Chromosome": chromosome,
        "Person": person,
        "Gene_Count": gene_count,
        "All_Variants_Count": all_variants_count,
        "Intron_Variants_Count": intron_variants_count,
        "Exon_Variants_Count": exon_variants_count,
        "CDS_Variants_Count": cds_variants_count,
        "UTR_Variants_Count": utr_variants_count,
        "SNPs_Count": snps_count,
        "Insertions_Count": insertion_count,
        "Deletions_Count": deletion_count,
        "SNPs_Exon_Count": snps_exon_count,
        "SNPs_CDS_Count": snps_cds_count,
        "Deletion_Exon_Count": deletion_exon_count,
        "Deletion_CDS_Count": deletion_CDS_count,
        "Insertion_Exon_Count": insertion_exon_count,
        "Insertion_CDS_Count": insertion_CDS_count,
        "Frameshift_CDS_Count": frameshift_cds_count
    }
    # Create a DataFrame from the summary data
    summary_df = pd.DataFrame([summary_data])

    return summary_df, cds_variants, frameshift_cds_variants

        

if __name__ == "__main__":
    protein_coding_genes_df = read_protein_coding_genes_list("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genes/mart_export.txt.gz")

    
    parser = argparse.ArgumentParser(description="Create statistics for protein coding genes.")
    parser.add_argument("-c", "--chromosome", type=str, required=True, help="Chromosome to analyze (e.g., '2').")
    parser.add_argument("-p", "--person", type=str, required=True, help="Person to analyze (e.g., 'child').")
    args = parser.parse_args()
    chromosome = args.chromosome
    person = args.person

    statistics_df, cds_variant_df, frame_shift_df = make_statistics(protein_coding_genes_df, chromosome, person)
    
    cds_variant_folder = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/CDS_variants/{person}/"
    os.makedirs(cds_variant_folder, exist_ok=True)
    statistics_folder = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Statistics_no_dup/{person}/"
    os.makedirs(statistics_folder, exist_ok=True)

    statistics_df.to_csv(os.path.join(statistics_folder, f"stats_{chromosome}_{person}.tsv"), sep="\t", index=False)
    #cds_variant_df.to_csv(os.path.join(cds_variant_folder, f"cds_variants_{chromosome}_{person}.tsv"), sep="\t", index=False)
    #frame_shift_df.to_csv(os.path.join(cds_variant_folder, f"frameshift_cds_variants_{chromosome}_{person}.tsv"), sep="\t", index=False)
    