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

def make_statistics(protein_coding_genes_df, chromosome, gene_id_chunk_file, person):
    # /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genes/gene_id/gene_id_list_chr1_chunk1.txt
    chunk_name = os.path.basename(gene_id_chunk_file).replace(".txt", "").replace("gene_id_list_", "")

    def read_gene_id_chunk_list(gene_id_chunk_file):
        with open(gene_id_chunk_file, 'r') as f:
            gene_id_chunk_list = [line.strip() for line in f if line.strip()]
        return gene_id_chunk_list

    gene_chunk_list = read_gene_id_chunk_list(gene_id_chunk_file)
    # Filter the DataFrame for the given chromosome
    df_chromosome = protein_coding_genes_df[protein_coding_genes_df["Chromosome"] == chromosome]
    #print(df_chromosome.shape)
    # Continue filtering for the gene IDs in the gene_chunk_list
    df_chromosome = df_chromosome[df_chromosome["Gene_ID"].isin(gene_chunk_list)]
    #print(df_chromosome.shape)

    outer_folder = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results/chr{chromosome}/"

    all_variants_count = 0
    intron_variants_count = 0
    exon_variants_count = 0
    cds_variants_count = 0
    utr_variants_count = 0
    snps_count = 0
    insertion_count = 0
    deletion_count = 0
    deletion_exon_count = 0
    deletion_CDS_count = 0
    insertion_exon_count = 0
    insertion_CDS_count = 0
    snps_exon_count = 0
    snps_cds_count = 0
    frameshift_cds_count = 0

    all_cds_variants = []
    all_frameshift_cds_variants = []
    # iterate over the genes in the chromosome
    for _, row in df_chromosome.iterrows():
        gene_id = row["Gene_ID"]
        gene_name = row["Gene_Name"]
        variant_file = os.path.join(outer_folder, gene_id, person, "all_variants.tsv")
        if not os.path.exists(variant_file):
            print(f"Variant file does not exist for gene {gene_id} ({gene_name}) in chromosome {chromosome}. Skipping.")
            continue

        # Read the variant file
        df_variants = pd.read_csv(variant_file, sep="\t", header=0)

        # add all CDS variants to the list
        all_cds_variants.extend(df_variants[df_variants['Region'].str.contains("CDS", na=False)].values.tolist())
        # add all frameshift CDS variants to the list: CDS, insertions: len(alt-1) % 3 != 0, or deletions: len(ref-1) % 3 != 0
        frameshift_cds_variants = df_variants[
            (df_variants['Region'].str.contains("CDS", na=False)) &
            (((df_variants['Type'] == "Insertion") & (df_variants['Alt'].apply(lambda x: (len(x) - 1) % 3 != 0))) |
              (df_variants['Type'] == "Deletion") & (df_variants['Ref'].apply(lambda x: (len(x) - 1) % 3 != 0)))]

        all_frameshift_cds_variants.extend(frameshift_cds_variants.values.tolist())
        frameshift_cds_count += len(frameshift_cds_variants)


        # Gene	Chromosome	Position	Region	Person	Ref	Alt	Type	Genotype	AllelDepth	Imputation
        all_variants_count += len(df_variants)

        # Count variants in different regions
        intron_variants_count += df_variants['Region'].str.contains("intron", na=False).sum()
        exon_variants_count += df_variants['Region'].str.contains("exon", na=False).sum()
        cds_variants_count += df_variants['Region'].str.contains("CDS", na=False).sum()
        utr_variants_count += df_variants['Region'].str.contains("utr", na=False).sum()

        # Count SNPs, insertions, and deletions
        snps_count += (df_variants["Type"] == "SNP").sum()
        insertion_count += (df_variants["Type"] == "Insertion").sum()
        deletion_count += (df_variants["Type"] == "Deletion").sum()

        # Count deletions and insertions in introns, exons, and CDS
        deletion_exon_count += ((df_variants["Type"] == "Deletion") & (df_variants["Region"].str.contains("exon", na=False))).sum()
        deletion_CDS_count += ((df_variants["Type"] == "Deletion") & (df_variants["Region"].str.contains("CDS", na=False))).sum()
        insertion_exon_count += ((df_variants["Type"] == "Insertion") & (df_variants["Region"].str.contains("exon", na=False))).sum()
        insertion_CDS_count += ((df_variants["Type"] == "Insertion") & (df_variants["Region"].str.contains("CDS", na=False))).sum()

        # Count SNPs in exons and CDS
        snps_exon_count += ((df_variants["Type"] == "SNP") & (df_variants["Region"].str.contains("exon", na=False))).sum()
        snps_cds_count += ((df_variants["Type"] == "SNP") & (df_variants["Region"].str.contains("CDS", na=False))).sum()
    # Create a summary DataFrame
    summary_data = {
        "Chromosome": chromosome,
        "Chunk": chunk_name,
        "Person": person,
        "Gene_Count": len(df_chromosome),

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
    summary_df = pd.DataFrame([summary_data])
    all_cds_variants_df = pd.DataFrame(all_cds_variants, columns=["Gene", "Chromosome", "Position", "Region", "Person", "Ref", "Alt", "Type", "Genotype", "AllelDepth", "Imputation"])
    all_frameshift_cds_variants_df = pd.DataFrame(all_frameshift_cds_variants, columns=["Gene", "Chromosome", "Position", "Region", "Person", "Ref", "Alt", "Type", "Genotype", "AllelDepth", "Imputation"])
    return summary_df, all_cds_variants_df, all_frameshift_cds_variants_df

if __name__ == "__main__":
    protein_coding_genes_df = read_protein_coding_genes_list("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genes/mart_export.txt.gz")

    
    parser = argparse.ArgumentParser(description="Create statistics for protein coding genes.")
    parser.add_argument("-c", "--chromosome", type=str, required=True, help="Chromosome to analyze (e.g., '2').")
    parser.add_argument("-g", "--gene_list_chunk_file", type=str, required=True, help="Path to the gene ID chunk file.")
    parser.add_argument("-p", "--person", type=str, required=True, help="Person to analyze (e.g., 'child').")
    args = parser.parse_args()
    chromosome = args.chromosome
    gene_list_chunk_file = args.gene_list_chunk_file
    person = args.person
    chunk_name = os.path.basename(gene_list_chunk_file).replace(".txt", "").replace("gene_id_list_", "")

    
    statistics_df, cds_variant_df, frame_shift_df = make_statistics(protein_coding_genes_df, chromosome, gene_list_chunk_file, person)
    
    cds_variant_folder = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/CDS_variants/{person}/"
    os.makedirs(cds_variant_folder, exist_ok=True)
    statistics_folder = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Statistics/{person}/"
    os.makedirs(statistics_folder, exist_ok=True)

    statistics_df.to_csv(os.path.join(statistics_folder, f"stats_{chromosome}_{chunk_name}_{person}.tsv"), sep="\t", index=False)

""" cds_variant_df.to_csv(os.path.join(cds_variant_folder, f"cds_variants_{chromosome}_{chunk_name}_{person}.tsv"), sep="\t", index=False)
    frame_shift_df.to_csv(os.path.join(cds_variant_folder, f"frameshift_cds_variants_{chromosome}_{chunk_name}_{person}.tsv"), sep="\t", index=False)
"""