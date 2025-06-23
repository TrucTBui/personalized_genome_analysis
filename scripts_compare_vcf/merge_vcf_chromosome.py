import os
import pandas as pd
import argparse
import time

def merge_vcf_chromosome(chromosome, person):
    base_path = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results/{chromosome}"
    # iterate over all subdirectories 
    all_vcf_data = []
    for subdir in os.listdir(base_path):
        vcf_path = os.path.join(base_path, subdir, person, "all_variants_clinvar_integrated.vcf")
        if not os.path.exists(vcf_path):
            error_message = f"VCF file not found: {vcf_path}"
            print(error_message)
            continue
        all_vcf_data.append(pd.read_csv(vcf_path, sep="\t", header=0, index_col=None, skiprows=range(0, 12)))

    merged_vcf = pd.concat(all_vcf_data, ignore_index=True)
    #merged_vcf = merged_vcf.drop_duplicates(subset=["CHROM", "POS", "REF", "ALT"], keep="first")
    # sort by CHROM and POS. 
    # order to sort: 1,2,3,...,X,Y

    merged_vcf = merged_vcf.sort_values(by=["POS"]).reset_index(drop=True)
    if 'sample' in merged_vcf.columns:
        merged_vcf = merged_vcf.drop(columns=['sample'])
    output_dir = os.path.join("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_VCF", person)
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, f"{chromosome}_merged_variants.tsv.gz")
    merged_vcf.to_csv(output_path, sep="\t", index=False, compression="gzip")

def merge_vcf_chromosome_family(chromosome):
    FAMILY = [
    "grandfather_father", "grandmother_father",  # paternal grandparents
    "father", "child", "mother", "aunt",         # nuclear family
    "grandmother_mother", "grandfather_mother"   # maternal grandparents
    ]
    FAMILY_MALE = ["grandfather_father", "grandfather_mother", "father", "child"]
    for person in FAMILY:
        if chromosome == "chrY" and person not in FAMILY_MALE:
            continue
        start_time = time.time()
        merge_vcf_chromosome(chromosome, person)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Time taken to merge VCF for {person} on {chromosome}: {elapsed_time:.2f} seconds")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge VCF files for a specific chromosome and person.")
    parser.add_argument("-c","--chromosome", type=str, help="Chromosome to merge VCF files for (e.g., 'chr1', 'chr2', ...)")
    
    args = parser.parse_args()
    chromosome = args.chromosome
    
    merge_vcf_chromosome_family(chromosome)
    
