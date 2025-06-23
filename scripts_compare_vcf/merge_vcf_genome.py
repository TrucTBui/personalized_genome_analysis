import os
import pandas as pd

def merge_vcf_genome(person):
    FAMILY_MALE = ["grandfather_father", "grandfather_mother", "father", "child"]

    base_path = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_VCF/{person}"
    chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    all_vcf_data = []
    for chromosome in chromosomes:
        if chromosome == "chrY" and person not in FAMILY_MALE:
            continue
        vcf_path = os.path.join(base_path, f"{chromosome}_merged_variants.tsv.gz")
        if not os.path.exists(vcf_path):
            error_message = f"VCF file not found: {vcf_path}"
            print(error_message)
            continue
        vcf_data = pd.read_csv(vcf_path, sep="\t", header=0, index_col=None)
        all_vcf_data.append(vcf_data)
    merged_vcf = pd.concat(all_vcf_data, ignore_index=True)
    merged_vcf = merged_vcf.reset_index(drop=True)
    return merged_vcf

def write_final_vcf(merged_vcf, person):
    output_dir = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_VCF/{person}"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "all_variants.vcf")
    # Add these rows as headers:
    ##fileformat=VCFv4.2
##INFO=<ID=Gene,Number=1,Type=String,Description="Ensembl Gene ID">
##INFO=<ID=Imputation,Number=1,Type=String,Description="Imputation status">
##INFO=<ID=Variant_type,Number=1,Type=String,Description="Variant Type">
##INFO=<ID=Clinvar_ID,Number=1,Type=String,Description="ClinVar Variation ID">
##INFO=<ID=dbSNP_ID,Number=1,Type=String,Description="dbSNP RS ID">
##INFO=<ID=Condition,Number=1,Type=String,Description="Associated Condition">
##INFO=<ID=Molecular_consequence,Number=1,Type=String,Description="Molecular Consequence">
##INFO=<ID=Classification,Number=1,Type=String,Description="Variant Classification">
##INFO=<ID=Region_type,Number=1,Type=String,Description="Variant Region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=String,Description="Allele depth for each allele (REF, ALT1, ALT2, ...)">
    vcf_metadata_lines = ["##fileformat=VCFv4.2",
                          "##INFO=<ID=Gene,Number=1,Type=String,Description=\"Ensembl Gene ID\">",
                          "##INFO=<ID=Imputation,Number=1,Type=String,Description=\"Imputation status\">",
                          "##INFO=<ID=Variant_type,Number=1,Type=String,Description=\"Variant Type\">",
                          "##INFO=<ID=Clinvar_ID,Number=1,Type=String,Description=\"ClinVar Variation ID\">",
                          "##INFO=<ID=dbSNP_ID,Number=1,Type=String,Description=\"dbSNP RS ID\">",
                          "##INFO=<ID=Condition,Number=1,Type=String,Description=\"Associated Condition\">",
                          "##INFO=<ID=Molecular_consequence,Number=1,Type=String,Description=\"Molecular Consequence\">",
                          "##INFO=<ID=Classification,Number=1,Type=String,Description=\"Variant Classification\">",
                          "##INFO=<ID=Region_type,Number=1,Type=String,Description=\"Variant Region\">",
                          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                          "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depth for each allele (REF, ALT1, ALT2, ...)\">"]
    with open(output_path, 'w') as f:
            for line in vcf_metadata_lines:
                f.write(line.strip() + '\n')
            merged_vcf.to_csv(f, sep="\t", index=False, header=True, compression=None) 
    # bgzip the file
    os.system(f"bgzip -c {output_path} > {output_path}.gz")
    # delete the original file
    os.remove(output_path)
    # tabix index the file
    os.system(f"tabix -p vcf {output_path}.gz")

if __name__ == "__main__":
    FAMILY = [
    "grandfather_father", "grandmother_father",  # paternal grandparents
    "father", "child", "mother", "aunt",         # nuclear family
    "grandmother_mother", "grandfather_mother"   # maternal grandparents
    ]
    for person in FAMILY:
        merged_vcf = merge_vcf_genome(person)
        if merged_vcf.empty:
            print(f"No VCF data found for {person}. Skipping writing final VCF.")
            continue
        write_final_vcf(merged_vcf, person)
        print(f"Final VCF written for {person}.")