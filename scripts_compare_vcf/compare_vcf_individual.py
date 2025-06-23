# Unused code!
import pandas as pd
import subprocess
import argparse
import os

argparser = argparse.ArgumentParser(description="Compare variant table with pre-computed VCF files.")
argparser.add_argument("-i", type=str, required=True, help="Path to the variant table")
argparser.add_argument("-v", type=str, required=True, help="Path to the file containing VCF paths")
argparser.add_argument("-o", type=str, required=False, help="Path to the output folder")

args = argparser.parse_args()
variant_table_path = args.i
vcf_path = args.v
output_path = args.o

#variant_table_path = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results/chr17/ENSG00000141480/child/variants.tsv"
#vcf_path = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genomes/VCF/child.txt"
#output_path = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/test_results.tsv"
reference_gennome = "/mnt/raidinput/input/own/ReferenceGenomes/human_g1k_v37.fasta.gz"

def parse_vcf_path_file(vcf_path):
    # tsv file, 1st col: indel vcf, 2nd col: snp vcf
    with open(vcf_path, "r") as f:
        lines = f.readlines()
    vcf_list = []
    for line in lines:
        if line.startswith("#"):
            continue
        line = line.strip().split()
        #print(line)
        indel_vcf = line[0]
        snp_vcf = line[1]
        vcf_list.append((snp_vcf, indel_vcf))
    return vcf_list


def integrate_vcf_region(finding_table, snp_vcf, indel_vcf, region):
    """
    Read variants from a VCF region using bcftools and integrates insertions,
    while reporting conflicting SNPs with finding_table.
    """
    person = finding_table["Person"].iloc[0]
    def parse_vcf(vcf_path, region):
        # If the vcf_path starts with 56, then add chr to the region. Since it is the other set of seq data
        if os.path.basename(vcf_path).startswith("56"):
            region = f"chr{region}"
        
        command = ["/home/b/buit/miniconda3/envs/HiWi/bin/bcftools", "view", "-r", region, vcf_path]
        vcf_proc = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,  
            text=True
        )
        variants = {}
        for line in vcf_proc.stdout:
            if line.startswith("#"):
                continue
            #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  #SAMPLE_NAME
            chrom, pos, _, ref, alt, *_, genotype = line.strip().split("\t")
            pos = int(pos)
            GT = genotype.split(":")[0]
            if len(ref) < len(alt):
                variant_type = "insertion"
            elif len(ref) > len(alt):
                variant_type = "deletion"
            else:
                variant_type = "SNP"
            
            variants[pos] = {
                "Chromosome": chrom,
                "Position": pos,
                "VCF_REF": ref,
                "VCF_ALT": alt,
                "VCF_GT": GT,
                "VCF_TYPE": variant_type
            }
        vcf_proc.stdout.close()
        vcf_proc.wait()
        return variants

    # Parse both VCFs
    snp_variants = parse_vcf(snp_vcf, region)
    indel_variants = parse_vcf(indel_vcf, region)

    # Annotate finding_table using SNPs
    def compare_snps(row):
        var = snp_variants.get(row["Position"])
        
        #if row["Type"] == "Deletion":
        #    return pd.Series(["-", "-", "-" , False],
        #                     index=["VCF_REF", "VCF_ALT", "VCF_GT" , "Conflict"])   # Ignore deletion, TODO:review
        if not var:
            return pd.Series(["-", "-", "-" , True],
                         index=["VCF_REF", "VCF_ALT", "VCF_GT" , "Conflict"])
        match = (row["Ref"] == var["VCF_REF"]) and (row["Alt"] == var["VCF_ALT"])
        # The genotype also needs to match. Split by "/" or "|"
        vcf_gt = var["VCF_GT"].split("/") if "/" in var["VCF_GT"] else var["VCF_GT"].split("|")
        finding_gt = row["Genotype"].split("/") 
        # if they contain the same alleles, they are equal
        match = match and (set(vcf_gt) == set(finding_gt))
        conflict = not match
        return pd.Series([var["VCF_REF"], var["VCF_ALT"], var['VCF_GT'] , conflict],
                         index=["VCF_REF", "VCF_ALT", "VCF_GT" , "Conflict"])

    annotated = finding_table.copy()
    annotated[["VCF_REF", "VCF_ALT", "VCF_GT" , "Conflict"]] = annotated.apply(compare_snps, axis=1)

    def compare_indels(row):
        pass  # TODO
    """
    # Extract insertions from indel VCF not in table
    if len(indel_variants) > 0:
        insertions = [v for k, v in indel_variants.items() if v["VCF_TYPE"] == "insertion"]
    #Gene	Chromosome	Position	Region	Ref	Alt	Type	Genotype	AllelDepth	Imputation	Person
        insertion_df = pd.DataFrame({"Gene": finding_table["Gene"].iloc[0],
                                    "Chromosome": finding_table["Chromosome"].iloc[0],
                                    "Position": [v["Position"] for v in insertions],
                                    "Region": finding_table["Region"].iloc[0],
                                    "Ref": [v["VCF_REF"] for v in insertions],
                                    "Alt": [v["VCF_ALT"] for v in insertions],
                                    "Type": ["Insertion" for _ in insertions],
                                    "Genotype": [v["VCF_GT"] for v in insertions],
                                    "AllelDepth": ["-" for _ in insertions],
                                    "Imputation": ["-" for _ in insertions],
                                    "Person": finding_table["Person"].iloc[0]
        })
    else:
        insertion_df = pd.DataFrame(columns=["Gene", "Chromosome", "Position", "Region", "Ref", "Alt", "Type", "Genotype", "AllelDepth", "Imputation", "Person"])
                                
    """
    # Final outputs
    conflict_table = annotated[annotated["Conflict"] == True].copy()
    updated_table = pd.concat([annotated.drop(columns=["Conflict", "VCF_REF", "VCF_ALT", "VCF_GT"]), insertion_df], ignore_index=True)
    updated_table.sort_values(by=["Chromosome", "Position"], inplace=True)
    updated_table.reset_index(drop=True, inplace=True)

    return updated_table, conflict_table

def integrate_multiple_vcfs(finding_table, vcf_list, region):
    """
    Integrate variants from multiple VCF file pairs (snp_vcf, indel_vcf) for a given region.

    Returns:
    - final_merged_table: all integrated and merged findings
    - all_conflicts: all conflicting SNP rows from all replicates
    """

    merged_tables = []
    conflict_tables = []

    for (snp_vcf, indel_vcf) in vcf_list:
        updated_table, conflict_table = integrate_vcf_region(finding_table, snp_vcf, indel_vcf, region)
        merged_tables.append(updated_table)
        conflict_tables.append(conflict_table)

    final_merged_table = pd.concat(merged_tables, ignore_index=True)
    all_conflicts = pd.concat(conflict_tables, ignore_index=True)
    final_merged_table.sort_values(by=["Chromosome", "Position"], inplace=True)
    final_merged_table.drop_duplicates(subset=["Chromosome", "Position", "Type"], keep="first", inplace=True)
    final_merged_table.reset_index(drop=True, inplace=True)
    all_conflicts.sort_values(by=["Chromosome", "Position"], inplace=True)
    all_conflicts.drop_duplicates(subset=["Chromosome", "Position", "Type","VCF_GT"], keep="first", inplace=True)
    all_conflicts.reset_index(drop=True, inplace=True)

    return final_merged_table, all_conflicts
