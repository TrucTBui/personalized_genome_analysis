import pandas as pd
import subprocess
import argparse
import os

argparser = argparse.ArgumentParser(description="Integrate VCF files with variant table")
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

def clean_variant_table(df):
    #Before: Gene	Chromosome	Position	Type	Person	Reference_Base	Frequency_TSAB2165	Frequency_56001808050381	Coverage	Final_Base	Consistency	Special_Case	Imputation
    # After: Gene	Chromosome	Position	Region	Person Ref Alt Type Genotype AllelDepth Imputation
    if df.empty:
        return pd.DataFrame(columns=["Gene", "Chromosome", "Position", "Region", "Person", "Ref", "Alt", "Type", "Genotype", "AllelDepth", "Imputation"])
    df_new = df.rename(columns={
        "Reference_Base": "Ref",
        "Final_Base": "Alt",
        "Type": "Region",})
    df_new["Genotype"] = df_new["Special_Case"].apply(lambda x: "1/1" if x == "Homozygous!=Ref" else ("0/1" if x == "Heterozygous" else "-")) 
    # If the Special Case is heterozygous and Alt contains both allels different from the ref, then 1/2
    df_new["Genotype"] = df_new.apply(lambda x: "1/2" if x["Genotype"] == "0/1" and x["Ref"] not in x["Alt"] else x["Genotype"], axis=1)
    # AllelDepth: The number of alt and ref alleles in the sample. Have to propagete from all columns starting with Frequency, then count
    # E.g. G:6;A:9	G:18;A:10 => 6+18=24;9+10=19
    freq_cols = [col for col in df_new.columns if col.startswith("Frequency")]

    def get_allele_depth(row):
        ref = str(row["Ref"])
        ref_count = 0
        alt_count = 0

        for col in freq_cols:
            val = row[col]
            if isinstance(val, str):
                for part in val.split(";"):
                    if ":" in part:
                        allele, count = part.split(":")
                        try:
                            count = int(count)
                            if allele == ref:
                                ref_count += count
                            else:
                                alt_count += count
                        except ValueError:
                            continue
        return f"{ref_count},{alt_count}"
    
    def get_alt_allele(row):
        ref = str(row["Ref"])
        alt = str(row["Alt"])
        if "/" in alt:
            alts = alt.split("/")
            #print(alts)
            if ref in alts:
                alts.remove(ref)
                return alts[0]
            else:
                return ",".join(alts)
        else:
            return alt
    def assign_variant_type(row):
        alt = str(row["Alt"])
        if alt != "D":
            return "SNP"
        else:
            return "-"

    df_new["AllelDepth"] = df_new.apply(get_allele_depth, axis=1)
    df_new['Alt'] = df_new.apply(get_alt_allele, axis=1)
    df_new["Type"] = df_new.apply(assign_variant_type, axis=1)
    # If the table does not have the column "Imputation", then add it
    if "Imputation" not in df_new.columns:
        df_new["Imputation"] = "-"
    df_new = df_new[["Gene", "Chromosome", "Position", "Region", "Person", "Ref", "Alt", "Type", "Genotype", "AllelDepth", "Imputation"]]
    return df_new


def get_reference_bases(chromosome, start_pos, end_pos,
                         reference_genome="/mnt/raidinput/input/own/ReferenceGenomes/human_g1k_v37.fasta.gz"):
    """Fetch the reference bases for a given region (start_pos to end_pos) using samtools."""
    cmd = f"samtools faidx {reference_genome} {chromosome}:{start_pos}-{end_pos}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"samtools faidx failed: {result.stderr}")
    lines = result.stdout.strip().split("\n")[1:]  # Skip header
    return "".join(lines).replace("\n", "").replace("\r", "")

def merge_deletion_blocks(df, reference_genome):
    
    df_del = df[df["Alt"] == "D"]
    if df_del.empty:
        return pd.DataFrame(columns=["Gene", "Chromosome", "Position", "Region", "Person", "Ref", "Alt", "Type", "Genotype", "AllelDepth", "Imputation"])
    df_del = df_del.sort_values(by=["Chromosome", "Position"]).reset_index(drop=True)

    merged_rows = []
    current_block = []

    for _, row in df_del.iterrows():
        is_deletion = row["Alt"] == "D"

        if not is_deletion:
            if current_block:
                merged_rows.append(current_block)
                current_block = []
            continue

        if not current_block:
            current_block = [row]
        else:
            prev_pos = current_block[-1]["Position"]
            if row["Position"] == prev_pos + 1:
                current_block.append(row)
            else:
                merged_rows.append(current_block)
                current_block = [row]

    if current_block:
        merged_rows.append(current_block)

    # Build merged block rows
    merged_data = []
    for block in merged_rows:
        start = block[0]["Position"]
        end = block[-1]["Position"]
        chrom = block[0]["Chromosome"]

        # Get deleted sequence (from reference)
        ref_seq = get_reference_bases(chrom, start-1, end, reference_genome)

        # Get base before deletion
        if start > 1:
            final_base = get_reference_bases(chrom, start - 1, start - 1, reference_genome)
        else:
            final_base = "N"

        # Copy first row and update relevant fields
        merged_data.append({
            "Gene": block[0]["Gene"],
            "Chromosome": block[0]["Chromosome"],
            "Position": start-1,
            #"End": end,
            #"Length": end - start + 1,
            "Region": block[0]["Region"],
            "Ref": ref_seq,
            "Alt": final_base,
            "Type": "Deletion",
            "Genotype": "1/1",
            "AllelDepth": "-",
            "Imputation": block[int(len(block)/2)]["Imputation"],
            "Person": block[0]["Person"]
        })

    return pd.DataFrame(merged_data)

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
            stderr=subprocess.DEVNULL,  # this suppresses all stderr
            text=True
        )
        #vcf_proc = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
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
        
        if row["Type"] == "Deletion":
            return pd.Series(["-", "-", "-" , False],
                             index=["VCF_REF", "VCF_ALT", "VCF_GT" , "Conflict"])   # Ignore deletion, TODO:review
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
                                
    #print(insertion_df)
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

if __name__ == "__main__":


    vcf_list = parse_vcf_path_file(vcf_path)
    #print(vcf_list)
    vcf_list = vcf_list[:1]
    #If the variant table has #No variant positions found, then df = pd.DataFrame(), else read the file
     
    
    with open(variant_table_path, "r") as f:
        lines = f.readlines()
    if lines[0].startswith("#No variant positions found"):
        df = pd.DataFrame(columns=["Gene", "Chromosome", "Position", "Region", "Person", "Ref", "Alt", "Type", "Genotype", "AllelDepth", "Imputation"])
    else:
        df = pd.read_csv(variant_table_path, sep="\t")

    df = clean_variant_table(df)
    # Remove all rows with "D" in the "Final_Base" column
    df_no_del = df[df["Alt"] != "D"]
    df_del = merge_deletion_blocks(df, reference_gennome)
    df_concat = pd.concat([df_del, df_no_del])
    df_concat.sort_values(by=["Chromosome", "Position"], inplace=True)
    #df_concat.to_csv("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/test_results1.tsv", sep="\t", index=False)
    # Integrate VCF
    # region: chrom in Chromosome, pos1 is min Position, pos2 is max Position
    if not df_concat.empty:
        chromosome = df_concat["Chromosome"].iloc[0]
        pos1 = df_concat["Position"].min()
        pos2 = df_concat["Position"].max()
        region = f"{chromosome}:{pos1}-{pos2}"
        updated_table, conflict_table = integrate_multiple_vcfs(df_concat, vcf_list, region)

    else:
        updated_table = pd.DataFrame(columns=["Gene", "Chromosome", "Position", "Region", "Person", "Ref", "Alt", "Type", "Genotype", "AllelDepth", "Imputation"])
        conflict_table = pd.DataFrame(columns=["Gene", "Chromosome", "Position", "Region", "Person", "Ref", "Alt", "Type", "Genotype", "AllelDepth", "Imputation"])
        
    updated_table.to_csv(os.path.join(output_path, "all_variants.tsv"), sep="\t", index=False)
    conflict_table.to_csv(os.path.join(output_path, "conflicting_variants_vcf.tsv"), sep="\t", index=False)