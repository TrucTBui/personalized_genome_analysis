import pandas as pd
import subprocess
import argparse
import os

argparser = argparse.ArgumentParser(description="Transform a variant table into a format suitable for further analysis.")
argparser.add_argument("-i", type=str, required=True, help="Path to the variant table")

args = argparser.parse_args()
variant_table_path = args.i
output_dir = os.path.dirname(variant_table_path)

#variant_table_path = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results/chr17/ENSG00000141480/child/variants.tsv"
#output_path = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/test_results.tsv"
reference_genome = "/mnt/raidinput/input/own/ReferenceGenomes/human_g1k_v37.fasta.gz"

def clean_variant_table(df):
    #Before: Gene	Chromosome	Position	Type	Person	Reference_Base	Frequency_TSAB2165	Frequency_56001808050381	Coverage	Final_Base	Consistency	Special_Case	Imputation
    # After: Gene	Chromosome	Position	Region	Person Ref Alt Type Genotype AllelDepth Imputation
    if df.empty:
        return pd.DataFrame(columns=["Gene", "Chromosome", "Position", "Region", "Person", "Ref", "Alt", "Type", "Genotype", "AllelDepth", "Imputation"])
    df_new = df.rename(columns={
        "Reference_Base": "Ref",
        "Final_Base": "Alt",
        "Type": "Region",})
    
    # remove all rows where ref = N
    df_new = df_new[df_new["Ref"] != "N"]

    df_new["Genotype"] = df_new["Special_Case"].apply(lambda x: "1/1" if "Homozygous!=Ref" in x else ("0/1" if "Heterozygous" in x else "-")) 
    # If the Special Case is heterozygous and Alt contains both allels different from the ref, then 
    
    df_new["Genotype"] = df_new.apply(lambda x: "1/2" if x["Genotype"] == "0/1" and x["Ref"] not in x["Alt"] else x["Genotype"], axis=1)
        
    # AllelDepth: The number of alt and ref alleles in the sample. Have to propagete from all columns starting with Frequency, then count
    # E.g. G:6;A:9	G:18;A:10 => 6+18=24;9+10=19
    freq_cols = [col for col in df_new.columns if col.startswith("Frequency")]

    def get_allele_depth(row):
        if row["Genotype"] != "1/2":
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
                                elif allele == row["Alt"]:
                                    alt_count += count
                            except ValueError:
                                continue
            return f"{ref_count},{alt_count}"
        else:
            alt1,alt2 = row["Alt"].split(",")
            alt1_count = 0
            alt2_count = 0

            for col in freq_cols:
                val = row[col]
                if isinstance(val, str):
                    for part in val.split(";"):
                        if ":" in part:
                            allele, count = part.split(":")
                            try:
                                count = int(count)
                                if allele == alt1:
                                    alt1_count += count
                                elif allele == alt2:
                                    alt2_count += count
                            except ValueError:
                                continue
            return f"{alt1_count},{alt2_count}"

    
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
        ref = str(row["Ref"])
        type =[]
        alts = alt.split(",")
        for a in alts:
            if a == ref:
                continue
            if a == "D":
                type.append("Deletion")
            elif len(a) > len(ref):
                type.append("Insertion")
            else:
                type.append("SNP")
        if not type:
            return "NA"
        return ",".join(type)

    df_new['Alt'] = df_new.apply(get_alt_allele, axis=1)
    df_new["AllelDepth"] = df_new.apply(get_allele_depth, axis=1)
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
    df = df[~df.index.isin(df_del.index)]

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

        genotype_counts = {}
        allel_depths_count = {"ref": 0, "alt": 0}
        imputation_counts = {"family": 0, "pangenome": 0, "individual": 0, "propagation":0, "reference": 0, "-":0}
        for row in block:
            # determine genotype: take the most common genotype in the block: take the one with the most votes
            genotype = row["Genotype"]
            if genotype not in genotype_counts:
                genotype_counts[genotype] = 0
            genotype_counts[genotype] += 1

            # determine allel depth, take the average of the allel depth in the block (both ref and alt)
            allel_depths = row["AllelDepth"].split(",")
            if len(allel_depths) == 2:
                allel_depths_count["ref"] += int(allel_depths[0])
                allel_depths_count["alt"] += int(allel_depths[1])

            # determine imputation: take the one with the most votes
            if "|" in row["Imputation"]:
                imputation = row["Imputation"].split("|")[1].split(";")
            else:
                imputation = [row["Imputation"]]
            # add all imputation types to the imputation_counts dictionary
            for imp in imputation:
                imputation_counts[imp] = imputation_counts.get(imp, 0) + 1

        # Get the genotype with the most votes
        most_common_genotype = max(genotype_counts, key=genotype_counts.get)
        # If it is "-", then set it to "1/1"
        if most_common_genotype == "-":
            most_common_genotype = "1/1"
        
        # Calculate average allel depth
        allel_depth_ref = allel_depths_count["ref"] // len(block)
        allel_depth_alt = allel_depths_count["alt"] // len(block)
        allel_depth = f"{allel_depth_ref},{allel_depth_alt}"

        # Determine the most common imputation type
        most_common_imputation = max(imputation_counts, key=imputation_counts.get)
        second_most_common_imputation = sorted(imputation_counts.items(), key=lambda x: x[1], reverse=True)[1][0] if len(imputation_counts) > 1 else ""
        if most_common_imputation == "propagation":
            most_common_imputation += ";" + second_most_common_imputation
        
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
            "Genotype": most_common_genotype,
            "AllelDepth": allel_depth,
            "Imputation": most_common_imputation,
            "Person": block[0]["Person"]
        })

    df_del = pd.DataFrame(merged_data)
    # integrate back to the original df, delete the rows that were merged
    merged_df = pd.concat([df, df_del], ignore_index=True)
    # Sort by Chromosome and Position
    merged_df.sort_values(by=["Chromosome", "Position"], inplace=True)
    merged_df.reset_index(drop=True, inplace=True)

    return merged_df


if __name__ == "__main__":


    #vcf_list = parse_vcf_path_file(vcf_path)
    #print(vcf_list)
    #vcf_list = vcf_list[:1]
    #If the variant table has #No variant positions found, then df = pd.DataFrame(), else read the file
     
    
    with open(variant_table_path, "r") as f:
        lines = f.readlines()
    if lines[0].startswith("#No variant positions found"):
        df = pd.DataFrame(columns=["Gene", "Chromosome", "Position", "Region", "Person", "Ref", "Alt", "Type", "Genotype", "AllelDepth", "Imputation"])
    else:
        df = pd.read_csv(variant_table_path, sep="\t")

    df = clean_variant_table(df)
    df = merge_deletion_blocks(df, reference_genome)    
    df.to_csv(os.path.join(output_dir, "all_variants.tsv"), sep="\t", index=False)
        
