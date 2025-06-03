import pandas as pd
import subprocess
import re
import argparse



def read_variant_file(variant_file):
    #Columns: Gene	Chromosome	Position	Region	Ref	Alt	Type	Genotype	AllelDepth	Imputation	Person

    variants_df = pd.read_csv(variant_file, sep='\t')
    # change column name from AllelDepth to AlleleDepth
    if 'AllelDepth' in variants_df.columns:
        variants_df.rename(columns={'AllelDepth': 'AlleleDepth'}, inplace=True)
    return variants_df

def read_clinvar_region(clinvar_vcf, region):
    command = f"/home/b/buit/miniconda3/envs/HiWi/bin/bcftools view -r {region} {clinvar_vcf}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Error reading ClinVar VCF: {result.stderr}")
    # Parse the output into a DataFrame
    clinvar_data = []
    for line in result.stdout.splitlines():
        if line.startswith("#"):
            continue
        fields = line.strip().split()
        if len(fields) < 8:
            continue  # Skip malformed lines

        chrom = fields[0]  # Chromosome
        pos = fields[1]  # Position
        ref = fields[3]
        alt = fields[4]
        info = fields[7]  # INFO field

        # Extract relevant fields from INFO
        name_match = re.search(r"CLNHGVS=([^;]+)", info)
        name = name_match.group(1) if name_match else "NA"

        condition_match = re.search(r"CLNDN=([^;]+)", info)
        condition = condition_match.group(1) if condition_match else "not_provided"
        if condition == "not provided":
            condition = "not_provided"

        rs_match = re.search(r"RS=(\d+)", info)
        rs_id = f"rs{rs_match.group(1)}" if rs_match else ""

        mc_match = re.search(r"MC=([^;]+)", info)
        mc = mc_match.group(1) if mc_match else ""
        if "|" in mc:
            mc1 = mc.split("|")
            mc = mc1[1].strip()

        classification_match = re.search(r"CLNSIG=([^;]+)", info)
        classification = classification_match.group(1) if classification_match else "NA"

        classification_conflict_match = re.search(r"CLNSIGCONF=([^;]+)", info)
        if classification_conflict_match:
            # Add to classification if it exists
            classification += f"({classification_conflict_match.group(1)})"

        variant_match = re.search(r"CLNVC=([^;]+)", info)
        #variant_type = variant_match.group(1) if variant_match else "NA"
        clinvar_data.append({
            "Chromosome": chrom,
            "Position": pos,
            "Ref": ref,
            "Alt": alt,
            "Var_ID": name,
            "RS_ID": rs_id,
            "Condition": condition,
            "MC": mc,
            "Classification": classification
            #"Variant_Type": variant_type
        })
    if not clinvar_data:
        clinvar_df = pd.DataFrame(columns=["Chromosome", "Position", "Ref", "Alt", "Var_ID", "RS_ID", "Condition", "MC", "Classification"])
    else:
        clinvar_df = pd.DataFrame(clinvar_data)
    return clinvar_df

def integrate_clinvar(variants_df, clinvar_df):
    variants_df['Chromosome'] = variants_df['Chromosome'].astype(str)
    variants_df['Position'] = variants_df['Position'].astype(str)
    clinvar_df['Chromosome'] = clinvar_df['Chromosome'].astype(str)
    clinvar_df['Position'] = clinvar_df['Position'].astype(str)

    # Merge the two DataFrames on Chromosome and Position
    merged_df = pd.merge(variants_df, clinvar_df, on=["Chromosome", "Position", "Ref", "Alt"], how="left", suffixes=("", "_clinvar"))

    # Fill NaN values in the ClinVar columns with "NA"
    #clinvar_columns = ["Var_ID", "RS_ID", "Condition", "MC", "Classification", "Variant_Type"]
    clinvar_columns = ["Var_ID", "RS_ID", "Condition", "MC", "Classification"]
    for col in clinvar_columns:
        merged_df[col] = merged_df[col].fillna("NA")
    
    return merged_df

def write_vcf_from_dataframe(df, output_file):
    if not df.empty:
        sample_name = df['Person'].iloc[0] if 'Person' in df.columns else 'sample'
    else:
        sample_name = 'sample'
    with open(output_file, 'w') as f:
        # Write VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write('##INFO=<ID=Gene,Number=1,Type=String,Description="Ensembl Gene ID">\n')
        f.write('##INFO=<ID=Imputation,Number=1,Type=String,Description="Imputation status">\n')
        f.write('##INFO=<ID=Variant_type,Number=1,Type=String,Description="Variant Type">\n')
        f.write('##INFO=<ID=Clinvar_ID,Number=1,Type=String,Description="ClinVar Variation ID">\n')
        f.write('##INFO=<ID=dbSNP_ID,Number=1,Type=String,Description="dbSNP RS ID">\n')
        f.write('##INFO=<ID=Condition,Number=1,Type=String,Description="Associated Condition">\n')
        f.write('##INFO=<ID=Molecular_consequence,Number=1,Type=String,Description="Molecular Consequence">\n')
        f.write('##INFO=<ID=Classification,Number=1,Type=String,Description="Variant Classification">\n')
        f.write('##INFO=<ID=Region_type,Number=1,Type=String,Description="Variant Region">\n')

        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=AD,Number=R,Type=String,Description="Allele depth for each allele (REF, ALT1, ALT2, ...)">\n')

        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")

        for _, row in df.iterrows():
            chrom = row['Chromosome']
            pos = row['Position']
            ref = row['Ref']
            alt = row['Alt']

           
            gene = row.get('Gene', '.')
            region_type = row.get('Region', '.').replace(";",",")
            imputation = row.get('Imputation', '.')
            variant_type = row.get('Type', '.') 
            
            clinvar_id = row.get('Var_ID', '.')
            dbsnp_id = row.get('RS_ID', '.')
            condition = row.get('Condition', '.')
            mc = row.get('MC', '.')
            classification = row.get('Classification', '.')

            info_parts = [
                f"Gene={gene}",
                f"Region_type={region_type}",
                f"Imputation={imputation}",
                f"Variant_type={variant_type}",
                f"Clinvar_ID={clinvar_id}",
                f"dbSNP_ID={dbsnp_id}",
                f"Condition={condition}",
                f"Molecular_consequence={mc}",
                f"Classification={classification}"
            ]
            # Filter out parts with missing values if represented by '.', or ensure NA/None are handled
            info = ";".join(part for part in info_parts if part.split('=')[1] not in ['.', 'NA', 'None'])
            
            info = info.replace(' ', '_')
            if not info: # If info string becomes empty
                info = "."

            
            format_field = "GT:AD"
            sample_genotype = str(row.get('Genotype', './.')) 
            sample_allele_depth = str(row.get('AlleleDepth', '.'))
            sample_data_value = f"{sample_genotype}:{sample_allele_depth}"

            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}\t{format_field}\t{sample_data_value}\n")

if __name__ == "__main__":

    argparser = argparse.ArgumentParser(description="Integrate variant data with ClinVar information.")
    argparser.add_argument("-v","--variant_file", type=str, required=True, help="Path to the variant file (TSV format).")
    argparser.add_argument("-c","--clinvar_vcf", type=str, required=False, default='/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/clinvar/clinvar_20250521.vcf.gz', help="Path to the ClinVar VCF file (gzipped).")

    args = argparser.parse_args()
    variant_file = args.variant_file
    clinvar_vcf = args.clinvar_vcf
    vcf_output_file = args.variant_file.replace('.tsv', '_clinvar_integrated.vcf')


    variants_df = read_variant_file(variant_file)
    #print("Variants DataFrame:")
    #print(variants_df.head())
    # Region = {Chromosome}:{min_position}-{max_position} in the varians_df
    if not variants_df.empty:
        region = f"{variants_df['Chromosome'].iloc[0]}:{variants_df['Position'].min()}-{variants_df['Position'].max()}"

        clinvar_df = read_clinvar_region(clinvar_vcf, region)
        #print("ClinVar DataFrame:")
        #print(clinvar_df.head())

        integrated_df = integrate_clinvar(variants_df, clinvar_df)
        #print("Integrated DataFrame:")
        #print(integrated_df.head())
        # If there is pathogenic variant in the integrated DataFrame print that out
        if not integrated_df.empty:
            pathogenic_variants = integrated_df[integrated_df['Classification'].str.contains(r"Likely_pathogenic|Pathogenic", case=True, na=False)]
            if not pathogenic_variants.empty:
                # Print each pathogenic variant in a separate line, sep by tab
                for _, row in pathogenic_variants.iterrows():
                    print(f"{row['Gene']}\t{row['Chromosome']}\t{row['Position']}\t{row['Ref']}\t{row['Alt']}\t{row['Classification']}\t{row.get('Person', 'Unknown')}")

        # Save the integrated DataFrame to a new file
        #integrated_df.to_csv(output_file, sep='\t', index=False)
        #print(f"Integrated DataFrame saved to {output_file}")

        # Write the integrated DataFrame to a VCF file
        write_vcf_from_dataframe(integrated_df, vcf_output_file)
        #print(f"VCF file saved to {vcf_output_file}")
    else:
        write_vcf_from_dataframe(variants_df, vcf_output_file)


