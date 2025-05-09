import pandas as pd
import os
import argparse
import subprocess
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pandas")

parser = argparse.ArgumentParser(description="Impute genotype for the whole family.")
parser.add_argument("-g","--gene", type=str, required=True, help="Gene to process")
parser.add_argument("-c", "--chunk", type=str, required=True, help="which chunk to analyse")
args = parser.parse_args()

gene = args.gene
chunk = args.chunk
chrom = chunk.split("_")[0]
analysis_folder = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results/{chrom}/"
pangenome_folder = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Pangenome_analysis/{chrom}/"
FAMILY = [
    "grandfather_father", "grandmother_father",  # paternal grandparents
    "father", "child", "mother", "aunt",         # nuclear family
    "grandmother_mother", "grandfather_mother"   # maternal grandparents
]
FAMILY_RELATIONSHIP = [["father", "mother", "child"], ["grandfather_father", "grandmother_father", "father"],
                            ["grandfather_mother", "grandmother_mother", "mother"],
                            ["grandfather_mother", "grandmother_mother", "aunt"]]

family_relationships_males_chromY = [["child","father"], ["father","grandfather_father"]]
family_relationship_males_chromX = [["child","mother"], ["father","grandmother_father"]]
MALES = ["father", "child", "grandfather_father", "grandfather_mother"]

def impute_all():
    family_merged_result_file = os.path.join(analysis_folder, gene, "merged_results_filtered.tsv")
    if not chrom == "chrY":
        for person in FAMILY:
            ambiguous_file = os.path.join(analysis_folder, gene, person, "ambiguous_positions.tsv")
            cmd = f"/home/b/buit/miniconda3/envs/HiWi/bin/python  /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_imputation/impute_individual.py -a {ambiguous_file} -f {family_merged_result_file} -p {pangenome_folder}"
            # If the imputation file exists, skip it
            #imputation_file = os.path.join(analysis_folder, gene, person, "imputed_ambiguous_positions.tsv")
            #if os.path.exists(imputation_file) and os.path.exists(ambiguous_file):
            #    continue
            subprocess.run(cmd, shell=True, check=True)
    else:
        for person in MALES:
            ambiguous_file = os.path.join(analysis_folder, gene, person, "ambiguous_positions.tsv")
            cmd = f"/home/b/buit/miniconda3/envs/HiWi/bin/python  /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_imputation/impute_individual.py -a {ambiguous_file} -f {family_merged_result_file} -p {pangenome_folder}"
            # If the imputation file exists, skip it
            #imputation_file = os.path.join(analysis_folder, gene, person, "imputed_ambiguous_positions.tsv")
            #if os.path.exists(imputation_file) and os.path.exists(ambiguous_file):
            #    continue  
            subprocess.run(cmd, shell=True, check=True)

def merge_variants():
    # Find all files ending with variants.tsv, merge them into a df
    all_files = []
    for person in FAMILY: # TODO: family_chromY
        file_path = os.path.join(analysis_folder, gene, person, "variants.tsv")
        if os.path.exists(file_path):
            all_files.append(file_path)
        else:
            print(f"File not found: {file_path}")
    if not all_files:
        print("No files found to merge.")
        return None
    df_list = []
    for file in all_files:
        # if file contain "#No variant positions found.", skip it
        with open(file, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith("#No variant positions found."):
                continue
        df = pd.read_csv(file, sep="\t")
        df = df[['Gene', 'Chromosome', 'Position', 'Type', 'Person', 'Reference_Base',
                 'Final_Base']]

        df_list.append(df)
    if len(df_list) == 0:
        return None
    merged_df = pd.concat(df_list, ignore_index=True)
    merged_df = merged_df.drop_duplicates()
    merged_df = merged_df.reset_index(drop=True)
    return merged_df

def transform_pivot_variants(merged_df):
    # Create the pivot table
    pivoted_df = merged_df.pivot_table(
        index=["Gene", "Chromosome", "Position", "Type", "Reference_Base"],
        columns="Person",
        values="Final_Base",
        aggfunc="first"
    ).reset_index()

    # Ensure all FAMILY members are present as columns
    for person in FAMILY:
        if person not in pivoted_df.columns:
            pivoted_df[person] = pd.NA  

    # Fill NaNs in FAMILY person columns with the Reference_Base
    base_columns = ["Gene", "Chromosome", "Position", "Type", "Reference_Base"]
    person_cols = [col for col in pivoted_df.columns if col not in base_columns]
    pivoted_df[person_cols] = pivoted_df[person_cols].astype(object)
    
    if not chrom == "chrY":
        for person in FAMILY:
            pivoted_df[person] = pivoted_df[person].fillna(pivoted_df["Reference_Base"])
        final_columns = base_columns + FAMILY
        pivoted_df = pivoted_df[final_columns]
    else:
        for person in MALES:
            pivoted_df[person] = pivoted_df[person].fillna(pivoted_df["Reference_Base"])
        final_columns = base_columns + MALES
        pivoted_df = pivoted_df[final_columns]
    

    return pivoted_df


def find_mendelian_troublesome_positions(df, family_relationships):
    """Find positions that do not make sense genetically amond the family.
    E.g.: Father: G, Mother: G, Child: G/T"""
    if df.empty:
        return pd.DataFrame()

    trouble_positions = []
    continue_outer = False
    result = pd.DataFrame()

    # Iterate through each row
    for _, row in df.iterrows():

        for family in family_relationships:

            father_genotype = row[family[0]]
            mother_genotype = row[family[1]]
            child_genotype = row[family[2]]

            p1_alleles = father_genotype.split('/')  
            p2_alleles = mother_genotype.split('/')
            if "/" in child_genotype:
                child_alleles = child_genotype.split('/')
            else:
                child_alleles = [child_genotype, child_genotype]  # Homozygous case
            if any("N" in allels for allels in [p1_alleles, p2_alleles, child_alleles]):
                continue_outer = True
                break

            # Combine mother and father alleles, since child can only have those alleles
            valid_combinations = {tuple(sorted([fa, ma])) for fa in p1_alleles for ma in p2_alleles}
            child_combination = tuple(sorted(child_alleles))

            if (
                child_combination not in valid_combinations
                and "D" not in child_combination
                and all("D" not in valid_combination for valid_combination in valid_combinations)
            ):
                trouble_positions.append(row)
                continue_outer = True
                break                

        if continue_outer:
            continue_outer = False
            continue

    if len(trouble_positions) > 0:
            result = pd.DataFrame(trouble_positions)
            result = result.drop_duplicates()
    else:
        return pd.DataFrame()


    return result

def find_mendelian_troublesome_positions_males_sex_chrom(df, family_relationships):
    if df.empty:
        return pd.DataFrame()

    trouble_positions = []
    continue_outer = False
    result = pd.DataFrame()

    # Iterate through each row
    for _, row in df.iterrows():

        for family in family_relationships:

            child_genotype = row[family[0]]
            parent_genotype = row[family[1]]
            
            parent_alleles = parent_genotype.split('/')  

            if "/" in child_genotype:
                trouble_positions.append(row)  # Because in males there can be only homozygous genotypes for X and Y.
                continue_outer = True
                break

            if child_genotype == "D":
                continue

            if child_genotype not in parent_alleles:
                trouble_positions.append(row)
                continue_outer = True
                break                       

        if continue_outer:
            continue_outer = False
            continue

    if len(trouble_positions) > 0:
            result = pd.DataFrame(trouble_positions)
            result = result.drop_duplicates()
    else:
        return pd.DataFrame()


    return result



if __name__ == "__main__":
    print(f"Processing gene: {gene}")
    impute_all()
    merged_df = merge_variants()
    if merged_df is None or merged_df.empty:
        print(f"No variants found after imputation for gene {gene}.")
        exit(0)
    pivoted_df = transform_pivot_variants(merged_df)
    
    if not chrom in {"chrY", "chrX"}:
        mendelian_troublesome_positions = find_mendelian_troublesome_positions(pivoted_df, FAMILY_RELATIONSHIP)
        if not mendelian_troublesome_positions.empty:
            mendelian_troublesome_positions.to_csv(os.path.join(analysis_folder, gene, "mendelian_troublesome_positions_after_imputation.tsv"), sep="\t", index=False)
        else:   
            print(f"No mendelian troublesome positions found for gene {gene}.")
    else:
        if chrom == "chrY":
            mendelian_troublesome_positions = find_mendelian_troublesome_positions_males_sex_chrom(pivoted_df, family_relationships_males_chromY)
        else:
            mendelian_troublesome_positions = find_mendelian_troublesome_positions_males_sex_chrom(pivoted_df, family_relationship_males_chromX)

        if not mendelian_troublesome_positions.empty:
            mendelian_troublesome_positions.to_csv(os.path.join(analysis_folder, gene, "mendelian_troublesome_positions_after_imputation.tsv"), sep="\t", index=False)
        else:
            print(f"No mendelian troublesome positions found for gene {gene}.")    
    # Save the pivoted DataFrame to a file
    output_file = os.path.join(analysis_folder, gene, "variants_after_imputation_family.tsv")
    pivoted_df.to_csv(output_file, sep="\t", index=False)
