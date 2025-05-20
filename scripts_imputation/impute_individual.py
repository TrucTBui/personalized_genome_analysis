import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description="Impute genotype for an individual based on family genotypes and pangenomic analysis.")
parser.add_argument("-a","--ambiguous", type=str, required=True, help="Path to the ambiguous positions file.")
parser.add_argument("-p", "--pangenome_dir", type=str, required=True, help="Path to the pangenome directory.")
parser.add_argument("-f","--family", type=str, required=True, help="Path to the family table file.")
args = parser.parse_args()
ambiguous_file = args.ambiguous
pangenome_dir = args.pangenome_dir
family_table_file = args.family
imputed_file = os.path.join(os.path.dirname(ambiguous_file), "imputed_ambiguous_positions.tsv")
variant_file = os.path.join(os.path.dirname(ambiguous_file), "variants.tsv")
# /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results/chr3/ENSG00000007402/aunt/ambiguous_positions.tsv
chrom = ambiguous_file.split("/")[-4]

PEDIGREE = {
    'child': {'parents': ('mother', 'father'), 'children': [], 'siblings': []}, 
    'mother': {'parents': ('grandmother_mother', 'grandfather_mother'), 'children': ['child'], 'siblings': ['aunt']},
    'father': {'parents': ('grandmother_father', 'grandfather_father'), 'children': ['child'], 'siblings': []},
    'aunt': {'parents': ('grandmother_mother', 'grandfather_mother'), 'children': [], 'siblings': ['mother']}, 
    'grandmother_mother': {'parents': (None, None), 'children': ['mother', 'aunt'], 'siblings': []}, 
    'grandfather_mother': {'parents': (None, None), 'children': ['mother', 'aunt'], 'siblings': []}, 
    'grandmother_father': {'parents': (None, None), 'children': ['father'], 'siblings': []}, 
    'grandfather_father': {'parents': (None, None), 'children': ['father'], 'siblings': []}, 
}
ALL_INDIVIDUALS = list(PEDIGREE.keys())
MALES = ["father", "grandfather_father", "grandfather_mother", "child"]


def get_ambiguous_table(ambiguous_positions_file):
    #structure: Chromosome	Position	Type	Person	Reference_Base	Frequency_TSAB2165	Frequency_56001808050381	Coverage	Final_Base	Consistency	Special_Case
    #1	119840376	intron	child	T	NA	NA	0	N	no_counts	Unknown_Base

    #The first row is the header
    ambiguous_positions_df = pd.read_csv(ambiguous_positions_file, sep="\t", header=0)
    return ambiguous_positions_df

def get_variant_table(variant_file):
    # If The file contains #No variant positions found, then return an empty dataframe
    with open(variant_file, 'r') as f:
        first_line = f.readline().strip()
        if first_line == "#No variant positions found":
            return pd.DataFrame()
    
    variants_df = pd.read_csv(variant_file, sep="\t", header=0)
    return variants_df
    
def get_pangenomic_analysis(pangenome_dir, gene_id):
    # the directory contains the pangenomic analysis results for one chromosome. Each file is named with the gene ID.
    # Check if the gene is in the pangenome analysis
    
    gene_file = os.path.join(pangenome_dir, f"{gene_id}_pangenome_analysis_lifted.bed")
    if not os.path.exists(gene_file):
        #print(f"Pangenome analysis for gene {gene_file} does not exist.")
        return None
    
    # Load the bed file
    # Structure: chromosome	start	end	variant_type	region_type	allel_change	number_of_haplotypes
    pangenome_df = pd.read_csv(gene_file, sep="\t", header=None)
    pangenome_df.columns = ["chromosome", "start", "end", "variant_type", "region_type", "allele_change", "number_of_haplotypes"]
    return pangenome_df

def get_family_table(family_table_file):
    # Load the family table
    # Structure: Chromosome	Position	Gene	Type	Reference_Base	aunt	child	father	grandfather_father	grandfather_mother	grandmother_father	grandmother_mother	mother
    
    #if family file contain "#No differences with the reference found", then return an empty dataframe
    with open(family_table_file, 'r') as f:
        first_line = f.readline().strip()
        #if first_line == "#No differences with the reference found" or first_line empty:
        if first_line.startswith("#No differences with the reference found") or first_line == "":
            return pd.DataFrame()

    family_table_df = pd.read_csv(family_table_file, sep="\t")
    return family_table_df  

def impute_genotype_position(ambiguous_positions_df, pangenome_dir, family_genotype_table, pedigree=PEDIGREE):
    # Confidence in the family: parents > siblings > children

    resolved_ambiguous_positions_df = ambiguous_positions_df.copy()
    ambiguous_positions_df["Imputation"] = "-"
    ambiguous_positions = sorted(ambiguous_positions_df["Position"].unique())
    gene_id = ambiguous_positions_df["Gene"].values[0]
    pangenome_df = get_pangenomic_analysis(pangenome_dir, gene_id)
    
    
    # Iterate over the ambiguous positions
    for position in ambiguous_positions:
        if pangenome_df is not None:
            check_pangenome = True
        else:
            check_pangenome = False
        check_ref = True

        possible_genotypes_with_score = {}
        possible_genotypes_with_score["D"] = {'score': 0.0, 'method':set()}
        # Extract the ambiguos position information
        ambiguous_position = ambiguous_positions_df.loc[ambiguous_positions_df["Position"] == position]
        ref = ambiguous_position["Reference_Base"].values[0]
        geno_ref = f"{ref}/{ref}"
        individual = ambiguous_position["Person"].values[0]
        consistency_check = ambiguous_position["Consistency"].values[0]
        sex_special_cases = individual in MALES and (chrom =="chrX" or chrom == "chrY")
        
        # Look up family genotypes 
        if not family_genotype_table.empty:
            family_current_pos_df = family_genotype_table.loc[family_genotype_table["Position"] == position]
            if not family_current_pos_df.empty:  

                family_current_pos_df = family_current_pos_df.iloc[0]
                
                # Get the pedigree of the individual
                individual_pedigree = pedigree[individual]
                parents = individual_pedigree["parents"]
                children = individual_pedigree["children"]
                siblings = individual_pedigree["siblings"]
                

                if len(siblings) > 0:
                    check_siblings = True
                else:
                    check_siblings = False
                if len(children) > 0:
                    check_children = True
                else:
                    check_children = False 
                

                if sex_special_cases:
                    if chrom == "chrY":
                        # If the chromosome is Y, then only check the father
                        if parents[1] is not None:
                            father_genotype = family_current_pos_df[parents[1]]
                            if father_genotype not in ["N","D"]:
                                check_ref = False
                                p1_alleles = sorted(father_genotype.split('/'))  
                                if len(p1_alleles) == 1:  # Otherwise something is wrong since there is only one copy of Y
                                    # The "/" is only temporary, it will be removed later in the infer_inputed_genotype function
                                    possible_genotypes_with_score[f"{p1_alleles[0]}/{p1_alleles[0]}"] = {'score':1.0, 'method':{"family"}}  

                    elif chrom == "chrX":
                        # If the chromosome is X, then check the mother, since the mother is the one that can pass the X chromosome
                        if parents[0] is not None:
                            mother_genotype = family_current_pos_df[parents[0]]
                            if mother_genotype not in ["N","D"]:
                                check_ref = False
                                p1_alleles = sorted(mother_genotype.split('/'))
                                # The "/" is only temporary, it will be removed later in the infer_inputed_genotype function  
                                if len(p1_alleles) == 1:
                                    possible_genotypes_with_score[f"{p1_alleles[0]}/{p1_alleles[0]}"] = {'score':1.0, 'method':{"family"}}
                                else: # 50% of each allele
                                    possible_genotypes_with_score[f"{p1_alleles[1]}/{p1_alleles[1]}"] = {'score':0.5, 'method':{"family"}}
                                    possible_genotypes_with_score[f"{p1_alleles[0]}/{p1_alleles[0]}"] = {'score':0.5, 'method':{"family"}}
                else:

                    # Check if the individual has parents
                    if parents[0] is not None and parents[1] is not None:
                        father_genotype = family_current_pos_df[parents[1]]
                        mother_genotype = family_current_pos_df[parents[0]]
                        
                        # Check if both the parents have genotypes
                        if not (father_genotype in ["N","D"] or mother_genotype in ["N","D"]):  # Best case, confidence factor = 1.0 
                            check_children = False
                            check_siblings = False
                            check_ref = False

                            factor = 1.0

                            p1_alleles = sorted(father_genotype.split('/'))  
                            p2_alleles = sorted(mother_genotype.split('/'))

                            # If both parents are homozygous (no "/"), the genotype of the child is unique. Score 1.0
                            if len(p1_alleles) == 1 and len(p2_alleles) == 1:
                                if p1_alleles[0] == p2_alleles[0]:  # Homozygous
                                    possible_genotypes_with_score[f"{p1_alleles[0]}/{p1_alleles[0]}"] = {'score':1.0 * factor, 'method':{"family"}}
                                else: # Heterozygous
                                    # if p1 is ref and p2 not, then p1/p2, otherwise p2/p1
                                    if p1_alleles[0] == ref:
                                        possible_genotypes_with_score[f"{p1_alleles[0]}/{p2_alleles[0]}"] = {'score':1.0 * factor, 'method':{"family"}}
                                    else:
                                        possible_genotypes_with_score[f"{p2_alleles[0]}/{p1_alleles[0]}"] = {'score':1.0 * factor, 'method':{"family"}}
                            
                            # If both parents are heterozygous (with "/"), the genotype of the child is not unique
                            elif len(p1_alleles) == 2 and len(p2_alleles) == 2:
                                # If parents have the same alleles pairs, then the probability is 25 50 25
                                # E.g: A/G, A/G -> 0.25 A/A, 0.5 A/G, 0.25 G/G
                                if p1_alleles[0] == p2_alleles[0] and p1_alleles[1] == p2_alleles[1]:  # It works because they are sorted
                                    possible_genotypes_with_score[f"{p1_alleles[0]}/{p2_alleles[0]}"] = {'score': 0.25 * factor, 'method':{"family"}}
                                    possible_genotypes_with_score[f"{p1_alleles[0]}/{p2_alleles[1]}"] = {'score': 0.5 * factor, 'method':{"family"}}
                                    possible_genotypes_with_score[f"{p1_alleles[1]}/{p2_alleles[1]}"] = {'score': 0.25 * factor, 'method':{"family"}}
                                # If parents have different alleles pairs, then the probability is 25 25 25 25
                                # E.g: A/G, C/T -> 0.25 A/C, 0.25 A/T, 0.25 G/C, 0.25 G/T
                                else:
                                    possible_genotypes_with_score[f"{p1_alleles[0]}/{p2_alleles[0]}"] = {'score': 0.25 * factor, 'method':{"family"}}
                                    possible_genotypes_with_score[f"{p1_alleles[0]}/{p2_alleles[1]}"] = {'score': 0.25 * factor, 'method':{"family"}}
                                    possible_genotypes_with_score[f"{p1_alleles[1]}/{p2_alleles[0]}"] = {'score': 0.25 * factor, 'method':{"family"}}
                                    possible_genotypes_with_score[f"{p1_alleles[1]}/{p2_alleles[1]}"] = {'score': 0.25 * factor, 'method':{"family"}}

                            # If one parent is homozygous and the other is heterozygous, the genotype of the child is not unique
                            # E.g: A/A, A/G -> 0.5 A/A, 0.5 A/G
                            else:
                                if len(p1_alleles) == 1:
                                    homozygous_genotype = p1_alleles[0]
                                    heterozygous_genotype = p2_alleles
                                else:
                                    homozygous_genotype = p2_alleles[0]
                                    heterozygous_genotype = p1_alleles

                                possible_genotypes_with_score[f"{homozygous_genotype}/{heterozygous_genotype[0]}"] = {'score': 0.5 * factor, 'method':{"family"}}
                                possible_genotypes_with_score[f"{homozygous_genotype}/{heterozygous_genotype[1]}"] = {'score': 0.5 * factor, 'method':{"family"}}

                        # If both parents are missing (N), then either they experience sequencing error, or there is actually an indel (insertion)
                        # In this case, the indel hypothesis is very probable
                        elif father_genotype == "N" and mother_genotype == "N":
                            if consistency_check == "no_counts": #TODO: review this threshold
                                possible_genotypes_with_score["D"]['score'] += 0.4 
                                possible_genotypes_with_score["D"]['method'].add("family")

                        # Both parents have deletion
                        elif father_genotype == "D" and mother_genotype == "D":
                            possible_genotypes_with_score["D"] = {'score': 1.0, 'method':{"family"}}
                        
                        # If one parent is missing (N), then either that parent experiences sequencing error, or there is actually an indel (deletion)
                        else: 
                            factor = 0.8
                            if father_genotype == "N":
                                parent_genotypes = sorted(mother_genotype.split("/"))    
                            else:
                                parent_genotypes = sorted(father_genotype.split("/"))
                            if len(parent_genotypes) == 1:
                                    possible_genotypes_with_score[f"{parent_genotypes[0]}/{parent_genotypes[0]}"] = {'score': 0.5 * factor, 'method':{"family"}}
                            else:  #TODO: review this threshold
                                possible_genotypes_with_score[f"{parent_genotypes[0]}/{parent_genotypes[1]}"] = {'score': 0.5 * factor, 'method':{"family"}}
                                possible_genotypes_with_score[f"{parent_genotypes[0]}/{parent_genotypes[0]}"] = {'score': 0.25 * factor, 'method':{"family"}}
                                possible_genotypes_with_score[f"{parent_genotypes[1]}/{parent_genotypes[1]}"] = {'score': 0.25 * factor, 'method':{"family"}}
                            
                            # If a parent has N and the child also has N, then there is a probability that the child has deletion
                            if consistency_check == "no_counts": #TODO: review this threshold
                                possible_genotypes_with_score["D"]['score'] += 0.2
                                possible_genotypes_with_score["D"]['method'].add("family")
                    
                    # Check siblings
                    if check_siblings:
                        factor = 0.4
                        for sibling in siblings:
                            sibling_genotype = family_current_pos_df[sibling]
                            if sibling_genotype == "N":
                                continue
                            sibling_alleles = sorted(sibling_genotype.split("/"))
                            if len(sibling_alleles) == 1:
                                genotype = f"{sibling_alleles[0]}/{sibling_alleles[0]}"
                            else:
                                genotype = f"{sibling_alleles[0]}/{sibling_alleles[1]}"
                            if genotype not in possible_genotypes_with_score:
                                possible_genotypes_with_score[genotype] = {'score': factor, 'method':{"family"}}
                            else:
                                possible_genotypes_with_score[genotype]['score'] += factor
                    
                    # Check children
                    if check_children:
                        factor = 0.3
                        for child in children:
                            child_genotype = family_current_pos_df[child]
                            if child_genotype == "N":
                                continue
                            child_alleles = sorted(child_genotype.split("/"))
                            if len(child_alleles) == 1:
                                genotype = f"{child_alleles[0]}/{child_alleles[0]}"
                            else:
                                genotype = f"{child_alleles[0]}/{child_alleles[1]}"
                            if genotype not in possible_genotypes_with_score:
                                possible_genotypes_with_score[genotype] = {'score': factor, 'method':{"family"}}
                            else:
                                possible_genotypes_with_score[genotype]['score'] += factor

        
        if check_ref:
            if geno_ref not in possible_genotypes_with_score:
                possible_genotypes_with_score[geno_ref] = {'score': 0.3, 'method':{"reference"}}
            else:
                possible_genotypes_with_score[geno_ref]['score'] += 0.3
                possible_genotypes_with_score[geno_ref]['method'].add("reference")

        # Check pangenome
        if check_pangenome:
            # Take the subset of the pangenome analysis for the position, start <= pos <= end
            pangenome_current_pos_df = pangenome_df.loc[(pangenome_df["start"] <= position) & (pangenome_df["end"] >= position)]
            # Structure: chromosome	start	end	variant_type	region_type	allele_change	number_of_haplotypes    chr1	3392326	3392327	SNP	CDSINTRON	A->G	1/89

            # Only keep these columns: chromosome, start, end, variant_type, allele_change
            pangenome_current_pos_df = pangenome_current_pos_df[["chromosome", "start", "end", "variant_type", "allele_change"]]
            # Remove duplicates
            pangenome_current_pos_df = pangenome_current_pos_df.drop_duplicates()

            if not pangenome_current_pos_df.empty:
                # Iterate over the pangenome_current_pos_df
                for index, row in pangenome_current_pos_df.iterrows():
                    # Get the variant type
                    variant_type = row["variant_type"]
                    # Get the alleles
                    alleles = row["allele_change"]
                    # Process SNP and Deletion
                    if variant_type.lower() == "snp":
                        factor = 0.15
                        allels = alleles.split("->")
                        alt = allels[1]
                        ref_pangenome = allels[0]
                        sorted_alleles = sorted([ref_pangenome, alt])
                        if ref_pangenome == ref:
                            geno_alt_homo = f"{alt}/{alt}"
                            if not sex_special_cases:
                                geno_mix = f"{sorted_alleles[0]}/{sorted_alleles[1]}"
                            else:
                                geno_mix = geno_alt_homo # Only one allel, cannot be heterozygous

                            if geno_alt_homo not in possible_genotypes_with_score:
                                possible_genotypes_with_score[geno_alt_homo] = {'score': factor, 'method':{"pangenome"}}
                            else:
                                possible_genotypes_with_score[geno_alt_homo]['score'] += factor
                                possible_genotypes_with_score[geno_alt_homo]['method'].add("pangenome")
                            if geno_mix not in possible_genotypes_with_score:
                                possible_genotypes_with_score[geno_mix] = {'score': factor, 'method':{"pangenome"}}
                            else:   
                                possible_genotypes_with_score[geno_mix]['score'] += factor
                                possible_genotypes_with_score[geno_mix]['method'].add("pangenome")
                    elif variant_type.lower() == "deletion" and consistency_check == "no_counts":  
                        # If the ambiguous position also has no counts, then very likely that it is deletion
                        if "pangenome" not in possible_genotypes_with_score["D"]['method']:
                            possible_genotypes_with_score["D"]['score'] += 1.1
                            possible_genotypes_with_score["D"]['method'].add("pangenome")
        
            # If the directly previous position is a deletion, then the current position has a probability of being a deletion too if it has no counts
            prev = position - 1
            if prev in ambiguous_positions and consistency_check == "no_counts" and (resolved_ambiguous_positions_df.loc[resolved_ambiguous_positions_df["Position"] == prev, "Special_Case"] == "Deletion").any() and "pangenome" not in possible_genotypes_with_score["D"]['method']:
                    possible_genotypes_with_score["D"]['score'] += 1.1
                    possible_genotypes_with_score["D"]['method'].add("propagation")

        # Integrate the genotype determination of the individual in case of inconsistency
        if consistency_check == "inconsistent":
            assumed_genotype = ambiguous_position["Final_Base"].values[0]
            if "/" not in assumed_genotype: # Homozygous
                assumed_genotype = f"{assumed_genotype}/{assumed_genotype}"
                if assumed_genotype not in possible_genotypes_with_score:
                    possible_genotypes_with_score[assumed_genotype] = {'score': 0.2, 'method':{"individual"}}
                else:
                        possible_genotypes_with_score[assumed_genotype]['score'] += 0.2
                        possible_genotypes_with_score[assumed_genotype]['method'].add("individual")
            else: # Heterozygous
                if not sex_special_cases:  # Sex special cases cannot be heterozygous
                    if assumed_genotype not in possible_genotypes_with_score:
                        possible_genotypes_with_score[assumed_genotype] = {'score': 0.2, 'method':{"individual"}}
                    else:
                        possible_genotypes_with_score[assumed_genotype]['score'] += 0.2
                        possible_genotypes_with_score[assumed_genotype]['method'].add("individual")
           
                

        inferred_genotype, score, method = infer_inputed_genotype(possible_genotypes_with_score, ref=f"{ref}/{ref}")
        if score > 0:
            # Update the ambiguous positions dataframe. final_base is the inferred genotype, imputed column
            resolved_ambiguous_positions_df.loc[resolved_ambiguous_positions_df["Position"] == position, "Final_Base"] = inferred_genotype
            resolved_ambiguous_positions_df.loc[resolved_ambiguous_positions_df["Position"] == position, "Imputation"] = f"{score}|{';'.join(method)}"
            # new family table with the imputed genotype
            #family_genotype_table.loc[family_genotype_table["Position"] == position, individual] = inferred_genotype

            # Change the "special_case" columns
            if inferred_genotype == "D":
                resolved_ambiguous_positions_df.loc[resolved_ambiguous_positions_df["Position"] == position, "Special_Case"] = "Deletion"
            elif "/" in inferred_genotype:
                resolved_ambiguous_positions_df.loc[resolved_ambiguous_positions_df["Position"] == position, "Special_Case"] = "Heterozygous"
            else:
                if inferred_genotype != ref:
                    resolved_ambiguous_positions_df.loc[resolved_ambiguous_positions_df["Position"] == position, "Special_Case"] = "Homozygous!=Ref"
                else:
                    resolved_ambiguous_positions_df.loc[resolved_ambiguous_positions_df["Position"] == position, "Special_Case"] = "-"
    resolved_ambiguous_positions_df = propagate_positions_before_deletion(resolved_ambiguous_positions_df)

    return resolved_ambiguous_positions_df

def propagate_positions_before_deletion(resolved_ambiguous_positions_df):
    """
    Recursively mark upstream 'no_counts' positions as 'Deletion' in 'Special_Case'
    if followed by a position already marked as 'Deletion'.
    """
    df = resolved_ambiguous_positions_df.copy()
    df = df.sort_values("Position").reset_index(drop=True)

    deletion_indices = df.index[df["Special_Case"] == "Deletion"].tolist()

    for deletion_idx in deletion_indices:
        current_pos = df.at[deletion_idx, "Position"]
        i = deletion_idx - 1
        while i >= 0:
            prev_pos = df.at[i, "Position"]

            # Stop if not consecutive
            if current_pos - prev_pos != 1:
                break

            if df.at[i, "Special_Case"] == "Deletion":
                break  # already a deletion, stop

            if df.at[i, "Consistency"] == "no_counts":
                df.at[i, "Special_Case"] = "Deletion"
                df.at[i, "Final_Base"] = "D"
                df.at[i, "Imputation"] = "1.1|propagation"
                current_pos = prev_pos
                i -= 1
            else:
                break  # no longer no_counts, stop

    return df


def infer_inputed_genotype(possible_genotypes_with_score, ref):
    if len(possible_genotypes_with_score) == 0:
        return "N", 0, "-"
    
    # Sort the genotypes by score in descending order
    sorted_genotypes = sorted(
        possible_genotypes_with_score.items(), 
        key=lambda x: x[1]['score'], 
        reverse=True
    )

    # Get the highest score
    highest_score = sorted_genotypes[0][1]['score']
    
    # Get all genotypes with the highest score
    best_candidates = [item for item in sorted_genotypes if item[1]['score'] == highest_score]
    
    # Prefer genotypes with the ref allele
    best_genotype_item = next((item for item in best_candidates if ref in item[0]), best_candidates[0])

    best_genotype = best_genotype_item[0]
    best_score = best_genotype_item[1]['score']
    method = best_genotype_item[1]['method']

    if "/" in best_genotype:
        alleles = best_genotype.split("/")
        if len(set(alleles)) == 1:  # Homozygous
            return_geno = f"{alleles[0]}"
        else:  # Heterozygous
            return_geno = f"{alleles[0]}/{alleles[1]}"
    else:
        return_geno = best_genotype 
    best_score = round(best_score, 2)   
    
    return return_geno, best_score, method

                        
def process_imputed_positions(imputed_positions_df, variants_df):
    #print (ambiguous_file)

    variants_df = variants_df.drop_duplicates(subset=["Chromosome", "Position"], keep="first")
    #print(f"Variants df shape: {variants_df.shape}")
    # Get the subset of the imputed positions that have "Heterozygous" or "Homozygous!=Ref" in the Special_Case column
    imputed_positions_df_subset = imputed_positions_df.loc[imputed_positions_df["Special_Case"].isin(["Heterozygous", "Homozygous!=Ref", "Deletion"])]
    #print("Imputed subset positions:", imputed_positions_df_subset["Position"].tolist())
    #print("Variants positions:", variants_df["Position"].tolist())

    # Remove from the variants df the positions that are already in the imputed positions df
    variants_df = variants_df.loc[~variants_df["Position"].isin(imputed_positions_df["Position"])]
    #print(f"New variants df shape: {variants_df.shape}")

    new_variants_df = pd.concat([variants_df, imputed_positions_df_subset], ignore_index=True)
    # Sort the new_variants_df by Chromosome and Position
    new_variants_df = new_variants_df.sort_values(by=["Chromosome", "Position"])

    # If there is a Deletion in Special_Case, print the gene, and the person
    if "Deletion" in new_variants_df["Special_Case"].values:
        gene = new_variants_df["Gene"].values[0]
        person = new_variants_df["Person"].values[0]
        print(f"Deletion\t{gene}\t{person}")

    return new_variants_df

    

if __name__ == "__main__":
    
    if not os.path.exists(ambiguous_file):
        exit(0)
    
    # If the file contains #No ambiguous positions found, then also exit
    with open(ambiguous_file, 'r') as f:
        first_line = f.readline().strip()
        if first_line.startswith("#No ambiguous positions found"):
            exit(0)

    ambiguous_positions_df = get_ambiguous_table(ambiguous_file)
    

    family_genotype_table = get_family_table(family_table_file)
    

    imputed_table =  impute_genotype_position(ambiguous_positions_df, pangenome_dir, family_genotype_table)
    imputed_table.fillna("NA", inplace=True)
    imputed_table.to_csv(imputed_file, sep="\t", index=False)
    
    variants_df = get_variant_table(variant_file)
    if not variants_df.empty:
        variants_df['Imputation'] = "-"
        new_variants_df = process_imputed_positions(imputed_table, variants_df)
        #print("\nColumns of new_variants_df after process_imputed_positions:")
        #print(new_variants_df.columns)
        new_variants_df.fillna("NA", inplace=True)
        new_variants_df.to_csv(os.path.join(os.path.dirname(ambiguous_file), "variants.tsv"), sep="\t", index=False)

            





