import subprocess
import argparse
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gzip

parser = argparse.ArgumentParser()

parser.add_argument("-g", "--gene", type=str, required=True, help="which gene to analyse")
parser.add_argument("-c", "--chunk", type=str, required=True, help="which chunk to analyse")

args = parser.parse_args()
gene_id = args.gene
chunk = args.chunk
chrom = chunk.split("_")[0]


def make_line_plot_cov(df, output_line_plot_path):
    """
    Creates a line plot of coverage for each person across genomic positions,
    highlighting positions with coverage < 5.
    """

    # Create the faceted line plot
    g = sns.FacetGrid(df, col="Person", col_wrap=2, height=4, aspect=1.5)
    g.map(sns.lineplot, "Position", "Coverage", errorbar=None)

    # Highlight low coverage positions
    def highlight_low_coverage(data, **kwargs):
        low_cov_data = data[data['Coverage'] < 5]
        if not low_cov_data.empty:
            plt.scatter(low_cov_data['Position'], low_cov_data['Coverage'],
                        color='red', marker='x', s=50, zorder=5, label="Coverage < 5")

    g.map_dataframe(highlight_low_coverage)  # Use map_dataframe instead of map

    # Improve plot aesthetics
    for ax in g.axes.flat:
        ax.grid(False)
        for line in ax.lines:
            line.set_alpha(0.7)

    g.set_titles("{col_name}")

    # Add legend manually
    plt.legend(["Coverage < 5"], loc='upper right')

    g.fig.suptitle(f"Coverage for Gene {gene_id} by Person", y=1.02)
    g.fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    g.savefig(output_line_plot_path)
    plt.close()

def plot_stats_from_multiple_files(stats_files, output_path_gene, gene_name):
    output_path = os.path.join(output_path_gene, "plots")

    if not os.path.exists(output_path):
        os.makedirs(output_path)
        print(f"Directory '{output_path}' created.")

    total_positions = []
    inconsistent_positions = []
    no_counts_positions = []
    haplotype_positions = []
    not_equal_ref = []

    # Read each stats file and extract values
    for person, file in stats_files.items():
        with open(file, "r") as f:
            lines = f.readlines()
            total_positions.append(int(lines[0].split(":")[1].strip()))
            inconsistent_positions.append(int(lines[1].split(":")[1].split("(")[0].strip()))
            no_counts_positions.append(int(lines[2].split(":")[1].split("(")[0].strip()))
            haplotype_positions.append(int(lines[3].split(":")[1].strip()))
            not_equal_ref.append(int(lines[4].split(":")[1].strip()))

    # Create a DataFrame from the stats
    df_stats = pd.DataFrame({
        'Person': list(stats_files.keys()),
        'Total Positions': total_positions,
        'Inconsistent Positions': inconsistent_positions,
        'No Count Positions': no_counts_positions,
        'Haplotype Positions': haplotype_positions,
        'Reference != Final Base': not_equal_ref
    })

    # Bar plot for inconsistent positions
    plt.figure(figsize=(8, 6))
    df_stats.set_index('Person')['Inconsistent Positions'].plot(kind='bar', color='red')
    plt.title(f'Inconsistent Positions - {gene_name}')
    plt.ylabel('Number of Positions')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{output_path}/inconsistent_positions.png")


    # Bar plot for no count positions
    plt.figure(figsize=(8, 6))
    df_stats.set_index('Person')['No Count Positions'].plot(kind='bar', color='blue')
    plt.title(f'No Count Positions - {gene_name}')
    plt.ylabel('Number of Positions')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{output_path}/no_count_positions.png")


    # Bar plot for ref != final base positions
    plt.figure(figsize=(8, 6))
    df_stats.set_index('Person')['Reference != Final Base'].plot(kind='bar', color='orange')
    plt.title(f'Reference != Final Base - {gene_name}')
    plt.ylabel('Number of Positions')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{output_path}/reference.png")


    # Bar plot for haplotype positions
    plt.figure(figsize=(8, 6))
    df_stats.set_index('Person')['Haplotype Positions'].plot(kind='bar', color='green')
    plt.title(f'Haplotype Positions - {gene_name}')
    plt.ylabel('Number of Positions')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{output_path}/haplotype_positions.png")


    plt.close()

def make_combine_results(gene_name, list_results_files):
    df_results_list = []
    for file in list_results_files:
        # Read the file into a DataFrame
        df = pd.read_csv(file, compression='gzip')

        # Add a "Gene" column to keep track of the associated gene
        #df["Gene"] = gene_name

        df = df[['Chromosome', 'Position', 'Gene', 'Type', 'Person', 'Reference_Base',
                 'Coverage', 'Final_Base', 'Consistency', 'Special_Case']]

        df_results_list.append(df)

    if not df_results_list:
        return pd.DataFrame()

    merged_results_df = pd.concat(df_results_list, ignore_index=True)
    return merged_results_df

def transform_results_df(merged_results_df):
    if merged_results_df.empty:
        return pd.DataFrame()
    
    pivoted_df = merged_results_df.pivot_table(
        index=["Chromosome", "Position", "Gene", "Type", "Reference_Base"],
        columns="Person",
        values="Final_Base",
        aggfunc="first"  # Choose the first occurrence if there are duplicates
    )

    pivoted_df = pivoted_df.sort_values(by=["Gene", "Position"], ascending=[True, True])
    pivoted_df.reset_index(inplace=True)

    return pivoted_df

def find_diff_in_results(result_df):
    if result_df.empty:
        return pd.DataFrame()

    differences=[]
    final_base_cols = [col for col in result_df.columns if col not in ["Chromosome", "Position", "Gene", "Type", "Reference_Base"]]

    for index, row in result_df.iterrows():
        ref = row.Reference_Base  # Reference base

        alts = [row[col] for col in final_base_cols]  # Collect all bases from the family


        if len(set(alts)) > 1 or any('/' in alt for alt in alts) or alts[0] is not ref:  # Diff in read
           differences.append(row)
    if len(differences) == 0:
        return pd.DataFrame()

    diff_df = pd.DataFrame(differences)

    return diff_df

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

            if child_combination not in valid_combinations:
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

def gzip_tsv_file(tsv_filepath, gz_filename):
    """
    Compresses a TSV file using gzip.

    Args:
        tsv_filepath (str): The full path to the TSV file.
        gz_filename (str): The desired name for the output gz file (e.g., 'output.tsv.gz').
    """
    try:
        with open(tsv_filepath, 'rb') as tsv_file:
            with gzip.open(gz_filename, 'wb') as gz_file:
                gz_file.writelines(tsv_file)
        print(f"Successfully gzipped '{tsv_filepath}' to '{gz_filename}'")
    except FileNotFoundError:
        print(f"Error: TSV file not found at '{tsv_filepath}'")
    except Exception as e:
        print(f"An error occurred during gzipping: {e}")

def run_analysis(gene_ID:str, base_output_folder:str):

    family_relationships = [["father", "mother", "child"], ["grandfather_father", "grandmother_father", "father"],
                            ["grandfather_mother", "grandmother_mother", "mother"],
                            ["grandfather_mother", "grandmother_mother", "aunt"]]
    

    list_results_files = []
    results_folder = os.path.join(base_output_folder, gene_ID)

    # Find all "*.results_merged.tsv" files inside the gene folder (recursively)
    result_files = glob.glob(os.path.join(results_folder, "**", "results.tsv.gz"), recursive=True)
    for result_file in result_files:
        list_results_files.append(result_file)

    merged_results_df = make_combine_results(gene_ID, list_results_files)
    ## Make a line plot of coverage for each person across genomic positions
    #make_line_plot_cov(merged_results_df, os.path.join(results_folder, "coverage_plot.png"))

    transformed_results_df = transform_results_df(merged_results_df)
    transformed_results_df.to_csv(os.path.join(results_folder, "results_all.tsv.gz"),
                               sep="\t", index=False, compression="gzip")
    
    diff_df = find_diff_in_results(transformed_results_df)
    if not diff_df.empty:
        diff_df.to_csv(os.path.join(results_folder, "merged_results_filtered.tsv"), sep="\t", index=False)
    else:
        print(f"No differences with the reference found for gene {gene_ID}.")
        with open(os.path.join(results_folder, "merged_results_filtered.tsv"), "w") as f:
            f.write("#No differences with the reference found.")
    
    #trouble_pos = find_mendelian_troublesome_positions(diff_df, family_relationships)
    #if not trouble_pos.empty:
    #    trouble_pos.to_csv(os.path.join(base_output_folder, gene_name, "mendelian_troublesome_results.tsv"), sep="\t",
    #                       index=False)
    #else:
    #    with open(os.path.join(base_output_folder, gene_name, "mendelian_troublesome_results.tsv"), "w") as f:
    #        f.write("#No troublesome results found")



run_analysis(gene_id, base_output_folder = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results/{chrom}")