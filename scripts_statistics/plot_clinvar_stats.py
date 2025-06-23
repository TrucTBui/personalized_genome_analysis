import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def read_clinvar_output(file_path):
    
    try:
        df = pd.read_csv(file_path, header=None, sep='\t')
        # ./scripts_integrate_clinvar/logs/run_chr10_chunk2.out:ENSG00000196177	10	124797364	G	A	Conflicting_classifications_of_pathogenicity(Pathogenic(5)|Likely_pathogenic(5)|Uncertain_significance(1))	grandmother_father
        # turn into a df: GeneID,Chromosome,Position,Ref,Alt,VariantType,Classification,Person
        df.columns = ['GeneID', 'Chromosome', 'Position', 'Ref', 'Alt', 'Classification', 'Person']
        # convert the GeneID: ./scripts_integrate_clinvar/logs/run_chr10_chunk2.out:ENSG00000196177 -> ENSG00000196177
        df['GeneID'] = df['GeneID'].str.split(':').str[1]
        # VariantType: if len(alt) > len(ref) then insertion, if len(alt) < len(ref) then deletion, else snp
        df['VariantType'] = df.apply(lambda row: 'Insertion' if len(row['Alt']) > len(row['Ref']) else ('Deletion' if len(row['Alt']) < len(row['Ref']) else 'SNP'), axis=1)
        # Filter out duplicates: on Chromome, Position, Ref, Alt, Person
        df.drop_duplicates(subset=['Chromosome', 'Position', 'Ref', 'Alt', 'Person'], inplace=True)
        return df[['GeneID', 'Chromosome', 'Position', 'Ref', 'Alt', 'VariantType', 'Classification', 'Person']]
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None

def count_shared_variants(df, output_plot_file):
    #Count the number of shared variants between different persons.
    if df is None:
        return None
    # Group by GeneID, Chromosome, Position, Ref, Alt and count unique persons
    shared_counts = df.groupby(['GeneID', 'Chromosome', 'Position', 'Ref', 'Alt', 'VariantType', 'Classification']).agg({'Person': 'nunique'}).reset_index()
    shared_counts.rename(columns={'Person': 'SharedCount'}, inplace=True)
    # sort by SharedCount in descending order
    shared_counts.sort_values(by='SharedCount', ascending=False, inplace=True)
    shared_counts.reset_index(drop=True, inplace=True)

    # plot the number of shared variants: how many variants are shared by how many persons. x-axis: SharedCount, y-axis: number of variants
     # 1. Prepare the data for the plot.
    # We count how many times each 'SharedCount' value appears.
    # For example, how many variants have a SharedCount of 2, 3, etc.
    plot_data = shared_counts['SharedCount'].value_counts().reset_index()
    plot_data.columns = ['NumberOfPersons', 'VariantCount']
    plot_data.sort_values(by='NumberOfPersons', inplace=True)

    # 2. Create the bar plot.
    plt.figure(figsize=(12, 7))
    sns.set_style("whitegrid")
    
    barplot = sns.barplot(
        x='NumberOfPersons', 
        y='VariantCount', 
        data=plot_data,
        palette='viridis'
    )
    
    # 3. Add labels and titles for clarity.
    plt.title('Distribution of Shared Variants', fontsize=16, fontweight='bold')
    plt.xlabel('Number of Persons Sharing a Variant', fontsize=12)
    plt.ylabel('Count of Unique Variants', fontsize=12)
    
    # Add the count on top of each bar
    for p in barplot.patches:
        barplot.annotate(format(p.get_height(), '.0f'),
                         (p.get_x() + p.get_width() / 2., p.get_height()),
                         ha='center', va='center',
                         xytext=(0, 9),
                         textcoords='offset points')
    # Also add the total number of variants
    total_variants = shared_counts.shape[0]
    plt.figtext(0.99, 0.01, f'Total Variants: {total_variants}', horizontalalignment='right', fontsize=10)

    # 4. Save and show the plot.
    plt.tight_layout()
    plt.savefig(output_plot_file)
    plt.show()
    
    print(f"Plot saved to {output_plot_file}")
    

    return shared_counts

def find_shared_variants(df, branches):
    shared_variants = {}
    for branch in branches:
        # Filter the DataFrame for the current branch
        branch_df = df[df['Person'].isin(branch)]
        # Group by GeneID, Chromosome, Position, Ref, Alt and count unique persons
        shared_counts = branch_df.groupby(['Chromosome', 'Position', 'Ref', 'Alt']).agg({'Person': 'nunique'}).reset_index()
        shared_counts.rename(columns={'Person': 'SharedCount'}, inplace=True)
        # Filter to keep only variants shared by all members of the branch
        shared_counts = shared_counts[shared_counts['SharedCount'] == len(branch)]
        # reset index for better readability
        shared_counts.reset_index(drop=True, inplace=True)
        shared_variants[tuple(branch)] = shared_counts
    return shared_variants
    

df = read_clinvar_output("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/clinvar_out.txt")
#print(df.head())
#print(df.shape)
# Inspect: howmany variants shared by 
branch1 = ['child', 'mother', 'aunt', 'grandmother_mother']
branch2 = ['child', 'mother', 'aunt', 'grandfather_mother']
branch3 = ['child', 'father', 'grandmother_father']
branch4 = ['child', 'father', 'grandfather_father']
branches = [branch1, branch2, branch3, branch4]
# Find out if any of the branches have the same variants
shared_variants = find_shared_variants(df, branches)
for branch, variants in shared_variants.items():
    print(f"Shared variants in branch {branch}:")
    print(variants.head(10))  # Print the first 10 shared variants for each branch

shared_counts = count_shared_variants(df, "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/plots/clinvar_shared_variants_plot.png")

