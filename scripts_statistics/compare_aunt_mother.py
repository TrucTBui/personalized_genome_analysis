#Compare CDS variants between aunt and mother
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

def compare_aunt_mother(aunt_file, mother_file, output_file):
    # Load the aunt and mother data
    aunt_df = pd.read_csv(aunt_file, sep='\t')
    mother_df = pd.read_csv(mother_file, sep='\t')

    # Remove duplicates based on Chromosome, Position, Ref, Alt, and Person
    aunt_df.drop_duplicates(subset=['Chromosome', 'Position', 'Ref', 'Alt', 'Person'], inplace=True)
    mother_df.drop_duplicates(subset=['Chromosome', 'Position', 'Ref', 'Alt', 'Person'], inplace=True)
    # Reset index
    aunt_df.reset_index(drop=True, inplace=True)
    mother_df.reset_index(drop=True, inplace=True)

    # Ensure both DataFrames have the same columns
    if set(aunt_df.columns) != set(mother_df.columns):
        raise ValueError("Aunt and mother files must have the same columns.")

    # Compare the two DataFrames, create a df with differences
    # Gene	Chromosome	Position	Region	Person	Ref	Alt	Type	Genotype	AllelDepth	Imputation
    # compare the position and alt
    comparison_df = pd.merge(aunt_df, mother_df, on=['Gene', 'Chromosome', 'Position', 'Region', 'Ref', 'Alt'], 
                             suffixes=('_aunt', '_mother'), how='outer', indicator=True)
    
    for col in comparison_df.columns:
        if isinstance(comparison_df[col].dtype, pd.CategoricalDtype):
            comparison_df[col] = comparison_df[col].astype(object)
    comparison_df.fillna('-', inplace=True)

    # Filter for differences
    #comparison_df = comparison_df[comparison_df['_merge'] != 'both']

    # drop column person and _merge
    comparison_df.drop(columns=['Person_aunt', 'Person_mother'], inplace=True)

    # Sort by chromosome and position
    comparison_df.sort_values(by=['Chromosome', 'Position'], inplace=True)
    # Reset index
    comparison_df.reset_index(drop=True, inplace=True)

    comparison_df.to_csv(output_file, sep='\t', index=False)

def create_plot_from_df(df, output_file):
    """
    Plots the number of mutual and exclusive variants per chromosome
    using an explicit stacking method to ensure all categories are visible.
    """
    print("Generating plot with corrected stacking logic...")

    # --- Step 1: Categorize each variant ---
    df['Variant_Category'] = df['_merge'].map({
        'both': 'Mutual',
        'left_only': 'Aunt_Only',
        'right_only': 'Mother_Only'
    })

    # --- Step 2: Count variants and print totals (before pivoting) ---
    summary_counts = df['Variant_Category'].value_counts()
    mutual_count = summary_counts.get('Mutual', 0)
    aunt_only_count = summary_counts.get('Aunt_Only', 0)
    mother_only_count = summary_counts.get('Mother_Only', 0)
    print(f"Mutual Variants: {mutual_count}")
    print(f"Aunt Only Variants: {aunt_only_count}")
    print(f"Mother Only Variants: {mother_only_count}")

    # --- Step 3: Pivot the data for plotting ---
    # Group by Chromosome and Category, then unstack the category to create columns
    plot_data = df.groupby(['Chromosome', 'Variant_Category']).size().unstack(fill_value=0)

    # Ensure all three columns exist, even if a category has zero variants
    for cat in ['Mutual', 'Aunt_Only', 'Mother_Only']:
        if cat not in plot_data.columns:
            plot_data[cat] = 0

    # Reorder columns to control the stacking order
    plot_data = plot_data[['Mutual', 'Aunt_Only', 'Mother_Only']]

    # --- Step 4: Sort chromosomes for proper plotting order ---
    plot_data['chr_num'] = plot_data.index.astype(str).str.replace('chr', '')
    plot_data['chr_num'] = plot_data['chr_num'].replace({'X': 23, 'Y': 24, 'M': 25, 'MT': 25}).astype(int)
    plot_data = plot_data.sort_values(by='chr_num')
    plot_data = plot_data.drop(columns='chr_num')

    # --- Step 5: Create the plot with manual stacking ---
    plt.figure(figsize=(16, 9))
    sns.set_style("whitegrid")

    chromosomes = plot_data.index

    # Create the bottom layer (Mutual variants)
    plt.bar(chromosomes, plot_data['Mutual'], label='Mutual', color='darkgray')

    # Create the middle layer (Aunt_Only), starting on top of the Mutual layer
    plt.bar(chromosomes, plot_data['Aunt_Only'], bottom=plot_data['Mutual'],
            label='Aunt_Only', color='skyblue')

    # Create the top layer (Mother_Only), starting on top of Mutual + Aunt_Only
    plt.bar(chromosomes, plot_data['Mother_Only'], bottom=plot_data['Mutual'] + plot_data['Aunt_Only'],
            label='Mother_Only', color='lightcoral')


    plt.title('Comparison of CDS Variants by Chromosome between Aunt and Mother', fontsize=18, fontweight='bold')
    plt.xlabel('Chromosome', fontsize=16)
    plt.ylabel('Number of Variants', fontsize=16)
    plt.xticks(rotation=45)
    plt.legend(title='Variant Source')
    # change size of the xticks and yticks
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()

    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")

# Create plot of all CDS variants, regardless of chromosome
def create_plot_all_variants(df, output_file):
    total_variants = len(df)
    mutual_variants = len(df[df['_merge'] == 'both'])
    aunt_only_variants = len(df[df['_merge'] == 'left_only'])
    mother_only_variants = len(df[df['_merge'] == 'right_only'])
    print(f"Total Variants: {total_variants}")
    print(f"Mutual Variants: {mutual_variants}")
    print(f"Aunt Only Variants: {aunt_only_variants}")
    print(f"Mother Only Variants: {mother_only_variants}")
    # Create a DataFrame for the plot
    plot_data = pd.DataFrame({
        'Category': ['Mutual', 'Aunt_Only', 'Mother_Only'],
        'Count': [mutual_variants, aunt_only_variants, mother_only_variants]
    })
    plt.figure(figsize=(12, 8))
    sns.barplot(x='Category', y='Count', data=plot_data, palette='Set2')
    # Print the total number of variants on each bar
    for index, row in plot_data.iterrows():
        plt.text(index, row['Count'] + 0.5, str(row['Count']), ha='center', va='bottom', fontsize=14)
    plt.title('Comparison of CDS Variants (All Chromosomes) between Aunt and Mother', fontsize=16, fontweight='bold')
    plt.xlabel('Variant Category', fontsize=16)
    plt.ylabel('Number of Variants', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig(output_file)


if __name__ == "__main__":
    aunt_file="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/CDS_variants/aunt/aunt_merged.tsv"
    mother_file="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/CDS_variants/mother/mother_merged.tsv"
    output_file="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/CDS_variants/aunt_mother_comparison.tsv"
    compare_aunt_mother(aunt_file, mother_file, output_file)
    print("Comparison between aunt and mother completed.")
    # Create a plot from the comparison DataFrame
    comparison_df = pd.read_csv(output_file, sep='\t')
    plot_output_file = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/CDS_variants/aunt_mother_comparison_plot.png"
    create_plot_from_df(comparison_df, plot_output_file)
    # Create a plot of all variants, regardless of chromosome
    all_variants_plot_file = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/CDS_variants/aunt_mother_all_variants_plot.png"
    create_plot_all_variants(comparison_df, all_variants_plot_file)
    print("Plot for all variants created.")