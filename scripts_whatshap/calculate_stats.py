import subprocess
import matplotlib.pyplot as plt

FAMILY = [
    "grandfather_father", "grandmother_father",  # paternal grandparents
    "father", "child", "mother", "aunt",         # nuclear family
    "grandfather_mother", "grandmother_mother"    # maternal grandparents
]

family_sample_map = {"grandfather_father": "TSAB2086", "grandmother_father": "56001811224164",
                         "father": "TSAB1838", "child": "TSAB2165",
                         "mother": "TSAB1874", "aunt": "56001811224448",
                         "grandfather_mother": "TSAB1829", "grandmother_mother": "TSAB2138" }

WHATSHAP_RESULTS_DIR = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_whatshap"

def count_number_of_phased_variants(bgzipped_vcf, sample_ID):
    # exclude chromosome Y
    cmd = f"bcftools view -v snps -s {sample_ID} {bgzipped_vcf} | bcftools query -f '%CHROM\t%POS\t[%GT]\n' - | grep -v '^Y\t' | grep '|' | cut -f1,2 | sort -u | wc -l"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Error counting phased variants: {result.stderr.strip()}")
    
    return int(result.stdout.strip())

def count_number_of_unique_hetero_variants(bgzipped_vcf):
    # exclude chromosome Y
    cmd = f"bcftools view -v snps -H {bgzipped_vcf} | grep -v '^Y\t' | cut -f1,2,10 | grep -E '(0/1|1/0|0\|1|1\|0|1/2|2/1|1\|2|2\|1)' | sort -u | wc -l"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Error counting unique variants: {result.stderr.strip()}")
    return int(result.stdout.strip())

def create_bar_plot(is_read_based, output):
    stats_map = {}
    for person in FAMILY:
        bgzipped_vcf = f"{WHATSHAP_RESULTS_DIR}/{person}/phased_variants.vcf.gz"
        all_variants = count_number_of_unique_hetero_variants(bgzipped_vcf)
        if is_read_based:
            phased_variants = count_number_of_phased_variants(bgzipped_vcf, family_sample_map[person])
        else:
            bgzipped_vcf = f"{WHATSHAP_RESULTS_DIR}/correctly_phased_combined_variants.vcf.gz"
            phased_variants = count_number_of_phased_variants(bgzipped_vcf, family_sample_map[person])
        
        stats = {
            'phased_variants': phased_variants,
            'all_heterozygous_variants': all_variants
        }

        stats_map[person] = stats
        
        print(f"{person}: {stats['phased_variants']} phased variants, {stats['all_heterozygous_variants']} total unique variants")
    # stacked bar plot

    all_variants = [stats_map[p]['all_heterozygous_variants'] for p in FAMILY]
    phased_variants = [stats_map[p]['phased_variants'] for p in FAMILY]

    fig, ax = plt.subplots(figsize=(12, 7))

    bar1 = ax.bar(FAMILY, all_variants, color='lightgray', label='All Unique Heterozygous Variants')

    bar2 = ax.bar(FAMILY, phased_variants, color='skyblue', label='Phased Variants')

    # Add labels and title
    ax.set_ylabel('Number of Variants')
    if is_read_based:
        ax.set_title('SNPs Phasing Statistics per Family Member (Read-Based Phasing, Chromosme Y Excluded)')
    else:
        ax.set_title('SNPs Phasing Statistics per Family Member (Pedigree Phasing, Chromosme Y Excluded)')
    ax.set_xticks(range(len(FAMILY)))
    ax.set_xticklabels(FAMILY, rotation=45, ha='right')
    ax.legend(loc="center right", bbox_to_anchor=(1.15, 0.7))

    for i, count in enumerate(phased_variants):
        if count > 0: 
           # add count and (percentage):
            percentage = (count / all_variants[i]) * 100 if all_variants[i] > 0 else 0
            ax.text(i, count, f"{count} ({percentage:.1f}%)", ha='center', va='bottom', fontsize=9, color='black')
    for i, count in enumerate(all_variants):
        if count > 0:
            ax.text(i, count, str(count), ha='center', va='top', fontsize=9, color='black')


    plt.tight_layout() 
    plt.savefig(output)


create_bar_plot(True, "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/plots_no_dup/phased_variants_read_based_stats.png")
create_bar_plot(False, "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/plots_no_dup/phased_variants_pedigree_stats_redo.png")
