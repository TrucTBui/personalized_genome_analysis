import os
import subprocess
import matplotlib.pyplot as plt
import pandas as pd

def get_statistics_from_vcf(vcf_dir):
    f0000 = os.path.join(vcf_dir, "0000.vcf.gz")
    f0002 = os.path.join(vcf_dir, "0002.vcf.gz")

    cmd1 = f"/home/b/buit/miniconda3/envs/HiWi/bin/bcftools query -f '%CHROM\\t%POS\\n' {f0000} | sort -u | wc -l"
    cmd2 = f"/home/b/buit/miniconda3/envs/HiWi/bin/bcftools query -f '%CHROM\\t%POS\\n' {f0002} | sort -u | wc -l"

    private_stat = subprocess.check_output(cmd1, shell=True).strip()
    shared_stat = subprocess.check_output(cmd2, shell=True).strip()
    private_stat = int(private_stat)
    shared_stat = int(shared_stat)

    return private_stat, shared_stat

def get_stats_family():
    FAMILY = [
    "grandfather_father", "grandmother_father",  # paternal grandparents
    "father", "child", "mother", "aunt",         # nuclear family
    "grandmother_mother", "grandfather_mother"   # maternal grandparents
    ]
    stats_map = {}
    for person in FAMILY:
        vcf_dir_56 = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_VCF/{person}/comparison_56/snps"
        stats_map[f"{person}_56"] = get_statistics_from_vcf(vcf_dir_56)
        if person not in ["aunt", "grandmother_father"]:
            vcf_dir_TSA = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_VCF/{person}/comparison_TSA/snps"
            stats_map[f"{person}_TSA"] = get_statistics_from_vcf(vcf_dir_TSA)
    return stats_map

def create_plot_data(stats_map):
    plot_data = []
    for person, (private, shared) in stats_map.items():
        plot_data.append({
            "Person": person,
            "Newly Identified Variants": private,
            "Shared Variants": shared
        })
    return plot_data

def plot_statistics(plot_data, output_path):
    df = pd.DataFrame(plot_data)
    df.set_index("Person", inplace=True)
    ax = df.plot(kind='bar', stacked=True, figsize=(10, 6))
    # Add values on top of the bars
    for p in ax.patches:
        ax.annotate(f'{int(p.get_height())}', (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='top', fontsize=10, color='black')
        
    #for container in ax.containers:
    #    ax.bar_label(container, label_type='center', fmt='{:,.0f}', color='black')

    ax.set_ylabel('Number of Variants')
    ax.set_title('VCF Comparison Statistics by Person (SNPs)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

if __name__ == "__main__":
    stats_map = get_stats_family()
    plot_data = create_plot_data(stats_map)
    output_path = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_VCF/variant_statistics_snps.png"
    plot_statistics(plot_data, output_path)
    print(f"Statistics plot saved to {output_path}")
    