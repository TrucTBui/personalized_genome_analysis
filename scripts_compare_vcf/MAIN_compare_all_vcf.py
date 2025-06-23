import subprocess
import os

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

def run_isec_normalize(pre_computed_vcf, vcf_to_compare, output_folder, person):
    bcftools_path = "/home/b/buit/miniconda3/envs/HiWi/bin/bcftools"
    normalize_script_path = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_compare_vcf/normalize_vcf.txt"
    temp_vcf_name = f"/tmp/{person}_normalized.vcf.gz"

    try:
        annotate_args = [
            bcftools_path, "annotate",
            "--rename-chrs", normalize_script_path,
            pre_computed_vcf,
            "-Oz", 
            "-o", temp_vcf_name 
        ]
        subprocess.run(annotate_args, check=True)

        index_args = [bcftools_path, "index", temp_vcf_name]
        subprocess.run(index_args, check=True)

        isec_args = [
            bcftools_path, "isec",
            "-c", "none",
            "-O", "z",
            "-W",
            "-p", output_folder,
            vcf_to_compare,
            temp_vcf_name 
        ]
        subprocess.run(isec_args, check=True)
    finally:
        if os.path.exists(temp_vcf_name):
            os.remove(temp_vcf_name)
        temp_index_name = temp_vcf_name + ".csi"
        if os.path.exists(temp_index_name):
            os.remove(temp_index_name)


FAMILY = [
    "grandfather_father", "grandmother_father",  # paternal grandparents
    "father", "child", "mother", "aunt",         # nuclear family
    "grandmother_mother", "grandfather_mother"   # maternal grandparents
    ]

vcf_map = {}

for person in FAMILY:
    pre_computed_vcf_details = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genomes/VCF/{person}.txt"
    vcf_list = parse_vcf_path_file(pre_computed_vcf_details)
    vcf_map[person] = vcf_list

# Split VCF files into SNPs and indels for each person in the family
for person in FAMILY:
    full_vcf_path = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_VCF/{person}/all_variants.vcf.gz"
    subprocess.run(["/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_compare_vcf/split_vcf_to_snps_and_indels.sh", full_vcf_path])

# Compare vcfs files with the pre-computed VCF files
for person in FAMILY:
    
    for snp_vcf, indel_vcf in vcf_map[person]:
        # if the basename of the snp_vcf starts with 56:
        if os.path.basename(snp_vcf).startswith("56"):
            output_folder_snp = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_VCF/{person}/comparison_56/snps/"
            output_folder_indel = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_VCF/{person}/comparison_56/indels/"
        else:
            output_folder_snp = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_VCF/{person}/comparison_TSA/snps/"
            output_folder_indel = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_VCF/{person}/comparison_TSA/indels/"
        
        os.makedirs(output_folder_snp, exist_ok=True)
        os.makedirs(output_folder_indel, exist_ok=True)

        # remove the old files in the output folders
        for file in os.listdir(output_folder_snp):
            os.remove(os.path.join(output_folder_snp, file))
        for file in os.listdir(output_folder_indel):
            os.remove(os.path.join(output_folder_indel, file))

        vcf_snp_to_compare = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_VCF/{person}/all_snps.vcf.gz"
        vcf_indel_to_compare = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_VCF/{person}/all_indels.vcf.gz"

        if os.path.basename(snp_vcf).startswith("TSA"):
            #continue
            subprocess.run(["/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_compare_vcf/compare_vcf.sh", vcf_snp_to_compare, snp_vcf, output_folder_snp])
            subprocess.run(["/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_compare_vcf/compare_vcf.sh", vcf_indel_to_compare, indel_vcf, output_folder_indel])
        else:
            run_isec_normalize(snp_vcf, vcf_snp_to_compare, output_folder_snp, person)
            run_isec_normalize(indel_vcf, vcf_indel_to_compare, output_folder_indel, person)



