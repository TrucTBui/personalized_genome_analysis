#!/bin/bash
vcf_to_split=$1
#extract the directory of the vcf file
vcf_dir=$(dirname "$vcf_to_split")

/home/b/buit/miniconda3/envs/HiWi/bin/bcftools view -v snps -Oz -o "$vcf_dir/all_snps.vcf.gz" "$vcf_to_split"
tabix -p vcf "$vcf_dir/all_snps.vcf.gz"

/home/b/buit/miniconda3/envs/HiWi/bin/bcftools view -v indels -Oz -o "$vcf_dir/all_indels.vcf.gz" "$vcf_to_split"
tabix -p vcf "$vcf_dir/all_indels.vcf.gz"