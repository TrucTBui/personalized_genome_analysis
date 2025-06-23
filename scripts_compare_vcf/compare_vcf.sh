#!/bin/bash
vcf_to_compare=$1
pre_compute_vcf=$2
output_dir=$3

/home/b/buit/miniconda3/envs/HiWi/bin/bcftools isec -c none -O z -W -p $output_dir/ $vcf_to_compare $pre_compute_vcf
#/home/b/buit/miniconda3/envs/HiWi/bin/bcftools isec -c none -p $output_dir/ $vcf_to_compare $pre_compute_vcf

