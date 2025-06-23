#!/bin/bash
set -x

whatshap_path="/home/b/buit/miniconda3/envs/HiWi/bin/whatshap"
reference_fasta="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genomes/human_g1k_v37.fasta"
combined_vcf="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_whatshap/combined_variants.vcf.gz"
output_vcf="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_whatshap/phased_combined_variants.vcf.gz"
ped_file="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_whatshap/family.ped"
log_dir="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_whatshap/logs"

mkdir -p "$log_dir"

BAM_FILES=(
  "/mnt/raidinput/input/own/Genomes/RZ2/22-FA1811137/TSAB2086.bam"
  "/mnt/raidinput/input/own/Genomes/RZ/4164_JA2012040/56001811224164.bam"
  "/mnt/raidinput/input/own/Genomes/RZ2/SA3018361/TSAB1838.bam"
  "/mnt/raidinput/input/own/Genomes/RZ2/22-UA010903/TSAB2165.bam"
  "/mnt/raidinput/input/own/Genomes/RZ2/CH3118662/TSAB1874.bam"
  "/mnt/raidinput/input/own/Genomes/RZ/4448_TX3118662/56001811224448.bam"
  "/mnt/raidinput/input/own/Genomes/RZ2/IH1810836/TSAB1829.bam"
  "/mnt/raidinput/input/own/Genomes/RZ2/22-UH1911438/TSAB2138.bam"
)

"$whatshap_path" phase \
  --reference "${reference_fasta}" \
  -o "${output_vcf}" \
  --ped "${ped_file}" \
  "${combined_vcf}" \
  "${BAM_FILES[@]}" 