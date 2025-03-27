#!/bin/bash

# Input: BAM files, reference genome with flags
# Output: Individual VCFs for each BAM and a merged VCF (if multiple BAMs)
#         in a specified output directory.

# --- Default Configuration ---
reference_genome=""
bam_string="" # Store BAM files as a space-separated string
region=""
output_dir=""

# --- Helper function to print usage ---
usage() {
  echo "Usage: $0 -f <reference.fasta.gz> -r <chr:start-end> -b \"<bam_file1> [<bam_file2> ...]\" -o <output_directory>"
  echo "Options:"
  echo "  -f <reference.fasta.gz>  Path to the reference genome FASTA file (required)."
  echo "  -r <chr:start-end>       Genomic region to analyze (required)."
  echo "  -b \"<bam_file1> [<bam_file2> ...]\" One or more BAM files to process (required, space-separated within quotes)."
  echo "  -o <output_directory>    Directory to store output files (required)."
  echo "  -i <person_id>  Person id in the family (personal uses)"
  exit 1
}

# --- Parse command-line arguments ---
while getopts "f:r:b:o:i:" opt; do
  case "$opt" in
    f)
      reference_genome="$OPTARG"
      ;;
    r)
      region="$OPTARG"
      ;;
    b)
      bam_string="$OPTARG" # Store all BAM files in a single string
      ;;
    o)
      output_dir="$OPTARG"
      ;;
    i)
      id="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
  esac
done
shift $((OPTIND - 1))

# --- Check for required arguments ---
if [ -z "$reference_genome" ] || [ -z "$region" ] || [ -z "$bam_string" ] || [ -z "$output_dir" ]; then
  echo "Error: Missing required arguments." >&2
  usage
fi

# --- Split the BAM file string into an array ---
IFS=' ' read -r -a bam_files <<< "$bam_string"

# --- Create the output directory if it doesn't exist ---
mkdir -p "$output_dir"
output_unmerged_dir="${output_dir}/unmerged_vcf"
mkdir -p "${output_unmerged_dir}"
output_merged_dir="${output_dir}/merged_vcf"
mkdir -p "${output_merged_dir}"
output_bam_dir="${output_dir}/modified_bam"
mkdir -p "${output_bam_dir}"

# --- Function to subset BAM file based on region ---
subset_bam() {
  local bam_file="$1"
  local bam_base=$(basename "$bam_file" .bam)
  local subset_bam="${output_bam_dir}/${bam_base}_ss.bam"

  samtools view -b "$bam_file" "$region" > "$subset_bam"
  samtools index "$subset_bam"

  echo "$subset_bam"
}

# --- Function to process a single BAM file for VCF ---
process_single_bam() {
  local bam_file="$1"
  local bam_base=$(basename "$bam_file" .bam)
  local output_vcf="${output_unmerged_dir}/${bam_base}.vcf.gz"
  local index_file="$bam_file.bai"
  local temp_all_variants="${output_unmerged_dir}/tmp_${bam_base}.vcf.gz"

  echo "Processing single BAM file: $bam_file"

  # Index the BAM file (index will be in the same dir as BAM)
  if [ ! -f "$index_file" ]; then
    echo "Indexing: $bam_file.bai"
    samtools index "$bam_file"
    if [ $? -ne 0 ]; then
      echo "Error indexing: $bam_file"
      return 1
    fi
  fi

  # Call variants
  echo "Calling variants into: $output_vcf"
  bcftools mpileup -Ou -f "$reference_genome" -r "$region" -a AD,DP "$bam_file" | \
  bcftools call -mv -Oz -o "$temp_all_variants" -
  bcftools view -v snps -Oz -o "$output_vcf" "$temp_all_variants"

  if [ $? -ne 0 ]; then
    echo "Error calling variants for: $bam_file"
    return 1
    else
      rm "$temp_all_variants"
  fi

  echo "Successfully processed: $bam_file -> $output_vcf"
  echo "Output VCF file: $output_vcf"
  return 0
}

# --- Subset all BAM files ---
declare -a subset_bam_files
for bam_file in "${bam_files[@]}"; do
  subset_bam_files+=("$(subset_bam "$bam_file")")
done
bam_files=("${subset_bam_files[@]}")

# --- Process based on the number of BAM files ---
if [ ${#bam_files[@]} -ge 1 ]; then
  # Process each BAM file individually
  echo "Generating individual VCF files..."
  for bam_file in "${bam_files[@]}"; do
    process_single_bam "$bam_file"
  done

  echo "Processing multiple BAM files - merging..."
  declare -a no_rg_bams
  declare -a bam_basenames

  # --- Process each input BAM file to remove @RG lines and store basenames ---
  echo "Removing @RG lines from input BAM files..."
  for bam_file in "${bam_files[@]}"; do
    base_name=$(basename "$bam_file" .bam)
    bam_basenames+=("$base_name")
    no_rg_bam="${output_bam_dir}/${base_name}_no_rg.bam"

    if [ ! -f "$no_rg_bam" ]; then
      samtools view -h "$bam_file" | grep -v '@RG' | samtools view -b -o "$no_rg_bam"
      if [ ! $? -eq 0 ]; then
        echo "Error removing @RG lines from: $bam_file"
        exit 1
      fi
    fi

    no_rg_bams+=("$no_rg_bam")

  done

  # --- Construct the merged filename ---
  merged_basename="merged_${bam_basenames[0]}"
  for i in "${!bam_basenames[@]}"; do
    if [[ "$i" -gt 0 ]]; then
      merged_basename="${merged_basename}_${bam_basenames[$i]}"
    fi
  done
  merged_bam="${output_bam_dir}/${merged_basename}.bam"
  merged_vcf="${output_merged_dir}/${id}_${merged_basename}.vcf.gz"
  temp_all_variants="${output_merged_dir}/tmp_${id}_${merged_basename}.vcf.gz"


  # --- Merge the no-RG BAM files ---
  echo "Merging the no-RG BAM files into: $merged_bam"
  samtools merge -f "$merged_bam" "${no_rg_bams[@]}"
  if [ $? -eq 0 ]; then
    echo "Successfully merged BAM files."
  else
    echo "Error merging BAM files."
    exit 1
  fi

  # --- Index the merged BAM file ---
  echo "Indexing the merged BAM file: $merged_bam.bai"
  samtools index "$merged_bam"
  if [ $? -eq 0 ]; then
    echo "Successfully indexed merged BAM file."
  else
    echo "Error indexing merged BAM file."
    exit 1
  fi


  # --- Call variants on the merged BAM ---
  echo "Calling variants on the merged BAM into: $merged_vcf"

  bcftools mpileup -Ou -f "$reference_genome" -r "$region" -a AD,DP "$merged_bam" | \
  bcftools call -mv -Oz -o "$temp_all_variants" -
  bcftools view -v snps -Oz -o "$merged_vcf" "$temp_all_variants"


  if [ $? -eq 0 ]; then
    echo "Successfully called variants on merged BAM."
    rm "$temp_all_variants"
  else
    echo "Error calling variants on merged BAM."
    exit 1
  fi

  echo "Merged VCF file generated in: $merged_vcf"
fi

exit 0