reference_genome="/mnt/raidinput/input/own/ReferenceGenomes/human_g1k_v37.fasta.gz"
region="17:7565096-7590855"
output_dir="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/VCF/playground"

# Define BAM file pairs for each run
bam_files=(
  "/mnt/raidinput/input/own/Genomes/RZ2/22-FA1811137/TSAB2086.bam /mnt/raidinput/input/own/Genomes/RZ/4412_FA1811137/56001811224412.bam"
  "/mnt/raidinput/input/own/Genomes/RZ/4164_JA2012040/56001811224164.bam"
  "/mnt/raidinput/input/own/Genomes/RZ2/SA3018361/TSAB1838.bam /mnt/raidinput/input/own/Genomes/RZ/0324_SA3018361/56001808050324.bam"
  "/mnt/raidinput/input/own/Genomes/RZ2/22-UA010903/TSAB2165.bam /mnt/raidinput/input/own/Genomes/RZ/0381_UA010903/56001808050381.bam"
  "/mnt/raidinput/input/own/Genomes/RZ2/CH3118662/TSAB1874.bam /mnt/raidinput/input/own/Genomes/RZ/0343_CH3118662/56001808050343.bam"
  "/mnt/raidinput/input/own/Genomes/RZ/4448_TX3118662/56001811224448.bam"
  "/mnt/raidinput/input/own/Genomes/RZ2/22-UH1911438/TSAB2138.bam /mnt/raidinput/input/own/Genomes/RZ/4074_UH1911438/56001811224074.bam"
  "/mnt/raidinput/input/own/Genomes/RZ2/IH1810836/TSAB1829.bam /mnt/raidinput/input/own/Genomes/RZ/4452_IH1810836/56001811224452.bam"
)


# Loop through BAM file pairs and execute the make_vcf.sh command
i=1  # Initialize the ID counter
for bam_pair in "${bam_files[@]}"; do
  # Run make_vcf.sh and capture the VCF output path, passing the ID as -i
  bash make_vcf.sh -f "$reference_genome" -r "$region" -b "$bam_pair" -o "$output_dir" -i "$i"
  # Increment the ID for the next person
  ((i++))
done

vcf_dir="${output_dir}/merged_vcf"
merged_vcf="${output_dir}/family_output.vcf.gz"
merged_vcf_snv="${output_dir}/family_output.vcf.snp.gz"


for vcf_file in "$vcf_dir"/*.vcf.gz; do
  # Index the VCF file if it hasn't been indexed already
  #if [ ! -f "${vcf_file}.tbi" ]; then
    echo "Indexing VCF file with tabix: $vcf_file"
    tabix -p vcf "$vcf_file"
    if [ $? -eq 0 ]; then
      echo "Successfully indexed: $vcf_file"
    else
      echo "Error indexing: $vcf_file"
      exit 1
    fi
  #fi
done

echo "Merging all VCF files..."
bcftools merge -o "$merged_vcf" -Oz "$vcf_dir"/*.vcf.gz

if [ $? -eq 0 ]; then
  echo "Successfully merged VCF files into: $merged_vcf"
else
  echo "Error merging VCF files."
  exit 1
fi

tabix -p vcf "$merged_vcf"
if [ $? -eq 1 ]; then
  echo "Error indexing: $merged_vcf"
  exit 1
fi


