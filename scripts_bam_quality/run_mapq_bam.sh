bam_file=$1
base_name=$(basename "$bam_file" .bam)
output_file=$2

# Check if the BAM file exists
if [ ! -f "$bam_file" ]; then
    echo "Error: BAM file $bam_file does not exist."
    exit 1
fi

echo "Processing BAM file: $bam_file"

samtools view "$bam_file" | cut -f5 | sort | uniq -c | awk -v sample="$base_name" '{print sample "\t" $2 "\t" $1}' > "$output_file" 
