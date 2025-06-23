import os
import csv
import subprocess
import tempfile
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Convert pangenomic analysis to BED format and liftOver.")
parser.add_argument("-i", "--input_file", type=str, required=True, help="Input pangenomic analysis file.")
#parser.add_argument("-o", "--output_file", type=str, required=True, help="Output lifted BED file.")
args = parser.parse_args()
input_file = args.input_file
#output_lifted_file = args.output_file

def liftover_pangenomic_analysis(input_file, output_lifted_file, chrom, chain_file, liftOver_executable):
    # Create a temporary BED file
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".bed") as temp_bed:
        bed_path = temp_bed.name
        writer = temp_bed
        writer.write("#chrom\tstart\tend\tvariant_type\tregion_type\tallele_change\thaplotypes\n")

        with open(input_file, newline='') as infile:
            reader = csv.DictReader(infile, delimiter='\t')
            for row in reader:
                genomic_start = int(row["GenomicStart"])
                variant_type = row["DifferenceType"]
                diff_info = row["DifferenceInfo"]
                haplotypes = row["Haplotypes"]
                n_haplotypes = haplotypes.count("|") + 1
                region_type = row["GenomicRegionType"]

                parts = diff_info.split("|")
                pos = int(parts[0])
                bed_start = genomic_start + pos
                bed_end = bed_start + 1
                allele_change = ""

                # Types: BaseMutation, Deletion, Insertion, SNP, StructuralVariant

                if variant_type.lower() == "deletion":
                    deletion_length = int(parts[1])
                    bed_end = bed_start + deletion_length
                    allele_change = str(deletion_length)
                elif variant_type.lower() == "snp":
                    allele_change = parts[-1]
                elif variant_type.lower() == "insertion":
                    insertion_length = int(parts[3])
                    allele_change = str(insertion_length)
                else: # BaseMutation, StructuralVariant
                    mutation_length = int(parts[1])
                    allele_change = str(mutation_length)
                    bed_end = bed_start + mutation_length


                bed_line = f"{chrom}\t{bed_start}\t{bed_end}\t{variant_type}\t{region_type}\t{allele_change}\t{n_haplotypes}/89\n"
                writer.write(bed_line)

    # Create temporary files for liftOver output
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".unmapped.bed") as temp_unmapped:
        unmapped_path = temp_unmapped.name

    # Run liftOver
    os.makedirs(os.path.dirname(output_lifted_file), exist_ok=True)

    subprocess.run([
        liftOver_executable, bed_path, chain_file, output_lifted_file, unmapped_path, "-bedPlus=3"
    ], check=True)

    # Clean up temp files
    os.remove(bed_path)
    if os.path.exists(unmapped_path):
        os.remove(unmapped_path)

    print(f"Lifted file written to: {output_lifted_file}")

if __name__ == "__main__":
    # extract chromosome from input file path. 
    # The input file format: /mnt/raidproj/proj/projekte/personalizedmed/pangenome/analysis_schertler_stroebele_ws24/[chrom]/[gene_id]/[gene_id]_pangenome_analysis_grouped.tsv 
    chrom = os.path.basename(os.path.dirname(os.path.dirname(input_file))).replace("om","")    
    chain_file = os.path.expanduser("~/tools/LiftOver/hg38ToHg19.over.chain.gz")
    liftOver_executable = os.path.expanduser("~/tools/LiftOver/liftOver")
    output_folder = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Pangenome_analysis/{chrom}/"
    #check if the output folder exists, if not create it
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # Leave out the version suffix from the input file name
    output_lifted_file = os.path.join(output_folder, os.path.basename(input_file).split(".")[0] + "_pangenome_analysis_lifted.bed")
    print(f"Output lifted file: {output_lifted_file}")

    liftover_pangenomic_analysis(input_file, output_lifted_file, chrom, chain_file, liftOver_executable)
