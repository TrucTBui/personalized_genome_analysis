"""
Whatshap usage:
whatshap phase \
    -o phase_vcf (gz) \
    --reference ref_fasta \
    unphased_vcf (gz) \
    bam \
"""

import pysam
import os
import tempfile
import subprocess
import argparse

def change_sample_name_with_bcftools(input_vcf, output_vcf, new_sample_name):
    
    bcftools_path = "/home/b/buit/miniconda3/envs/HiWi/bin/bcftools"
    samples_file = ""
    try:
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".txt") as f:
            f.write(new_sample_name + '\n')
            samples_file = f.name

        # Construct the command
        reheader_args = [
            bcftools_path,
            "reheader",
            "--samples", samples_file,
            "-o", output_vcf,
            input_vcf
        ]
        
        # Run the command
        subprocess.run(reheader_args, check=True, capture_output=True, text=True)
        print(f"Successfully created new VCF: '{output_vcf}'")

    except FileNotFoundError:
        print(f"Error: '{bcftools_path}' not found. Please ensure bcftools is installed and in your PATH.")
    except subprocess.CalledProcessError as e:
        print("Error during bcftools execution:")
        print(' '.join(str(arg) for arg in e.args))
        print("Stderr:", e.stderr)
    finally:
        # Clean up the temporary samples file
        if os.path.exists(samples_file):
            os.remove(samples_file)

        # tabix the output VCF file
        if os.path.exists(output_vcf):
            try:
                subprocess.run(["/home/b/buit/miniconda3/envs/HiWi/bin/tabix",output_vcf], check=True)
                print(f"Tabix indexing completed for '{output_vcf}'.")
            except subprocess.CalledProcessError as e:
                print("Error during tabix indexing:")
                print(' '.join(str(arg) for arg in e.args))
                print("Stderr:", e.stderr)

def get_sample_name_from_bam_file(bam_file):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        return bam.header['RG'][0]['SM']  

def phase_haplotype_individual(output_vcf, unphased_vcf, bam):
    whatshap_path = "/home/b/buit/miniconda3/envs/HiWi/bin/whatshap"
    reference_fasta="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genomes/human_g1k_v37.fasta"
    
    try:
        sample_name = get_sample_name_from_bam_file(bam)
        unphased_vcf_renamed = unphased_vcf.replace(".vcf.gz", f"_unphased.vcf.gz")
        change_sample_name_with_bcftools(unphased_vcf, unphased_vcf_renamed, sample_name)

        cmd = [whatshap_path, "phase",
           "-o", output_vcf,
           "--reference", reference_fasta,
           unphased_vcf_renamed,
           bam]
        
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"Phasing completed successfully for sample {sample_name}. Output VCF: '{output_vcf}'")
    except FileNotFoundError:
        print(f"Error: '{whatshap_path}' not found. Please ensure whatshap is installed and in your PATH.")
    except subprocess.CalledProcessError as e:
        print("Error during whatshap execution:")
        print(' '.join(str(arg) for arg in e.args))
        print("Stderr:", e.stderr)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phase haplotypes for an individual using Whatshap.")
    parser.add_argument("-o", "--output_vcf", required=True, help="Output phased VCF file (gzipped).")
    parser.add_argument("-v", "--unphased_vcf", required=True, help="Input unphased VCF file (gzipped).")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file for the individual.")

    args = parser.parse_args()

    phase_haplotype_individual(args.output_vcf, args.unphased_vcf, args.bam)