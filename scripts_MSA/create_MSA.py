import subprocess
import os
import tempfile

def call_msa_pangenome_script(chromosome, start, end, strand, gfa, output_dir):
    reference_hg38_fasta = "/mnt/raidbio/biosoft/Data/GENOMIC/UNMODIFIED/PAN/homo_sapiens/hg38/dna/Homo_sapiens.GRCh38.dna.toplevel.fa"
    fasta_index = "/mnt/raidbio/biosoft/Data/GENOMIC/UNMODIFIED/PAN/homo_sapiens/hg38/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.fai"
    java_path = "/usr/lib64/jvm/java-21/bin/java"
    jar_path = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_MSA/msaPangenome.jar"

    # call jar file:
    command = [
        java_path, "-jar", jar_path,
        "-chr", str(chromosome),
        "-start", str(start),
        "-end", str(end),
        "-str", f"\"{strand}\"",
        "-fasta", reference_hg38_fasta,
        "-fai", fasta_index,
        "-gfa", gfa,
        "-out", output_dir
    ]
    try:
        subprocess.run(command, capture_output=True, text=True, check=True)
        output_file = os.path.join(output_dir, "msa_extracted.tsv")
        # read the output file and return the content
        content = None
        if os.path.exists(output_file):
            with open(output_file, 'r') as file:
                content = file.read()

        if content:
            return content
    except subprocess.CalledProcessError as e:
        print(f"Error running MSA Pangenome script: {e}")
        print("STDOUT:")
        print(e.stdout)
        print("STDERR:")
        print(e.stderr)
    except FileNotFoundError:
        print(f"Error: Java executable not found at {java_path} or JAR file not found at {jar_path}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def get_subset_vcf(vcf_file, chromosome, start, end, output_file):
    command = [
        "/home/b/buit/miniconda3/envs/HiWi/bin/bcftools", "view",
        "-r", f"{chromosome}:{start}-{end}",
        vcf_file, 
        "-o", output_file,
        "-O", "z"  
    ]
    try:
        subprocess.run(command, text=True, check=True)
        os.system(f"/home/b/buit/miniconda3/envs/HiWi/bin/tabix {output_file}") 
        
    except subprocess.CalledProcessError as e:
        print(f"Error extracting VCF subset: {e}")
        print("STDOUT:")
        print(e.stdout)
        print("STDERR:")
        print(e.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def find_haplotype_sequence(vcf, chromosome, start, end, strand, sample_id=None):

    def find_sample_id(vcf):
        # Extract the sample ID from the VCF header. the vcf is bgzipped.
        # the last entry in the #CHROM line is the sample ID
        try:
            with subprocess.Popen(["/home/b/buit/miniconda3/envs/HiWi/bin/bcftools", "query", "-l", vcf], stdout=subprocess.PIPE, text=True) as proc:
                sample_id = proc.stdout.readline().strip()
                if not sample_id:
                    raise ValueError("No sample ID found in the VCF file.")
                return sample_id
        except subprocess.CalledProcessError as e:
            print(f"Error reading VCF file: {e}")
            print("STDOUT:")
            print(e.stdout)
            print("STDERR:")
            print(e.stderr)
    
    reference_hg19_fasta = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genomes/human_g1k_v37.fasta"
    haplotypes = []
    try:
        if not sample_id:
            sample_id = find_sample_id(vcf)    

        # create a temporary VCF file for the subset using tempfile
        with tempfile.NamedTemporaryFile(delete=False, suffix=".vcf.gz") as temp_vcf:
            subset_vcf = temp_vcf.name
            get_subset_vcf(vcf, chromosome, start, end, subset_vcf)

        for i in range(1,3):
            # create fasta temporary files for each haplotype
            with tempfile.NamedTemporaryFile(delete=True, suffix=".fa") as temp_fasta:
                if i == 1:
                    fasta1 = temp_fasta.name
                else:
                    fasta2 = temp_fasta.name
                command_string = (
                    f"/usr/bin/samtools faidx {reference_hg19_fasta} {chromosome}:{start}-{end} "
                    f"| /home/b/buit/miniconda3/envs/HiWi/bin/bcftools consensus "
                    f"-s {sample_id} " 
                    f"-H {i} {subset_vcf} "
                    f"-o {fasta1 if i == 1 else fasta2}"
                )

                subprocess.run(command_string, shell=True, capture_output=True, text=True, check=True)

                # read the fasta file to get the sequence, add it to the haplotypes list
                with open(fasta1 if i == 1 else fasta2, 'r') as fasta_file:
                    sequence = "".join(line.strip() for line in fasta_file if not line.startswith('>'))
                    if strand == "-":
                        sequence = get_reverse_complement(sequence)
                    haplotypes.append(sequence)

        # delete the temporary VCF file
        os.remove(subset_vcf)
        os.remove(f"{subset_vcf}.tbi")  
        return haplotypes

    except subprocess.CalledProcessError as e:
        print(f"Error extracting haplotype sequence: {e}")
        print("STDOUT:")
        print(e.stdout)
        print("STDERR:")
        print(e.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    

def lift_over_coordinates(chromosomehg19, starthg19, endhg19):
    # Lift over coordinates from hg19 to hg38 using liftOver tool.
    try:
        # create a temporary bed file for the coordinates, both hg19 and hg38
        with tempfile.NamedTemporaryFile(delete=False, suffix=".bed") as temp_bed_hg19:
            bed_hg19 = temp_bed_hg19.name
            with open(bed_hg19, 'w') as f:
                f.write(f"{chromosomehg19}\t{starthg19}\t{endhg19}\n")
        with tempfile.NamedTemporaryFile(delete=False, suffix=".bed") as temp_bed_hg38:
            bed_hg38 = temp_bed_hg38.name
        # run liftOver to change the coordinates from hg19 to hg38
        liftOver_path = "/home/b/buit/tools/LiftOver/liftOver"
        chain_file = "/home/b/buit/tools/LiftOver/hg19ToHg38.over.chain.gz"
        liftOver_command = [liftOver_path, bed_hg19, chain_file, bed_hg38, f"{bed_hg19}_unmapped.bed", "-bedPlus=3"]
        subprocess.run(liftOver_command, check=True, capture_output=True, text=True)

        # read the lifted over coordinates from the bed file
        with open(bed_hg38, 'r') as f:
            # the file should contain one line with the lifted over coordinates
            line = f.readline().strip()
            if line:
                parts = line.split('\t')
                if len(parts) >= 3:
                    chromosomehg38 = parts[0].replace("chr", "")  # remove 'chr' prefix if present
                    starthg38 = int(parts[1])
                    endhg38 = int(parts[2])
                    print(f"Lifted over coordinates in hg38: {chromosomehg38}:{starthg38}-{endhg38}")
                else:
                    print("Error: Lifted over coordinates file is malformed.")
                    return None, None, None
        # delete the temporary bed files
        os.remove(bed_hg19)
        os.remove(bed_hg38)
        os.remove(f"{bed_hg19}_unmapped.bed")

        return chromosomehg38, starthg38, endhg38
    except subprocess.CalledProcessError as e:
        print(f"Error during liftOver execution: {e}")
        print("STDOUT:")
        print(e.stdout)
        print("STDERR:")
        print(e.stderr)
    except FileNotFoundError as e:  
        print(f"Error: {e}. Please ensure the liftOver tool and chain file are correctly specified.")

def find_strand_from_gtf(gene_id, gtf_file = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genes/isar.ensembl-75.gtf.gz"):
    # zcat the file, grep for the gene_id and extract the strand information on the 7th column
    try:
        command = f"zcat {gtf_file} | grep -w '{gene_id}' | awk '{{print $7}}' | head -n 1"
        strand = subprocess.check_output(command, shell=True, text=True).strip()
        if strand and strand in ['+', '-']:
            return strand
        else:
            print(f"No strand information found for gene ID: {gene_id}")
            return None
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        print("STDOUT:")
        print(e.stdout)
        print("STDERR:")
        print(e.stderr)

def find_gene_in_pangenome(gene_id, chromosome, pangenome_dir="/mnt/raidbio/biosoft/projekte/pangenome/mini-gfas"):
    # find if the directory contains any file with the gene_id in the name and .gfa extension
    try:
        pangenome_dir_chromosome = os.path.join(pangenome_dir, f"chrom{str(chromosome)}")
        gfa_files = [f for f in os.listdir(pangenome_dir_chromosome) if f.endswith('.gfa') and gene_id in f]
        if not gfa_files:
            print(f"No GFA files found for gene ID: {gene_id} in directory: {pangenome_dir}")
            return None
        
        # if there are multiple files, return the first one
        gfa_file = gfa_files[0]
        gfa_path = os.path.join(pangenome_dir_chromosome, gfa_file)
        return gfa_path
    except FileNotFoundError as e:
        print(f"Error: {e}. Please ensure the pangenome directory exists and is accessible.")

def get_reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))

def apply_pangenomic_structure(hap_seq, pangenome_template_line):
    # Add tab separated segments to the haplotype sequence based on the pangenomic structure
    def get_pangenome_segment_structure(pangenome_template_line):
        segments = pangenome_template_line.strip().split('\t')
        segment_lengths = [len(segment) for segment in segments]
        return segment_lengths
    
    segments_lengths = get_pangenome_segment_structure(pangenome_template_line)
    hap_seq_segments = []
    start = 0
    for length in segments_lengths:
        end = start + length
        hap_seq_segments.append(hap_seq[start:end])
        start = end
    # Add any remaining sequence after the last segment 
    if start < len(hap_seq):
        hap_seq_segments.append(hap_seq[start:])
    return '\t'.join(hap_seq_segments)
    
def msa_all(gene_id, chromosome, start, end, output_dir):
    has_msa_pangenome = False
    strand = find_strand_from_gtf(gene_id)
    
    # Extract MSA for the gene in the pangenome
    gfa = find_gene_in_pangenome(gene_id, chromosome)
    chromosomehg38, starthg38, endhg38 = lift_over_coordinates(f"chr{chromosome}", start, end)
    if gfa and chromosomehg38 and starthg38 and endhg38:
        msa_apngenome = call_msa_pangenome_script(chromosomehg38, starthg38, endhg38, strand, gfa, output_dir)
        if msa_apngenome:
            has_msa_pangenome = True
            print(f" MSA Pangenome for gene {gene_id} in chromosome {chromosome} from {start} to {end} has been created.")

    # Extract haplotype sequences from VCF for the samples

    FAMILY = [
        "grandfather_father", "grandmother_father",  # paternal grandparents
        "father", "child", "mother", "aunt",         # nuclear family
        "grandmother_mother", "grandfather_mother"   # maternal grandparents
    ]

    base_dir = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_whatshap"
    hap_map = {}
    for person in FAMILY:
        vcf = f"{base_dir}/{person}/phased_variants.vcf.gz"
        # tabix the VCF file if it is not indexed
        if not os.path.exists(f"{vcf}.tbi"):
            try:
                subprocess.run(["/home/b/buit/miniconda3/envs/HiWi/bin/tabix", vcf], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error indexing VCF file for {person}: {e}")
                continue
        hap_list = find_haplotype_sequence(vcf, chromosome, start, end, strand)
        if hap_list:
            hap_map[person] = hap_list
            print(f"Haplotype sequences for {person} extracted successfully.")
        else:
            print(f"Failed to extract haplotype sequences for {person}.")
    
    # Create a concatenated MSA file with the pangenome MSA and haplotype sequences
    output_file = os.path.join(output_dir, f"{gene_id}_{start}_{end}_msa.txt")
    if not has_msa_pangenome:
        print(f"No MSA Pangenome found for gene {gene_id} in chromosome {chromosome} from {start} to {end}.")
    with open(output_file, 'w') as f:
        if has_msa_pangenome:
            first_line_pangenome = msa_apngenome.split('\n')[0]
            f.write(">Pangenome\n")
            f.write(msa_apngenome)
            f.write("\n")     
            # Write haplotype sequences to the MSA file
            for person, haplotypes in hap_map.items():
                #f.write(f">{person}_hap1\n{haplotypes[0]}\n")
                #f.write(f">{person}_hap2\n{haplotypes[1]}\n")
                f.write(f">{person}_hap1\n{apply_pangenomic_structure(haplotypes[0], first_line_pangenome)}\n")
                f.write(f">{person}_hap2\n{apply_pangenomic_structure(haplotypes[1], first_line_pangenome)}\n")
        else:
            for person, haplotypes in hap_map.items():
                f.write(f">{person}_hap1\n{haplotypes[0]}\n")
                f.write(f">{person}_hap2\n{haplotypes[1]}\n")

    print(f"MSA file created at {output_file}")

def msa_all_pedigree(gene_id, chromosome, start, end, output_dir):
    has_msa_pangenome = False
    strand = find_strand_from_gtf(gene_id)
    
    # Extract MSA for the gene in the pangenome
    gfa = find_gene_in_pangenome(gene_id, chromosome)
    chromosomehg38, starthg38, endhg38 = lift_over_coordinates(f"chr{chromosome}", start, end)
    if gfa and chromosomehg38 and starthg38 and endhg38:
        msa_apngenome = call_msa_pangenome_script(chromosomehg38, starthg38, endhg38, strand, gfa, output_dir)
        if msa_apngenome:
            has_msa_pangenome = True
            print(f" MSA Pangenome for gene {gene_id} in chromosome {chromosome} from {start} to {end} has been created.")

    # Extract haplotype sequences from VCF for the samples

    FAMILY = [
        "grandfather_father", "grandmother_father",  # paternal grandparents
        "father", "child", "mother", "aunt",         # nuclear family
        "grandmother_mother", "grandfather_mother"   # maternal grandparents
    ]

    family_sample_map = {"grandfather_father": "TSAB2086", "grandmother_father": "56001811224164",
                         "father": "TSAB1838", "child": "TSAB2165",
                         "mother": "TSAB1874", "aunt": "56001811224448",
                         "grandfather_mother": "TSAB1829", "grandmother_mother": "TSAB2138" }

    combined_phased_vcf = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_whatshap/phased_combined_variants.vcf.gz"
    hap_map = {}
    for person in FAMILY:
        # tabix the VCF file if it is not indexed
        if not os.path.exists(f"{combined_phased_vcf}.tbi"):
            try:
                subprocess.run(["/home/b/buit/miniconda3/envs/HiWi/bin/tabix", combined_phased_vcf], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error indexing VCF file for {person}: {e}")
                continue
        hap_list = find_haplotype_sequence(combined_phased_vcf, chromosome, start, end, strand, family_sample_map[person])
        if hap_list:
            hap_map[person] = hap_list
            print(f"Haplotype sequences for {person} extracted successfully.")
        else:
            print(f"Failed to extract haplotype sequences for {person}.")
    
    # Create a concatenated MSA file with the pangenome MSA and haplotype sequences
    output_file = os.path.join(output_dir, f"{gene_id}_{start}_{end}_msa.txt")
    if not has_msa_pangenome:
        print(f"No MSA Pangenome found for gene {gene_id} in chromosome {chromosome} from {start} to {end}.")
    with open(output_file, 'w') as f:
        if has_msa_pangenome:
            first_line_pangenome = msa_apngenome.split('\n')[0]
            f.write(">Pangenome\n")
            f.write(msa_apngenome)
            f.write("\n")     
            # Write haplotype sequences to the MSA file
            for person, haplotypes in hap_map.items():
                #f.write(f">{person}_hap1\n{haplotypes[0]}\n")
                #f.write(f">{person}_hap2\n{haplotypes[1]}\n")
                f.write(f">{person}_hap1\n{apply_pangenomic_structure(haplotypes[0], first_line_pangenome)}\n")
                f.write(f">{person}_hap2\n{apply_pangenomic_structure(haplotypes[1], first_line_pangenome)}\n")
        else:
            for person, haplotypes in hap_map.items():
                f.write(f">{person}_hap1\n{haplotypes[0]}\n")
                f.write(f">{person}_hap2\n{haplotypes[1]}\n")

    print(f"MSA file created at {output_file}")


# Example usage:

# Using VCF from individually phased samples
#msa_all("ENSG00000001460", 1, 24683724, 24684895, "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/test_pangenom_2")
msa_all("ENSG00000004975", 17, 7129840, 7133162, "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/test_pangenome_out")

# Using VCF from pedigree-phased samples
#msa_all_pedigree("ENSG00000001460", 1, 24683724, 24684895, "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/test_pangenom_2")




 