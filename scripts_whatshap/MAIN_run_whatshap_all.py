import os
import subprocess

def parse_BAM(path, person):
    genomes = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            genomes.append(line.strip())
    # print base names of the genomes
    if person in ["aunt", "grandmother_father"]:
        # take the bam path which has the base name starting with 56
        genomes = [genome for genome in genomes if os.path.basename(genome).startswith("56")]
    else:
        # take the bam path which has the base name starting with TSA
        genomes = [genome for genome in genomes if os.path.basename(genome).startswith("TSA")]

    if len(genomes) != 1:
        raise ValueError(f"Expected exactly one BAM file for {person}, found {len(genomes)}: {genomes}")
    #print(f"Found BAM file for {person}: {genomes[0]}")
    return genomes[0]


FAMILY = [
    "grandfather_father", "grandmother_father",  # paternal grandparents
    "father", "child", "mother", "aunt",         # nuclear family
    "grandmother_mother", "grandfather_mother"   # maternal grandparents
]


for person in FAMILY:
    try:
        print(f"Processing {person}...")
        unphased_vcf = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_VCF/{person}/all_variants.vcf.gz"
        bam = parse_BAM(f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genomes/BAM/{person}.txt", person)
        output_dir = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results_whatshap/{person}"
        os.makedirs(output_dir, exist_ok=True)
        output_vcf = f"{output_dir}/phased_variants.vcf.gz"

        log_dir = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_whatshap/logs/"
        os.makedirs(log_dir, exist_ok=True)
        output_log_file = f"{log_dir}{person}.out"
        error_log_file = f"{log_dir}{person}.err"
        # delete old log files if it exists
        if os.path.exists(output_log_file):
            os.remove(output_log_file)
        if os.path.exists(error_log_file):
            os.remove(error_log_file)
        
        script = "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_whatshap/phase_haplotype_individual.py"
    
        qsub_cmd = [
                    "qsub",
                    "-N", f"hap-{person}",
                    "-b", "y",
                    "-q", "dell.q",
                    "-o", output_log_file,
                    "-e", error_log_file,
                    "-l", "vf=0.5G",
                    "/home/b/buit/miniconda3/bin/python", script,
                    "-b", bam,
                    "-v", unphased_vcf,
                    "-o", output_vcf
                ]
        print(f"Submitting job for {person} with command: {' '.join(qsub_cmd)}")
        subprocess.run(qsub_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error submitting job for {person}: {e}")
        print("Command output:", e.output)
        print("Command error:", e.stderr)
    except Exception as e:
        print(f"Unexpected error for {person}: {e}")
        print("Please check the command and try again.")
    
    
