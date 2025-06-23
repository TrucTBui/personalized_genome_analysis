import subprocess
import argparse
import time

FAMILY = [
    "grandfather_father", "grandmother_father",  # paternal grandparents
    "father", "child", "mother", "aunt",         # nuclear family
    "grandmother_mother", "grandfather_mother"   # maternal grandparents
]

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Run statistics for a family of samples.")
    argparser.add_argument("-c", "--chromosome", type=str, required=True, help="Chromosome to run statistics for.")
    args = argparser.parse_args()
    chrom = args.chromosome
    analysis_folder = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results/{chrom}/"
    print(f"Running statistics for chromosome: {chrom}")
    for person in FAMILY:
        cmd = [
            "/home/b/buit/miniconda3/envs/HiWi/bin/python",
            "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_statistics_no_dup/create_statistics_no_dup.py",
            "-c", chrom.replace("chr", ""),
            "-p", person        
            ]
        start_time = time.time()
        subprocess.run(cmd, check=True)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Runtime: {elapsed_time:.2f} seconds for {person} on chromosome {chrom}.")