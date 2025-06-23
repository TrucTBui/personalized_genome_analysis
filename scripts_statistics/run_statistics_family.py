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
    argparser.add_argument("-c", "--chunk", type=str, required=True, help="Chunk to run statistics for.")
    args = argparser.parse_args()
    chunk = args.chunk
    chrom = chunk.split("_")[0]
    analysis_folder = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/Results/{chrom}/"
    gene_chunk_file = f"/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genes/gene_id/gene_id_list_{chunk}.txt"
    print(f"Running statistics for chunk: {chunk} on chromosome: {chrom}")
    for person in FAMILY:
        cmd = [
            "/home/b/buit/miniconda3/envs/HiWi/bin/python",
            "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_statistics/create_statistics.py",
            "-c", chrom.replace("chr", ""),
            "-p", person,
            "-g", gene_chunk_file
        ]
        start_time = time.time()
        subprocess.run(cmd, check=True)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Statistics for {person} in chunk {chunk} completed in {elapsed_time:.2f} seconds.")