import os
import gzip
from collections import defaultdict
import argparse
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--chromosome", type=str, required=True, help="chromosome to process")


args = parser.parse_args()
chromosome = args.chromosome


def parse_gtf_to_tsv(gtf_file, chromosome, output_dir, chunk_size=500):
    os.makedirs(output_dir, exist_ok=True)
    
    gene_id_dir = os.path.join(output_dir, "gene_id")
    os.makedirs(gene_id_dir, exist_ok=True)
    
    gene_id_with_name_dir = os.path.join(output_dir, "gene_id_with_name")
    os.makedirs(gene_id_with_name_dir, exist_ok=True)

    gtf_all_genes_dir = os.path.join(output_dir, "gtf_all_genes")
    os.makedirs(gtf_all_genes_dir, exist_ok=True)


    gene_lines = defaultdict(list)


    # Read GTF file
    with gzip.open(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split('\t')
            chr = fields[0]
            if chr != chromosome:
                continue
            
            attr_field = fields[8]
            try:
                attributes = {k: v.strip('"') for k, v in 
                              (attr.strip().split() for attr in attr_field.split(';') if attr.strip())}
            except ValueError:
                continue  # Skip malformed lines

            gene_id = attributes.get("gene_id", "Unknown")
            gene_lines[gene_id].append(line.strip())


    gene_id_list = []
    for gene_id, _ in gene_lines.items():
        gene_id_list.append(gene_id)
    
    chunks = [gene_id_list[i:i + chunk_size] for i in range(0, len(gene_id_list), chunk_size)]
    if chromosome.startswith("chr"):
        chromosome_name = chromosome
    else:
        chromosome_name = "chr" + chromosome

    for chunk_num, chunk in enumerate(chunks, 1):
        gtf_all_genes_chunk = os.path.join(gtf_all_genes_dir, f"{chromosome_name}_chunk{chunk_num}")
        os.makedirs(gtf_all_genes_chunk, exist_ok=True)

        for gene_id in chunk:
            gene_lines_chunk = gene_lines[gene_id]
            with open(os.path.join(gtf_all_genes_chunk, f"{gene_id}.tsv"), 'w') as out_file:
                out_file.write("\n".join(gene_lines_chunk))

        gene_id_file = os.path.join(gene_id_dir, f"gene_id_list_{chromosome_name}_chunk{chunk_num}.txt")
        with open(gene_id_file, "w") as f:
            f.write("\n".join(chunk))

        gene_id_with_name_file = os.path.join(gene_id_with_name_dir, f"gene_id_with_name_{chromosome_name}_chunk{chunk_num}.txt")
        cmd_map_gene_name = f"Rscript /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_process_isar/Ensembl2Gene.R -i {gene_id_file} -o {gene_id_with_name_file}"
        subprocess.run(cmd_map_gene_name, shell=True, check=True)

    print(f"Finished mapping gene names for chromosome {chromosome_name}")

    return gene_lines  # optionally return for inspection

parse_gtf_to_tsv("/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genes/isar.ensembl-75.gtf.gz", chromosome, "/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/input_genes/")