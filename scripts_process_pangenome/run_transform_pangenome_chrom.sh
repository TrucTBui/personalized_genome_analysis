#!/bin/bash

CHROM=$1
PANGENOME_DIR="/mnt/raidproj/proj/projekte/personalizedmed/pangenome/analysis_schertler_stroebele_ws24/$CHROM"

# Input will be all files ending with _pangenome_analysis_grouped.tsv
# E.g: chrom1/ENSG00000000457.14/ENSG00000000457.14_pangenome_analysis_grouped.tsv
# Run pangenome2bed.py for each gene

for dir in "$PANGENOME_DIR"/*; do
    if [ -d "$dir" ]; then
        for file in "$dir"/*_pangenome_analysis_grouped.tsv; do
            if [ -f "$file" ]; then
                echo "Processing $file"
                # Run pangenome2bed.py on the file
                /home/b/buit/miniconda3/envs/HiWi/bin/python /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_process_pangenome/pangenome2bed.py -i "$file" 
            fi
        done
    fi
done