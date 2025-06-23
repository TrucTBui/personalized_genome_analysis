#!/bin/bash

PANGENOME_DIR="/mnt/raidproj/proj/projekte/personalizedmed/pangenome/analysis_schertler_stroebele_ws24"

# Each sub directory in the PANGENOME_DIR is a pangenome analysis for a chromosome
# Run run_transform_pangenome_chrom.sh for each chromosome, using qsub

rm -r /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_process_pangenome/logs/run_transform_pangenome_*.out
rm -r /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_process_pangenome/logs/run_transform_pangenome_*.err

for CHROM in "$PANGENOME_DIR"/*; do
    if [ -d "$CHROM" ]; then
        CHROM_NAME=$(basename "$CHROM")
        echo "Submitting job for $CHROM_NAME"
        
        qsub -N pangenome_$CHROM_NAME -b y -q dell.q \
            -o /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_process_pangenome/logs/run_transform_pangenome_$CHROM_NAME.out \
            -e /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_process_pangenome/logs/run_transform_pangenome_$CHROM_NAME.err \
            -l vf=0.5G \
            /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_process_pangenome/run_transform_pangenome_chrom.sh "$CHROM_NAME"
    fi
done


