#!/bin/bash

START_TIME=$(date +"%m%d%H%M" -d "now + 2 minutes")


LOGDIR="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_compare_vcf/compare_logs"
SCRIPT="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_compare_vcf/MAIN_compare_all_vcf.py"

# remove log files if they exist
rm -f "$LOGDIR"/run_compare.out "$LOGDIR"/run_compare.err

qsub -a "$START_TIME" -N compare -b y -q "dell.q" \
        -o "$LOGDIR/run_compare.out" \
        -e "$LOGDIR/run_compare.err" \
        -l vf=0.5G \
        /home/b/buit/miniconda3/envs/HiWi/bin/python "$SCRIPT" 
