#!/bin/bash
#START_TIME=$(date +"%m%d%H%M" -d "now + 3 hours")

LOGDIR="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_compare_vcf/merge_logs"
CHUNKFILE="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/chromosome_list.txt"
SCRIPT="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_compare_vcf/merge_vcf_chromosome.py"

#rm -f "$LOGDIR"/run_*.out "$LOGDIR"/run_*.err
total_lines=$(wc -l < "$CHUNKFILE")
dell_q_end=73
dell256_q_end=$((73 + 39))  # 112
count=0

while IFS= read -r chunk; do
    rm -f "$LOGDIR"/run_$chunk.out "$LOGDIR"/run_$chunk.err

    if (( count < dell_q_end )); then
        queue="dell.q"
    elif (( count < dell256_q_end )); then
        queue="dell256.q"
    else
        queue="hpclient.q"
    fi

    echo "Submitting $chunk to $queue"

    #qsub -a "$START_TIME" -N mvcf-$chunk -b y -q "$queue" \
    qsub -N mvcf-$chunk -b y -q "$queue" \
        -o "$LOGDIR/run_$chunk.out" \
        -e "$LOGDIR/run_$chunk.err" \
        -l vf=0.5G \
        /home/b/buit/miniconda3/envs/HiWi/bin/python "$SCRIPT" -c "$chunk"

    ((count++))
done < "$CHUNKFILE"