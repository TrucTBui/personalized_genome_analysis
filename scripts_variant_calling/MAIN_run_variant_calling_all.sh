#!/bin/bash

LOGDIR="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_variant_calling/logs"
CHUNKFILE="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/chunk_list.txt"
SCRIPT="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_variant_calling/run_call_variant_family.py"

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

    qsub -N vc-$chunk -b y -q "$queue" \
        -o "$LOGDIR/run_$chunk.out" \
        -e "$LOGDIR/run_$chunk.err" \
        -l vf=0.5G \
        /home/b/buit/miniconda3/envs/HiWi/bin/python "$SCRIPT" -c "$chunk"

    ((count++))
done < "$CHUNKFILE"