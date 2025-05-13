#!/bin/bash

LOGDIR="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_base_calling_rep_assessment/logs"
CHUNKFILE="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/chunk_list_Y.txt"
SCRIPT="/mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_base_calling_rep_assessment/run_chunk.sh"

rm -f "$LOGDIR"/run_*.out "$LOGDIR"/run_*.err
total_lines=$(wc -l < "$CHUNKFILE")
dell_q_end=73
dell256_q_end=$((73 + 39))  # 112
count=0

while IFS= read -r chunk; do
    if (( count < dell_q_end )); then
        queue="dell.q"
    elif (( count < dell256_q_end )); then
        queue="dell256.q"
    else
        queue="hpclient.q"
    fi

    echo "Submitting $chunk to $queue"

    qsub -N bc-$chunk -b y -q "$queue" \
        -o "$LOGDIR/run_$chunk.out" \
        -e "$LOGDIR/run_$chunk.err" \
        -l vf=0.5G \
        "$SCRIPT" "$chunk"

    ((count++))
done < "$CHUNKFILE"