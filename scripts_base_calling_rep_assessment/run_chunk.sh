#!/bin/bash

chunk=$1
echo "Processing chunk $chunk"
python3 /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_base_calling_rep_assessment/run.py -c $chunk 
