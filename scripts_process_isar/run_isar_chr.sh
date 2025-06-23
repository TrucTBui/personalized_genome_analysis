#!/bin/bash

chr=$1
echo "Processing chromosome $chr"
python3 /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_process_isar/process_ISAR.py -c $chr 

