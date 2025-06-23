#!/bin/bash

rm -f /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_process_isar/logs/run_isar.*.out
rm -f /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_process_isar/logs/run_isar.*.err

for chr in {1..22} X Y; do
    echo "submitted $chr"
    qsub -N truc-$chr -b y -q dell.q -o /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_process_isar/logs/run_isar.$chr.out -e /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_process_isar/logs/run_isar.$chr.err -l vf=1G /mnt/raidproj/proj/projekte/personalizedmed/PPG/miRNAs/scripts_process_isar/run_isar_chr.sh $chr
done

