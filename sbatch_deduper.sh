#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --job-name=deduper
#SBATCH --out=./results/results_9/slurm/slurm-%j.out
#SBATCH --err=./results/results_9/slurm/slurm-%j.err

# Stop running on error
set -e

input_file=/projects/bgmp/jonat/bioinfo/Bi624/Deduper-tangjon1/deduper_part_3/input_sam_files/sorted/Dataset1.corrected.sorted.sam

num=9

input_file_basename=$(basename $input_file .sam)

/usr/bin/time -v /projects/bgmp/jonat/bioinfo/Bi624/Deduper-tangjon1/tang_deduper.py \
    -f $input_file \
    -u /projects/bgmp/jonat/bioinfo/Bi624/Deduper-tangjon1/STL96.txt \
    -o "/projects/bgmp/jonat/bioinfo/Bi624/Deduper-tangjon1/results/results_${num}/${input_file_basename}.dedup.sam" \
    -c 8 \
    -r quality_median
exit
