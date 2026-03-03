#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --job-name=sam_sort
#SBATCH --out=./input_sam_files/sorted/slurm/slurm-%j.out
#SBATCH --err=./input_sam_files/sorted/slurm/slurm-%j.err

# Stop running on error
#set -e

for sam_file in /projects/bgmp/jonat/bioinfo/Bi624/Deduper-tangjon1/deduper_part_3/input_sam_files/original/*.sam; do
    file_basename=$(basename $sam_file .sam)
    mamba run -p /projects/bgmp/shared/groups/2025/chimera/envs/samtools/ /usr/bin/time -v samtools sort \
        -o "/projects/bgmp/jonat/bioinfo/Bi624/Deduper-tangjon1/deduper_part_3/input_sam_files/sorted/${file_basename}.sorted.sam" \
        $sam_file
        
done

exit
