#!/bin/bash

#SBATCH --job-name HOMER
#SBATCH -N 1
#SBATCH --partition=128GB
#SBATCH -t 0-12:0:0
#SBATCH -o HOMER_motif.out
#SBATCH -e HOMER_motif.err
#SBATCH --mail-type ALL

module load homer/4.10.4

peakfile=$1
genome=$2
output=$3
size=$4

findMotifsGenome.pl $peakfile $genome $output $size

