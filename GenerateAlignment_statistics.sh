#!/bin/bash

#SBATCH --job-name Alignment_Stats
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --partition=32GB
#SBATCH -t 0-12:0:0
#SBATCH -o Alignment_stats.out
#SBATCH -e Alignment_stats.err

module load samtools/1.6

resultdir=$1
rawdir=$2
samplesheet=$3

python GenerateAlignment_statistics.py  -d $resultdir -r $rawdir -s $samplesheet