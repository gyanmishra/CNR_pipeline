#!/bin/bash

#SBATCH --job-name Alignment_Stats
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --partition=128GB
#SBATCH -t 0-12:0:0
#SBATCH -o mCA_vs_CNR.out
#SBATCH -e mCA_vs_CNR.err

module load samtools/1.6

echo "$SLURM_SUBMIT_DIR"
script_path="${SLURM_SUBMIT_DIR}/$(basename $0)"

echo "${script_path}"

resultdir=$1
rawdir=$2
samplesheet=$3

python ${script_path}/GenerateAlignment_statistics.py  -d $resultdir -r $rawdir -s $samplesheet