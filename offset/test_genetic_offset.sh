#!/bin/bash
#SBATCH --job-name test_genoff
#SBATCH --cpus-per-task=5
#SBATCH --mem=120G
#SBATCH --partition=bigmem
#SBATCH -A elai_most
#SBATCH --error /shared/projects/most_kmer/afterkiss/offset/log/slurm-%x_%j.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/offset/log/slurm-%x_%j.log


module load r/4.2.3

Rscript /shared/projects/most_kmer/afterkiss/offset/test_genetic_offset.R

## check job 33399990
