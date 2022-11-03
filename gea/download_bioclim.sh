#!/bin/bash
#SBATCH --job-name download_bioclim
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --partition=long
#SBATCH -A elai_most
#SBATCH --error /shared/projects/most_kmer/afterkiss/gea/slurm-%x.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/gea/slurm-%x.log

module unload r
module load r/4.1.0

Rscript /shared/projects/most_kmer/afterkiss/gea/download_bioclim.R