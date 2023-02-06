#!/bin/bash
#SBATCH --job-name download_bioclim_vn
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --partition=long
#SBATCH -A elai_most
#SBATCH --error /shared/projects/most_kmer/afterkiss/offset/slurm-%x.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/offset/slurm-%x.log

module unload r
module load r/4.1.0

Rscript /shared/projects/most_kmer/afterkiss/offset/download_bioclim_vn.R