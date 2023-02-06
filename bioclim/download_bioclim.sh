#!/bin/bash
#SBATCH --job-name download_bioclim_vn
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --partition=long
#SBATCH -A elai_most
#SBATCH --error /shared/projects/most_kmer/afterkiss/bioclim/slurm-%x.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/bioclim/slurm-%x.log

module unload r
module load r/4.1.0

Rscript /shared/projects/most_kmer/afterkiss/bioclim/download_bioclim_WC_VN.R
