#!/bin/bash
#SBATCH --job-name download_af_present
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --partition=fast
#SBATCH -A elai_most
#SBATCH --error /shared/projects/most_kmer/afterkiss/bioclim/slurm-%x.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/bioclim/slurm-%x.log

module unload r
module load r/4.1.1

cd /shared/projects/most_kmer/afterkiss/bioclim

Rscript ./download_bioclim_WC_AF.R wc1-4_present