#!/bin/bash
#SBATCH --job-name copy_kiss_out
#SBATCH --cpus-per-task=1
#SBATCH -p supermem
#SBATCH --error /data/projects/rob/scripts/kmer/gea/slurm-%x_%j.log
#SBATCH --output /data/projects/rob/scripts/kmer/gea/slurm-%x_%j.log


rsync -vaurP /scratch/kmer/af_kiss/results_ci2/3.TABLE2BED /data/projects/rob/scripts/kmer/gea/tmp
rsync -vaurP /scratch/kmer/af_kiss/results_ci2/5.RANGES /data/projects/rob/scripts/kmer/gea/tmp