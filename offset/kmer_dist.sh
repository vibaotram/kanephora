#!/bin/bash
#SBATCH --job-name kmer_dist_sub
#SBATCH --cpus-per-task=10
#SBATCH --mem=5G
#SBATCH --partition=fast
#SBATCH -A most_kmer
#SBATCH --error /shared/projects/most_kmer/afterkiss/offset/kmer_dist/log/slurm-%x_%a.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/offset/kmer_dist/log/slurm-%x_%a.log
#SBATCH --array=[1-1005]

module load r/4.1.1

Rscript --vanilla /shared/ifbstor1/projects/most_kmer/afterkiss/offset/kmer_dist.R \
100489652 \
100000 \
$SLURM_ARRAY_TASK_ID \
/shared/projects/most_kmer/afterkiss/offset/kmer_dist \
$SLURM_CPUS_PER_TASK 
