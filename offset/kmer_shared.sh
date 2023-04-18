#!/bin/bash
#SBATCH --job-name kmer_shared_sub
#SBATCH --cpus-per-task=10
#SBATCH --mem=5G
#SBATCH --partition=fast
#SBATCH -A most_kmer
#SBATCH --error /shared/projects/most_kmer/afterkiss/offset/kmer_shared/log/slurm-%x_%a.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/offset/kmer_shared/log/slurm-%x_%a.log
#SBATCH --array=[843-1039]

module load r/4.1.1

Rscript --vanilla /shared/ifbstor1/projects/most_kmer/afterkiss/offset/kmer_shared.R \
-t /shared/ifbstor1/projects/most_kmer/afterkiss/offset/signif_kmer_table.txt \
-i /shared/ifbstor1/projects/most_kmer/afterkiss/offset/all_samples_info.txt \
-n 103800513 \
-s 100000 \
-b $SLURM_ARRAY_TASK_ID \
-o /shared/projects/most_kmer/afterkiss/offset/kmer_shared \
-c $SLURM_CPUS_PER_TASK 
