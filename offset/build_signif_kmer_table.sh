#!/bin/bash
#SBATCH --job-name signif_kmer_table
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --partition=long
#SBATCH -A most_kmer
#SBATCH --error /shared/projects/most_kmer/afterkiss/offset/slurm-%x_%a.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/offset/slurm-%x_%a.log


module load singularity 

singularity shell -B /shared/projects/most_kmer /shared/projects/most_kmer/KIsS/kiss/containers/Singularity.kiss_tools.sif

echo "<$(date)> "
cd /shared/projects/most_kmer/afterkiss/offset
filter_kmers -t /shared/projects/most_kmer/kiss_all/output/2.KMERS_TABLE/kmers_table \
-k significant_kmers_bio.txt \
-o signif_kmer_table.txt