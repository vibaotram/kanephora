#!/bin/bash
#SBATCH --job-name kmer_table_sub
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition=long
#SBATCH -A most_kmer
#SBATCH --error /shared/projects/most_kmer/afterkiss/offset/kmer_table/log/slurm-%x_%a.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/offset/kmer_table/log/slurm-%x_%a.log
#SBATCH --array=[1-100]

module load singularity

cd /shared/projects/most_kmer/afterkiss/offset/kmer_table

kmer_list=$(echo $(ls x*.txt) | cut -d' ' -f $SLURM_ARRAY_TASK_ID)
singularity exec -B /shared/projects/most_kmer \
/shared/projects/most_kmer/KIsS/kiss/containers/Singularity.kiss_tools.sif \
filter_kmers -t /shared/projects/most_kmer/kiss_all/output/2.KMERS_TABLE/kmers_table \
-k $kmer_list \
-o signif_kmer_table_x05_${SLURM_ARRAY_TASK_ID}.txt