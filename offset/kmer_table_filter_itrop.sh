#!/bin/bash
#SBATCH --job-name kmer_table_sub
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition=highmem
#SBATCH --nodelist=node4
#SBATCH --error /scratch/most_kmer/filter_kmers/log/slurm-%x_%a.log
#SBATCH --output /scratch/most_kmer/filter_kmers/log/slurm-%x_%a.log
#SBATCH --array=[1-10]

module load singularity

cd /scratch/most_kmer/filter_kmers

kmer_list=$(echo $(ls x*.txt) | cut -d' ' -f $SLURM_ARRAY_TASK_ID)
singularity exec -B /shared/projects/most_kmer \
/data/projects/rob/scripts/KIsS/kiss/containers/Singularity.kiss_tools.sif \
filter_kmers -t kmers_table \
-k $kmer_list \
-o signif_kmer_table_x05_${SLURM_ARRAY_TASK_ID}.txt