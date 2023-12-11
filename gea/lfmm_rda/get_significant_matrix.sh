#!/bin/bash
#SBATCH --job-name signif_matrix
#SBATCH --cpus-per-task=5
#SBATCH --mem=20G
#SBATCH --partition=fast
#SBATCH -A elai_most
#SBATCH --error /shared/projects/most_kmer/afterkiss/gea/lfmm_rda/log/slurm-%x_%a.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/gea/lfmm_rda/log/slurm-%x_%a.log
#SBATCH --array=[0-18]


module load singularity

range=$SLURM_ARRAY_TASK_ID

bed=/shared/projects/most_kmer/kiss_af/output/3.TABLE2BED/output_file.${range}.bed

singularity exec -B /shared/home/baotram \
/shared/projects/vietcaf/rserver/Singularity.R_4-2-0-Rserver_2022.sif \
Rscript /shared/ifbstor1/projects/most_kmer/afterkiss/gea/lfmm_rda/get_significant_matrix.R $bed