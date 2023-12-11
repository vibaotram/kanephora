#!/bin/bash
#SBATCH --job-name get_submat
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
clim=/shared/projects/most_kmer/afterkiss/gea/explanatory_wc2.1.csv 
out=/shared/projects/most_kmer/afterkiss/gea/lfmm_rda/out

singularity exec -B /shared/home/baotram \
/shared/projects/vietcaf/Singularity.R_4-3-Rserver_2022.sif \
Rscript /shared/projects/most_kmer/afterkiss/gea/lfmm_rda/get_submat.R $bed $clim $out