#!/bin/bash
#SBATCH --job-name kiss_af
#SBATCH --partition=long
#SBATCH -A elai_most
#SBATCH --error /shared/projects/most_kmer/kiss_af/log/slurm-%x.log
#SBATCH --output /shared/projects/most_kmer/kiss_af/log/slurm-%x.log


module load singularity
module load python/3.9 
module load snakemake/7.7.0

# snakemake --unlock

kiss run_cluster -c /shared/projects/most_kmer/kiss_af/configs/workflow_config.yaml \
-cl /shared/projects/most_kmer/kiss_af/configs/cluster_config.yaml \
--singularity-args \" --bind /shared:/shared \" --nolock --dry-run
