#!/bin/bash
#SBATCH --job-name kiss_all
#SBATCH --partition=long
#SBATCH -A elai_most
#SBATCH --error /shared/projects/most_kmer/kiss_all/log/slurm-%x.log
#SBATCH --output /shared/projects/most_kmer/kiss_all/log/slurm-%x.log

# rsync -vaurP /shared/ifbstor1/projects/most_kmer/kiss_af/output/1.KMERS_GWAS /shared/projects/most_kmer/kiss_allput/

module load singularity
module load python/3.9 
module load snakemake/7.7.0

# snakemake --unlock

kiss run_cluster -c /shared/projects/most_kmer/kiss_all/configs/workflow_config.yaml \
-cl /shared/projects/most_kmer/kiss_all/configs/cluster_config.yaml \
--singularity-args \" --bind /shared:/shared \" --nolock \
--forcerun kmers_to_use 
