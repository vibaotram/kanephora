#!/bin/bash
#SBATCH --job-name kiss_af
#SBATCH --cpus-per-task 20
#SBATCH --partition=highmem
#SBATCH --mem 60G
#SBATCH --nodelist node4
#SBATCH --error /scratch/most_kmer/kiss_af/log/slurm-%x.log
#SBATCH --output /scratch/most_kmer/kiss_af/log/slurm-%x.log


module load system/singularity/3.6.0
module load system/python/3.8.12

# snakemake --unlock

kiss run_local -c /scratch/most_kmer/kiss_af/configs/workflow_config_itrop.yaml --dry-run -t $SLURM_CPUS_PER_TASK --forcerun kmer_position_from_bam
