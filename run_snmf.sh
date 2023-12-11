#!/bin/bash
#SBATCH --job-name run_snmf
#SBATCH --partition=long
#SBATCH --cpus-per-task=50
#SBATCH --mem=40G
#SBATCH --output /shared/projects/most_kmer/afterkiss/snmf/slurm-%x.log
#SBATCH --error /shared/projects/most_kmer/afterkiss/snmf/slurm-%x.log


module load r/4.2.1

Rscript /shared/ifbstor1/projects/most_kmer/afterkiss/run_snmf.R