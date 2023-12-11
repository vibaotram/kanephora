#!/bin/bash
#SBATCH --job-name assembly_signif_kmers_m15
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --partition=long
#SBATCH -A elai_most
#SBATCH --error /shared/projects/most_kmer/afterkiss/gea/assembly/slurm-%x.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/gea/assembly/slurm-%x.log

cd /shared/projects/most_kmer/afterkiss/gea/assembly/
./dekupl-mergeTags/mergeTags -k 31 -m 15 -n candidates_rda_lfmm.txt > mergeTags_significant_lfmm_rda.tsv
