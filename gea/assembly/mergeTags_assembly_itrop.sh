#!/bin/bash
#SBATCH --job-name assembly_signif_kmers_m15
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --partition=highmem
#SBATCH --nodelist=node4
#SBATCH --error /scratch/most_kmer/assembly/slurm-%x.log
#SBATCH --output /shared/most_kmer/assembly/slurm-%x.log

cd /scratch/most_kmer/assembly
./dekupl-mergeTags/mergeTags -k 31 -m 15 -n significant_kmers_lfmm.tsv > mergeTags_significant_kmers_lfmm_m15.tsv
