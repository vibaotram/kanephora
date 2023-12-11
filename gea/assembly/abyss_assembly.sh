#!/bin/bash
#SBATCH --job-name abyss25_signif_kmers
#SBATCH --cpus-per-task=40
#SBATCH --mem=60G
#SBATCH --partition=long
#SBATCH -A elai_most
#SBATCH --error /shared/projects/most_kmer/afterkiss/gea/assembly/slurm-%x.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/gea/assembly/slurm-%x.log

module load abyss/2.2.1

cd /shared/projects/most_kmer/afterkiss/gea/assembly/
# abyss-pe k=31 c=0 e=0 np=$SLURM_CPUS_PER_TASK --oversubscribe \
# name=abyss_contigs_k31 \
# in=/shared/ifbstor1/projects/most_kmer/afterkiss/gea/assembly/significant_kmers.fasta

mpirun -np 20 --oversubscribe ABYSS-P -k25 -q3 -e0 -c0 \
--coverage-hist=coverage.hist \
-s abyss_contigs_k25-bubbles.fa  \
-o abyss_contigs_k25.fa \
/shared/ifbstor1/projects/most_kmer/afterkiss/gea/assembly/significant_kmers.fasta