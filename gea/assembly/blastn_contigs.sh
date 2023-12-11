#!/bin/bash
#SBATCH --job-name blastn_contigs_m15_l500_CC1_8
#SBATCH --cpus-per-task=30
#SBATCH --mem=60G
#SBATCH --partition=long
#SBATCH -A most_kmer
#SBATCH --error /shared/projects/most_kmer/afterkiss/gea/annotation/slurm-%x.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/gea/annotation/slurm-%x.log

module load blast

echo "< $(date) > blastn against reference genome > tsv output (fmt6)"
blastn \
-subject /shared/projects/vietcaf/data/reference/CC1.8_v2_pseudomolecule_cat.fa \
-query /shared/projects/most_kmer/afterkiss/gea/assembly/contigs_m15_l500.fasta \
-out /shared/projects/most_kmer/afterkiss/gea/annotation/blast_contigs_m15_l500_CC1_8.tsv \
-outfmt "6 qacc sacc qstart qend sstart send bitscore score pident gaps mismatch length qlen" \
-evalue 1e6 \
-max_target_seqs 1 \
-num_threads $SLURM_CPUS_PER_TASK
echo "< $(date) > done"
# -outfmt "6 qacc sacc qstart qend sstart send pident gaps mismatch length"


echo "< $(date) > blastn against reference genome > xml output (fmt16)"
blastn \
-subject /shared/projects/vietcaf/data/reference/CC1.8_v2_pseudomolecule_cat.fa \
-query /shared/projects/most_kmer/afterkiss/gea/assembly/contigs_m15_l500.fasta \
-out /shared/projects/most_kmer/afterkiss/gea/annotation/blast_contigs_m15_l500_CC1_8.xml \
-outfmt 16 \
-evalue 1e6 \
-max_target_seqs 1 \
-num_threads $SLURM_CPUS_PER_TASK
echo "< $(date) > done"
