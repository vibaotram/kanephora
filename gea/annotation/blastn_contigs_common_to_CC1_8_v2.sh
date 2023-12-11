#!/bin/bash
#SBATCH --job-name blast_contigs_common_l400_to_CC1_8_v2
#SBATCH --cpus-per-task=50
#SBATCH --mem=20G
#SBATCH --partition=long
#SBATCH --error /shared/ifbstor1/projects/most_kmer/afterkiss/gea/annotation/slurm-%x.log
#SBATCH --output /shared/ifbstor1/projects/most_kmer/afterkiss/gea/annotation/slurm-%x.log

echo "jobid: $SLURM_JOB_ID"
echo "node: $SLURM_NODELIST"
echo "cpus: $SLURM_CPUS_PER_TASK"

module load blast/2.12.0

## work dir
workdir=/shared/projects/most_kmer/afterkiss/gea/annotation
cd $workdir

contig=/shared/projects/most_kmer/afterkiss/gea/assembly/contigs_common_l400.fasta

ref=CC1.8_v2_named_annotation_0.5_BTI.fa

blastn \
-query $contig \
-subject $ref \
-out blastn_contigs_common_l400_CC1.8_v2_results.txt \
-outfmt "6 qacc qstart qend qlen length sacc sstart send bitscore score evalue pident gaps mismatch" \
-evalue 1e-6
