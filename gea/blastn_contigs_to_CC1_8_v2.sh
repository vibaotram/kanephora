#!/bin/bash
#SBATCH --job-name blast_nt_ws4_contigs_m15_l400
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition=normal
#SBATCH --nodelist=node0
#SBATCH --error /data/projects/rob/tmp/slurm-%x.log
#SBATCH --output /data/projects/rob/tmp/slurm-%x.log

echo -e "jobid: $SLURM_JOB_ID"

module load bioinfo/blast/2.12.0+
blastn -version

wd=$(mktemp -d -p /scratch)
cd $wd

rsync -vauP baotram@core.cluster.france-bioinformatique.fr:/shared/projects/most_kmer/afterkiss/gea/assembly/contigs_m15_l400.fasta .

rsync -vauP baotram@core.cluster.france-bioinformatique.fr:/shared/projects/most_kmer/afterkiss/gea/annotation/CC1.8_v2_named_annotation_0.5_BTI.fa .

blastn \
-query contigs_m15_l400.fasta \
-subject CC1.8_v2_named_annotation_0.5_BTI.fa \
-out blastn_contigs_m15_l400_CC1.8_v2_results.txt \
-outfmt "6 qacc qstart qend qlen length sacc sstart send bitscore score evalue pident gaps mismatch" \
-evalue 1e6
# -num_threads $SLURM_CPUS_PER_TASK

# -outfmt "6 qacc sacc qstart qend sstart send pident gaps mismatch length"
