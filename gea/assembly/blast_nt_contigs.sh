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

blastn \
-db nt -remote \
-query /scratch/tmp.O4skjjxWYL_baotram/contigs_m15_l400.fasta \
-out /scratch/tmp.O4skjjxWYL_baotram/blast_nt_ws4_contigs_m15_l400_results.txt \
-outfmt "6 qacc sacc qstart qend sstart send bitscore score evalue pident gaps mismatch length qlen staxid ssciname sblastname stitle" \
-evalue 1e5 \
-word_size 4
# -num_threads $SLURM_CPUS_PER_TASK

# -outfmt "6 qacc sacc qstart qend sstart send pident gaps mismatch length"
