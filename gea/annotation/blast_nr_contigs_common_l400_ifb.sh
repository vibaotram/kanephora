#!/bin/bash
#SBATCH --job-name blast_nr_contigs_common_l400
#SBATCH --cpus-per-task=50
#SBATCH --mem=20G
#SBATCH --partition=fast
#SBATCH --error /shared/ifbstor1/projects/most_kmer/afterkiss/gea/annotation/slurm-%x.log
#SBATCH --output /shared/ifbstor1/projects/most_kmer/afterkiss/gea/annotation/slurm-%x.log

echo "jobid: $SLURM_JOB_ID"
echo "node: $SLURM_NODELIST"
echo "cpus: $SLURM_CPUS_PER_TASK"

module load blast/2.12.0

echo "Start blast"
## work dir
workdir=/shared/projects/most_kmer/afterkiss/gea/annotation/blastx
cd $workdir

## get taxid list of the Viridiplantae
echo "get taxid list of the Viridiplantae"
taxidfile=Viridiplantae.txids
# get_species_taxids.sh -t 58023 > $taxidfile

## query file
fafile=/shared/projects/most_kmer/afterkiss/gea/assembly/contigs_common_l400.fasta

## output name
outname=blast_nr_contigs_common_l400_besthit

## blastx output xml
# echo "blastx output xml format 14"
# blastx \
# -db nr \
# -query $fafile \
# -out ${outname}.xml \
# -outfmt 14 \
# -evalue 1e-6 \
# -max_target_seqs 1 \
# -taxidlist $taxidfile \
# -num_threads $SLURM_CPUS_PER_TASK 
# echo "< $(date) > done"


## blastx output tsv
echo "$(date) > blastx output tsv"
blastx \
-db nr \
-query $fafile \
-out ${outname}.tsv \
-outfmt "6 qacc qgi qstart qend sacc sgi sstart send bitscore score evalue pident gaps mismatch length qlen staxid ssciname sblastname stitle" \
-evalue 1e-6 \
-max_target_seqs 1 \
-taxidlist $taxidfile \
-num_threads $SLURM_CPUS_PER_TASK 
echo "< $(date) > done"
