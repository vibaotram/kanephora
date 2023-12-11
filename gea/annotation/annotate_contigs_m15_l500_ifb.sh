#!/bin/bash
#SBATCH --job-name blast_nr_ctg_m15_l500
#SBATCH --cpus-per-task=60
#SBATCH --mem=20G
#SBATCH --partition=long
#SBATCH --error /shared/projects/most_kmer/afterkiss/gea/annotation/slurm-%x.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/gea/annotation/slurm-%x.log

echo "jobid: $SLURM_JOB_ID"
echo "node: $SLURM_NODELIST"
echo "cpus: $SLURM_CPUS_PER_TASK"

module load blast/2.12.0

## work on a temp dir 
workdir=/shared/projects/most_kmer/afterkiss/gea/annotation
cd $workdir
echo "< $(date) > working on: ${workdir}"


## download nr database
echo "< $(date) > download nr database"
# update_blastdb.pl --decompress nr --num_threads $SLURM_CPUS_PER_TASK
echo "< $(date) > done"

## get taxid list of the Tracheophyta subkingdom
echo "< $(date) > get taxid list of the Tracheophyta subkingdom"
# get_species_taxids.sh -t 58023 > tracheophyta.txids
echo "< $(date) > done"


fafile=contigs_m15_l500.fasta
outname=contigs_m15_l500


## blastx output xml
echo "< $(date) > blastx output tsv"
blast_tsv=${outname}_blastx.tsv
blastx \
-db nr \
-query $fafile \
-out $blast_tsv \
-outfmt "6 qacc qgi qstart qend sacc sgi sstart send bitscore score evalue pident gaps mismatch length qlen staxid ssciname sblastname stitle" \
-evalue 1e6 \
-max_target_seqs 1 \
-show_gis \
-taxidlist tracheophyta.txids \
-num_threads $SLURM_CPUS_PER_TASK 
echo "< $(date) > done"


echo "< $(date) > blastx output xml"
blast_xml=${outname}_blastx.xml
blastx \
-db nr \
-query $fafile \
-out $blast_xml \
-outfmt 5 \
-evalue 1e6 \
-max_target_seqs 1 \
-show_gis \
-taxidlist tracheophyta.txids \
-num_threads $SLURM_CPUS_PER_TASK 
echo "< $(date) > done"
