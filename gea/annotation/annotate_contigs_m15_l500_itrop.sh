#!/bin/bash
#SBATCH --job-name annot_ctg_m15_l500
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --partition=highmem
#SBATCH --nodelist=node21
#SBATCH --error /data/projects/rob/tmp/slurm-%x.log
#SBATCH --output /data/projects/rob/tmp/slurm-%x.log

echo "jobid: $SLURM_JOB_ID"
echo "node: $SLURM_NODELIST"
echo "cpus: $SLURM_CPUS_PER_TASK"

module load bioinfo/blast/2.12.0+

## work on a temp dir 
nasdir=/data/projects/rob/tmp
# tmpdir=$(mktemp -d -p /scratch/ --suffix _blastx)
tmpdir=/scratch/tmp.V8bMGi7L38_blastx
cd $tmpdir
echo "< $(date) > working on temp dir: ${tmpdir}"


## download nr database
echo "< $(date) > download nr database"
# update_blastdb.pl --decompress nr --num_threads $SLURM_CPUS_PER_TASK
echo "< $(date) > done"

## get taxid list of the Tracheophyta subkingdom
echo "< $(date) > get taxid list of the Tracheophyta subkingdom"
# get_species_taxids.sh -t 58023 > tracheophyta.txids
echo "< $(date) > done"

## get query file from nas
echo "< $(date) > get query file from nas"
fafile=contigs_m15_l500.fasta
# rsync -vauP /data/projects/rob/tmp/$fafile .
echo "< $(date) > done"

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

## copy results to nas
echo "< $(date) > copy blast output to nas"
rsync -vauP ${outname}_* ${nasdir}/
echo "< $(date) > done"