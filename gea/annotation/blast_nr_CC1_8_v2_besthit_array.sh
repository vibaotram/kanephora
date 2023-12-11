#!/bin/bash
#SBATCH --job-name blast_nr_CC1_8_v2_besthit_fmt16
#SBATCH --cpus-per-task=20
#SBATCH --mem=10G
#SBATCH --partition=long
#SBATCH --exclude=cpu-node-33
#SBATCH -A most_kmer
#SBATCH --error /shared/projects/most_kmer/afterkiss/gea/annotation/log/slurm-%x_%a.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/gea/annotation/log/slurm-%x_%a.log
#SBATCH --array=[1061,1550,757]

echo "jobid: $SLURM_JOB_ID"
echo "node: $SLURM_NODELIST"
echo "cpus: $SLURM_CPUS_PER_TASK"
echo "array: $SLURM_ARRAY_TASK_ID"

a=$SLURM_ARRAY_TASK_ID
n=10 # number of seqs extracted in one fasta file
l=$(($a*$n*2))

module load blast/2.12.0

## work on a temp dir 
workdir=/shared/projects/most_kmer/afterkiss/gea/annotation/blastx
outdir=/shared/projects/most_kmer/afterkiss/gea/annotation/CC1.8_v2_blastx
cd $workdir
echo "< $(date) > working on: ${workdir}"
# outdir=$workdir/contigs_m15_l500_blastx
outname=CC1.8_v2_batch_${a}
echo "< $(date) > batch: ${a}"

## get fasta sequences from canephora gff
# echo "$(date) > get fasta sequences from canephora gff"
# cd $workdir
# infa=/shared/projects/vietcaf/data/reference/CC1.8_v2_pseudomolecule_cat.fa
# cangff=/shared/projects/most_kmer/afterkiss/gea/annotation/CC1.8_v2_named_annotation_0.5_BTI_unique.gff3
outfa=/shared/projects/most_kmer/afterkiss/gea/annotation/CC1.8_v2_named_annotation_0.5_BTI.fa
# bedtools getfasta -fi $infa -bed $cangff -fo $outfa


## extract contig
echo "$(date) > extracting $n sequences"
fafile=/shared/projects/most_kmer/afterkiss/gea/annotation/CC1.8_v2_named_annotation_0.5_BTI.fa
cat $fafile | head -n $l | tail -n $(($n*2)) > ${outname}.fasta
echo "< $(date) > done"

## download nr database
echo "< $(date) > download nr database"
# update_blastdb.pl --decompress nr --num_threads $SLURM_CPUS_PER_TASK
echo "< $(date) > done"

## get taxid list of the Tracheophyta subkingdom
echo "< $(date) > get taxid list of Viridiplantae"
taxidfile=Viridiplantae.txids
# get_species_taxids.sh -t 33090 > $taxidfile
echo "< $(date) > done"



## blastx output xml2
echo "< $(date) > blastx output xml format 16"
blast_xml=${outname}_blastx_besthit.xml
blastx \
-db nr \
-query ${outname}.fasta \
-out $blast_xml \
-outfmt 16 \
-evalue 1e6 \
-max_target_seqs 1 \
-show_gis \
-taxidlist $taxidfile \
-num_threads $SLURM_CPUS_PER_TASK
echo "< $(date) > output: $blast_xml on $workdir"
echo "< $(date) > moving $blast_xml to $outdir"
mv $blast_xml $outdir/
echo "< $(date) > done"

# echo "$(date) > blastx output tsv"
# output_tsv=${outname}_blastx_besthit.tsv
# blastx \
# -db nr \
# -query ${outname}.fasta \
# -out $output_tsv \
# -outfmt "6 qacc qgi qstart qend sacc sgi sstart send bitscore score evalue pident gaps mismatch length qlen staxid ssciname sblastname stitle" \
# -evalue 1e6 \
# -show_gis \
# -max_target_seqs 1 \
# -taxidlist $taxidfile \
# -num_threads $SLURM_CPUS_PER_TASK \
# -mt_mode 1
# echo "< $(date) > done"

