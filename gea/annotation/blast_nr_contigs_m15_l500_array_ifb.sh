#!/bin/bash
#SBATCH --job-name annot_ctg_m15_l500_besthit
#SBATCH --cpus-per-task=10
#SBATCH --mem=4G
#SBATCH --partition=long
#SBATCH -A most_kmer
#SBATCH --error /shared/projects/most_kmer/afterkiss/gea/annotation/slurm-%x_%a.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/gea/annotation/slurm-%x_%a.log
#SBATCH --array=[3385]

echo "jobid: $SLURM_JOB_ID"
echo "node: $SLURM_NODELIST"
echo "cpus: $SLURM_CPUS_PER_TASK"
echo "array: $SLURM_ARRAY_TASK_ID"

l=$SLURM_ARRAY_TASK_ID

module load blast/2.12.0

## work on a temp dir 
workdir=/shared/projects/most_kmer/afterkiss/gea/annotation/blastx
outdir=/shared/projects/most_kmer/afterkiss/gea/annotation/contigs_m15_l500_blastx
cd $workdir
echo "< $(date) > working on: ${workdir}"
# outdir=$workdir/contigs_m15_l500_blastx
outname=contig_${l}
echo "< $(date) > contig: ${outname}"

## extract contig
fafile=/shared/projects/most_kmer/afterkiss/gea/assembly/contigs_m15_l500.fasta
awk -v seq="contig_${l}" -v RS='>' '$1 == seq {print RS $0}' $fafile > ${outname}.fasta


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
-taxidlist $taxidfile 
-num_threads $SLURM_CPUS_PER_TASK \
-mt_mode 1
echo "< $(date) > output: $blast_xml on $workdir"
echo "< $(date) > moving $blast_xml to $outdir"
mv $blast_xml $outdir/
echo "< $(date) > done"
