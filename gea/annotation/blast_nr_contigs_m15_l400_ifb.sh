#!/bin/bash
#SBATCH --job-name blast_nr_contigs_m15_l400
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

## get taxid list of the Tracheophyta subkingdom
echo "get taxid list of the Tracheophyta subkingdom"
get_species_taxids.sh -t 58023 > tracheophyta.txids

## query file
fafile=/shared/projects/most_kmer/afterkiss/gea/assembly/contigs_m15_l400.fasta

## output name
outname=blast_nr_contigs_m15_l400

## blastx output xml
echo "blastx output xml"
blastx \
-db nr \
-query $fafile \
-out ${outname}.xml \
-outfmt 5 \
-evalue 1e6 \
-max_target_seqs 1 \
-taxidlist tracheophyta.txids \
-num_threads $SLURM_CPUS_PER_TASK 
