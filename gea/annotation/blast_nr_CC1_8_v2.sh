#!/bin/bash
#SBATCH --job-name blast_nr_CC1_8_v2
#SBATCH --cpus-per-task=50
#SBATCH --mem=20G
#SBATCH --partition=long
#SBATCH --error /shared/ifbstor1/projects/most_kmer/afterkiss/gea/annotation/slurm-%x.log
#SBATCH --output /shared/ifbstor1/projects/most_kmer/afterkiss/gea/annotation/slurm-%x.log

echo "jobid: $SLURM_JOB_ID"
echo "node: $SLURM_NODELIST"
echo "cpus: $SLURM_CPUS_PER_TASK"

workdir=/shared/projects/most_kmer/afterkiss/gea/annotation

## load modules
echo "$(date) > load modules"
# module load bedtools/2.30.0
# bedtools -version
module load blast/2.12.0
blastx -version

## get fasta sequences from canephora gff
# echo "$(date) > get fasta sequences from canephora gff"
cd $workdir
infa=/shared/projects/vietcaf/data/reference/CC1.8_v2_pseudomolecule_cat.fa
cangff=/shared/projects/most_kmer/afterkiss/gea/annotation/CC1.8_v2_named_annotation_0.5_BTI.gff3
outfa=CC1.8_v2_named_annotation_0.5_BTI.fa
# bedtools getfasta -fi $infa -bed $cangff -fo $outfa

## blastx canephora genes
echo "$(date) > blastx output xml"
output=CC1.8_v2_named_annotation_0.5_BTI_blastx.xml
blastx \
-db nr \
-query $outfa \
-out $output \
-outfmt 16 \
-evalue 1e6 \
-max_target_seqs 5 \
-taxidlist tracheophyta.txids \
-num_threads $SLURM_CPUS_PER_TASK 
