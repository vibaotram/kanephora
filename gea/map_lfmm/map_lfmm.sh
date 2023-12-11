#!/bin/bash
## SBATCH --job-name map_lfmm_bio10
#SBATCH --mem=30G
#SBATCH --partition=long
#SBATCH -A most_kmer
#SBATCH --error /shared/projects/most_kmer/afterkiss/gea/map_lfmm/log/slurm-%x_%a.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/gea/map_lfmm/log/slurm-%x_%a.log
#SBATCH --array=[0-18]

#### 1. merge lfmm pvalues from 100 batches of 1 file into a temp file
#### 2. map the pvalue from temp file to kmer position in the corresponding file
#### 3. seperate the kmer/position/pvalue file by chromosome


## file number
file=$SLURM_ARRAY_TASK_ID
## bioclim variable (1-19)
var=$1
## output dir
dir=/shared/projects/most_kmer/afterkiss/gea/map_lfmm/$file
mkdir -p $dir

## merge lfmm pvalues from 100 batches of 1 file into a temp file
date
echo "-> merge lfmm pvalues from 100 batches of file ${file} into a temp file\n"
lfmm_files=$(ls /shared/projects/most_kmer/afterkiss/gea/lfmm_rda/out/${file}_*/pvalue_bio_${var}.csv)
lfmm=$dir/pvalue_bio_${var}.csv
tmp=$dir/tmp_pvalue_bio_${var}.csv
cat $lfmm_files > $tmp
sort -k 1 $tmp > $lfmm
rm $tmp

## sort kmer sequences in the kmer position file
date
echo "-> sort kmer sequences in the kmer position file ${file}"
pos_file=/shared/projects/most_kmer/kiss_af/output/9.KMERPOSITION/output_file.${file}_vs_CC1.8_v2_pseudomolecule_cat_KMERPOSITION.txt
pos=$dir/output_file.${file}_vs_CC1.8_v2_pseudomolecule_cat_KMERPOSITION.txt
sort -k3 $pos_file > $pos

## map the pvalue from temp file to kmer position in the corresponding file
date
echo "-> map the pvalue from temp file to kmer position in the corresponding file ${file}"
map_lfmm=$dir/pvalue_bio_${var}_position.tsv

join -1 3 -2 1 --check-order $pos $lfmm -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 > $map_lfmm

## seperate the kmer/position/pvalue file by chromosome
date
echo "seperate the kmer/position/pvalue file by chromosome"
for chr in Contig CC1.8.Chr01 CC1.8.Chr02 CC1.8.Chr03 CC1.8.Chr04 CC1.8.Chr05 CC1.8.Chr06 CC1.8.Chr07 CC1.8.Chr08 CC1.8.Chr09 CC1.8.Chr10 CC1.8.Chr11
do
  echo "$chr"
  grep $chr $map_lfmm > $dir/pvalue_bio_${var}_$chr.tsv
done
