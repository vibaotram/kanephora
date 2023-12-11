#!/bin/bash
#SBATCH --job-name map_lfmm_bio2
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --partition=highmem
#SBATCH --nodelist node4
#SBATCH --error /scratch/most_kmer/map_lfmm/log/slurm-%x_%a.log
#SBATCH --output /scratch/most_kmer/map_lfmm/log/slurm-%x_%a.log
#SBATCH --array=[1-1856]

v=2

module load system/python/3.8.12

range=$SLURM_ARRAY_TASK_ID
if [[ $(($range%100)) == 0 ]]
then
r1=$(($range/100-1))
r2=100
else
r1=$(($range/100))
r2=$(($range%100))
fi
kmer_list=$kiss_dir/5.RANGES/output_file."$r1"/"$r2".txt
n=$r1"_"$r2

pval_lfmm=/scratch/most_kmer/lfmm_rda/out/${n}/pvalue_bio_$v.csv
echo "pvalue file: $pval_lfmm"

pos_file=/scratch/most_kmer/map_lfmm/tmp_position1.txt
echo "position file: $pos_file"

output=/shared/projects/most_kmer/afterkiss/gea/map_lfmm/out/${n}/pvalue_bio_$v.tsv
echo "output files: ${output}_*"
mkdir -p $(dirname $output)

python /scratch/most_kmer/map_lfmm/map_lfmm.py --pvalue_file $pval_lfmm --kmers_position $pos_file --output $output