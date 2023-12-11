#!/bin/bash
#SBATCH --job-name gea_af_all
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --partition=fast
#SBATCH -A elai_most
#SBATCH --error /shared/projects/most_kmer/afterkiss/gea/lfmm_rda/log/slurm-%x_%a.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/gea/lfmm_rda/log/slurm-%x_%a.log
#SBATCH --array=[1-1856%35]

module load singularity

# dist_dir=/shared/projects/most_kmer/afterkiss/gea/lfmm/out_P2
dist_dir=/shared/projects/most_kmer/afterkiss/gea/lfmm_rda/out
# dist_dir=/shared/projects/most_kmer/afterkiss/gea/lfmm/out_K5_chelsa
kiss_dir=/shared/projects/most_kmer/kiss_af/output

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
outlier_list=$kiss_dir/6.PCADAPT/output_file."$r1"_"$r2"_BH0.05.pcadapt_outliers.csv


check_file=$dist_dir/$n/candidates_elevation.csv
if [[ -e $check_file ]]
then
echo "skip lfmm as the results already exists"
exit 0
else
echo "run lfmm"
fi

echo "analyzing $kmer_list on $SLURM_JOB_NODELIST"
echo "working dir: $dir"


bed="$kiss_dir/3.TABLE2BED/output_file.$r1.bed"
clim=/shared/projects/most_kmer/afterkiss/gea/explanatory_wc2.1.csv
# clim=$(dirname $dist_dir)/lfmm_explanatory_af_chelsa.csv
out=$dist_dir/$n

## run script

## run lfmm with all 19 variables and K = 5
singularity exec -B /shared/home/baotram \
/shared/projects/vietcaf/rserver/Singularity.R_4-2-0-Rserver_2022.sif \
Rscript $(dirname $dist_dir)/lfmm_rda.R -b $bed -o $out -n $kmer_list -c $clim -K 5 -s $outlier_list -r 2 -f 0.01


## run lfmm with first 2 PCs and K = 5
# singularity exec -B /shared/home/baotram \
# /shared/projects/vietcaf/rserver/Singularity.R_4-2-0-Rserver_2022.sif \
# Rscript $(dirname $dist_dir)/lfmm.R -b $bed -o $out -n $kmer_list -c $clim -K 5 -p 2

