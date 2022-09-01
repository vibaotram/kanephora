#!/bin/bash
#SBATCH --job-name gea_af
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --partition=normal,long,highmem
#SBATCH --error /data/projects/rob/scripts/kmer/gea/lfmm/log/slurm-%x_%a.log
#SBATCH --output /data/projects/rob/scripts/kmer/gea/lfmm/log/slurm-%x_%a.log
#SBATCH --array=[1-1728]

module load system/python/3.8.12
module load system/singularity/3.6.0

dist_dir=/data/projects/rob/scripts/kmer/gea/lfmm/out_K5

range=$SLURM_ARRAY_TASK_ID
if [[ $(($range%100)) == 0 ]]
then
r1=$(($range/100-1))
r2=100
else
r1=$(($range/100))
r2=$(($range%100))
fi
kmer_list=/data/projects/rob/scripts/kmer/gea/tmp/5.RANGES/output_file."$r1"/"$r2".txt
n=$r1"_"$r2

check_file=$dist_dir/$n/candidates_bio1.csv
if [[ -e $check_file ]]
then
echo "skip lfmm as the results already exists"
exit 0
else
echo "run lfmm"
fi

## create temp dir
mkdir -p /scratch/baotram
dir=$(mktemp -p /scratch/baotram -d --suffix="_$n")


echo "analyzing $kmer_list on $SLURM_JOB_NODELIST"
echo "working dir: $dir"


## copy files from tmpdir on nas
plink=$(echo $kmer_list | sed 's/\/[0-9]*.txt/.*/' | sed 's/5.RANGES/3.TABLE2BED/')
rsync -vauP $plink $dir/
rsync -vauP $kmer_list $dir/
rsync -vauP /data/projects/rob/scripts/kmer/gea/climatic_variables.csv $dir/
rsync -vauP /data/projects/rob/scripts/kmer/gea/lfmm/lfmm.R $dir/

bed=$(ls $dir/*.bed)
clim=$dir/climatic_variables.csv
out=$dir/$n

## run script

singularity exec /home/baotram/singularity-container_myr_4-0-2_rstudio_1.3.sif Rscript $dir/lfmm.R -b $bed -o $out -n $kmer_list -c $clim -K 5

## transfer output to nas

mkdir -p $dist_dir
rsync -vaurP $dir/$n $dist_dir

rm -rf $dir