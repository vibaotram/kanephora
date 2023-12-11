#!/bin/bash
# SBATCH --job-name gengap_no45bi50_snp
#SBATCH --mem=80G
#SBATCH --partition=fast
#SBATCH --error /shared/projects/most_kmer/afterkiss/offset/log/slurm-%x_%a.log
#SBATCH --output /shared/projects/most_kmer/afterkiss/offset/log/slurm-%x_%a.log
# SBATCH --array=[1-640]


module load singularity


latfac_file=/shared/projects/most_kmer/afterkiss/gea/lfmm_rda/latent_factors_15_9.txt
matrix_dir=/shared/projects/most_kmer/afterkiss/gea/lfmm_rda
# climfile=/shared/projects/most_kmer/afterkiss/gea/explanatory_wc2.1.csv
climfile=/shared/projects/most_kmer/afterkiss/bioclim/wc1-4_present_Africa.csv

## run the genetic offset for vietnam climate
# new_env=/shared/projects/most_kmer/afterkiss/bioclim/wc2.1_30s_bioc_MRI-ESM2-0_ssp585_2041-2060_640occ_VN.csv
# l=$SLURM_ARRAY_TASK_ID
# # k=5
# name=$(basename $new_env | cut -d'_' -f4,5,6)
# outdir=/shared/projects/most_kmer/afterkiss/offset/${name}_18bc/${l}
# mkdir -p $outdir
# output=$(ls $outdir | wc -l)
# if [[ $output == 1 ]]; then echo "output exists"; exit; fi
# singularity exec -B /shared/projects/most_kmer /shared/projects/vietcaf/Singularity.R_4-3-Rserver_2022.sif \
# Rscript /shared/projects/most_kmer/afterkiss/offset/genetic_offset_extlatent.R \
# $latfac_file \
# $matrix_dir \
# $climfile \
# $new_env \
# $l \
# $outdir 


## run the genetic offset for african climate
new_env=/shared/projects/most_kmer/afterkiss/bioclim/$1"_Africa.csv"
# name=$(basename $new_env | cut -d'_' -f4,5,6)
name=$(basename $new_env | cut -d'_' -f1)
outdir=/shared/projects/most_kmer/afterkiss/offset/${name}
mkdir -p $outdir
output=$(ls $outdir | wc -l)
# if [[ $output == 1 ]]; then echo "output exists"; exit; fi
singularity exec -B /shared/projects/most_kmer /shared/projects/vietcaf/Singularity.R_4-3-Rserver_2022.sif \
Rscript /shared/projects/most_kmer/afterkiss/offset/genetic_offset_extlatent_africa.R \
$latfac_file \
$matrix_dir \
$climfile \
$new_env \
$outdir



# range=$SLURM_ARRAY_TASK_ID
# if [[ $(($range%100)) == 0 ]]
# then
# r1=$(($range/100-1))
# r2=100
# else
# r1=$(($range/100))
# r2=$(($range%100))
# fi
# b=${r1}_${r2}
# echo "batch $b" 

# k=5
# name=$(basename $new_env | cut -d'_' -f4,5,6)
# outdir=/shared/projects/most_kmer/afterkiss/offset/${name}_K${k}/${b}

# output=$(ls $outdir | wc -l)
# if [[ $output == 2 ]]; then echo "output exists"; exit; fi

# singularity exec -B /shared/projects/most_kmer /shared/projects/vietcaf/Singularity.R_4-3-Rserver_2022.sif \
# Rscript /shared/projects/most_kmer/afterkiss/offset/genetic_offset_average.R \
# $new_env \
# $b \
# $outdir \
# /shared/projects/most_kmer/afterkiss/gea/explanatory_wc2.1.csv 

# newclim_files=$(ls /shared/projects/most_kmer/afterkiss/bioclim/*bi50_Africa.csv)
# for newclim in $newclim_files
# do
#   id=$(basename $newclim | cut -d'_' -f1)
#   echo "sbatch --job-name gengap_${id}_snp /shared/ifbstor1/projects/most_kmer/afterkiss/offset/genetic_offset_ifb.sh $id"
# done