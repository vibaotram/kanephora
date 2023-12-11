#!/bin/bash
#SBATCH --job-name gengap_present
#SBATCH --mem=80G
#SBATCH --partition=highmemdell
#SBATCH --nodelist=node4
#SBATCH --error /scratch/most_kmer/offset/log/slurm-%x_%a.log
#SBATCH --output /scratch/most_kmer/offset/log/slurm-%x_%a.log
#SBATCH --array=[1-640]

range=$SLURM_ARRAY_TASK_ID
if [[ $(($range%100)) == 0 ]]
then
r1=$(($range/100-1))
r2=100
else
r1=$(($range/100))
r2=$(($range%100))
fi
b=${r1}_${r2}
echo "batch $b"

module load system/singularity/3.6.0

new_env=/scratch/most_kmer/bioclim/wc2.1_30s_bioc_present_1970-2000_600districts_VN.csv
# l=$SLURM_ARRAY_TASK_ID
k=5
name=$(basename $new_env | cut -d'_' -f4,5,6)
outdir=/scratch/most_kmer/offset/${name}_K${k}/${b}
mkdir -p $outdir

singularity exec -B /scratch/most_kmer ~/Singularity.R_4-3-Rserver_2022.sif \
Rscript /scratch/most_kmer/offset/genetic_offset_average.R \
$new_env \
$b \
$outdir \
/scratch/most_kmer/gea/explanatory_wc2.1.csv 
