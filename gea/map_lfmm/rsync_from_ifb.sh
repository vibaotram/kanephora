#!/bin/bash
#SBATCH --job-name rsync
#SBATCH --mem=60G
#SBATCH --partition=highmem
#SBATCH --nodelist=node4
#SBATCH --error /scratch/most_kmer/map_lfmm/log/slurm-%x_%j.log
#SBATCH --output /scratch/most_kmer/map_lfmm/log/slurm-%x_%j.log

rsync -vaurP baotram@core.cluster.france-bioinformatique.fr:/shared/projects/most_kmer/afterkiss/gea/lfmm_rda /scratch/most_kmer
