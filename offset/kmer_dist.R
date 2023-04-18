## calculate distance between individuals
## using presence/absence table of significant kmers
## input:
#### 1- number of total kmers (100489652)
#### 2- number of kmers each batch (10000)
#### 3- batch index ($SLURM_ARRAY_TASK_ID)
#### 4- outdir
#### 5- number of cores($SLURM_CPUS_PER_TASK)
## output:
#### a matrix containing number of different kmers between individuals

setwd("/shared/projects/most_kmer/afterkiss/offset")

library(tidyverse)
library(data.table)

args = commandArgs(trailingOnly = TRUE)
cat("arguments: ", args, "\n")

cat("read kmer table\n")
# kmer_table <- fread("./signif_kmer_table.txt", drop = 1, data.table = F, nThread = 30)
kmer_inds <- fread("./signif_kmer_table.txt", drop = 1, data.table = F, nThread = 30, nrows = 0)


tk <- as.integer(args[1])
nk <- as.integer(args[2])
batch <- as.integer(args[3])
ls <- seq((batch - 1)*nk + 1, min(batch*nk, tk))

cat("from line", ls[1], "to", max(ls), "\n")


kmer_sub <- fread("./signif_kmer_table.txt", drop = 1, data.table = F, nThread = 30,
                  skip = ls[1], 
                  nrows = length(ls))
colnames(kmer_sub) <- colnames(kmer_inds)
# kmer_sub <- kmer_table[which(lv == levels(lv)[i]),]
cat("count kmer difference pairwise\n")
kmer_sub <- t(kmer_sub)
dist_sub <- dist(kmer_sub, method = "manhattan", upper = T)

saveRDS(dist_sub, file = file.path(args[4], paste0("dist_sub_", batch, ".RDS")))