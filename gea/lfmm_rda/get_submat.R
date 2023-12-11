#!/usr/bin/env Rscript --vanilla

## estimate genetic offset in each kmer subset, then take the average
## when moving the african inds to Vietnam in the present/future
## input: 1- bedfile, 2- path to present clim matrix, 3- output dir 
## output: 1- genetic offet, rona, climate distance for each individual, 2- contribution of each bioclim variable


# .libPaths("/home/baotram/R/x86_64-pc-linux-gnu-library/4.3")
.libPaths("/shared/home/baotram/R/x86_64-pc-linux-gnu-library/4.3")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BEDMatrix))

# setwd("/scratch/most_kmer")
setwd("/shared/projects/most_kmer/afterkiss")

## input
arg <- commandArgs(trailingOnly = T)
bedfile <- arg[1]
clim_file <- arg[2]
outdir <- arg[3]
# matrix_dir <- arg[4]
# k <- as.numeric(arg[6])
dir.create(outdir, showWarnings = F)


## read present bioclim data
cat("read present bioclim data\n")
# clim_file <- "./bioclim/explanatory_wc2.1.csv"
# clim_file <- "./gea/explanatory_wc2.1.csv"
clim <- fread(clim_file, data.table = F)
clim <- clim %>% column_to_rownames(colnames(clim)[1])

# bedfile <- paste0("/shared/projects/most_kmer/kiss_af/output/3.TABLE2BED/output_file.", b, ".bed")
mat <- BEDMatrix(bedfile, simple_names = T)
b <- gsub("output_file.|.bed", "", basename(bedfile))
nb_sets <- length(list.files(paste0("/shared/projects/most_kmer/kiss_af/output/5.RANGES/output_file.", b), ".txt"))
cat("subset:")
for (s in c(1:nb_sets)) {
  ## kmer subset
  cat(s, "| ")
  rangefile <- paste0("/shared/projects/most_kmer/kiss_af/output/5.RANGES/output_file.", b, "/", s, ".txt")
  kmer <- scan(rangefile, quiet = T)
  submat <- mat[which(rownames(mat) %in% rownames(clim)), sort(kmer)]
  kmer_count <- colSums(submat)
  submat <- submat[, which(kmer_count < nrow(submat)*2 & kmer_count > 0)]
  colnames(submat) <- gsub("_\\d+", "", colnames(submat))
  write.table(submat, file = file.path(outdir, paste0(b, "_", s), "kmer_submat.tsv"),
              row.names = T, col.names = T, quote = F)
}

