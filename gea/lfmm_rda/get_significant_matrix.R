#!/usr/bin/env Rscript --vanilla

## get presence/absence matrix for significant kmers
## from a bed file
## input: bed file
## output: a matrix in tsv file

.libPaths("/shared/home/baotram/R/x86_64-pc-linux-gnu-library/4.0")
library(BEDMatrix)
library(data.table)
library(tidyverse)

setwd("/shared/projects/most_kmer/afterkiss")

## bed file
args <-  commandArgs(trailingOnly=TRUE)

bed_file <- args[1]
bedmat <- BEDMatrix(bed_file, simple_names = T)
batch <- gsub("output_file\\.|\\.bed", "", basename(bed_file))


## get significant kmers
com_files <- list.files("./gea/lfmm_rda/out", "candidates_rda_lfmm.txt", 
                        recursive = T, full.names = T)
com_batch <- grep(paste0("./gea/lfmm_rda/out/", batch, "_"), com_files, value = T)

kmer_seq <- unlist(sapply(com_batch, readLines, 
                          simplify = T, USE.NAMES = F))

clim_file <- "./gea/explanatory_wc2.1.csv"
clim <- fread(clim_file, data.table = F)
clim <- clim %>% column_to_rownames(colnames(clim)[1])

## extract matrix for the significant kmers
com_mat <- bedmat[rownames(bedmat) %in% rownames(clim), kmer_seq]
fwrite(as.data.frame(com_mat), paste0("./gea/lfmm_rda/matrix_significant_kmers_output_file_", batch, ".tsv"),
       row.names = T, sep = "\t")
