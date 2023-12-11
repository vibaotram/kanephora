#!/usr/bin/env Rscript --vanilla

## estimate genetic offset using all 18M significant kmers
## when moving the african inds to Vietnam in the present/future
## input: 1- new bioclim data file containing 640 points, 2- line (point) to test, 3- output dir, 
## 4- path containing genotype matrix files, 5- path of the original bioclim data, 6- number of latent factors
## output: 1- genetic offet, rona, climate distance for each individual, 2- contribution of each bioclim variable


# .libPaths("/home/baotram/R/x86_64-pc-linux-gnu-library/4.3")
.libPaths("/shared/home/baotram/R/x86_64-pc-linux-gnu-library/4.3")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(LEA))

# setwd("/scratch/most_kmer")
setwd("/shared/projects/most_kmer/afterkiss")

## input
arg <- commandArgs(trailingOnly = T)
new_env_file <- arg[1]
l <- as.numeric(arg[2])
outdir <- arg[3]
matrix_dir <- arg[4]
clim_file <- arg[5]
k <- as.numeric(arg[6])
# dir.create(outdir, showWarnings = F)

## read genotype
cat("read genotype\n")
# mat_files <- list.files("./lfmm_rda", "matrix_significant", full.names = T)
mat_files <- list.files(matrix_dir, "matrix_significant", full.names = T)
signif_mat <- do.call(cbind, sapply(mat_files, function(f) {
  m <- fread(f, data.table = F)
  m <- m %>% column_to_rownames(colnames(m)[1])
  }, simplify = T, USE.NAMES = F))

## read present bioclim data
cat("read present bioclim data\n")
# clim_file <- "./bioclim/explanatory_wc2.1.csv"
# clim_file <- "./gea/explanatory_wc2.1.csv"
clim <- fread(clim_file, data.table = F)
clim <- clim %>% column_to_rownames(colnames(clim)[1])
pres_clim <- matrix(unlist(clim[,1:19]), nrow = nrow(clim))

## read new bioclim data
cat("read new bioclim data\n")
new_clim <- fread(new_env_file, data.table = F)
new_clim <- new_clim %>% 
  relocate(str_sort(colnames(new_clim), numeric = T)) 
# colnames(fut_clim) <- gsub("bio", "bio_", colnames(fut_clim))
pred_clim <- matrix(rep(as.numeric(new_clim[l,1:19]), nrow(clim)), nrow = nrow(clim), byrow = T)

## estimate genetic gap
cat("estimate genetic gap\n")
gap <- genetic.gap(input = signif_mat,
                     env = pres_clim,
                     pred.env = pred_clim,
                     scale = TRUE,
                     K = k)


## estimate climate euclidean distance
env_dist <- sqrt(rowSums((pred_clim - pres_clim)^2))

cat("output results\n")
## save genetic offet, rona, climate distance for each individual
gap_ind <- data.frame(offset = gap$offset, 
                      rona = gap$distance,
                      env_dist = env_dist)
rownames(gap_ind) <- rownames(clim)
write.table(gap_ind, file.path(outdir, "genetic_gap_ind.tsv"), 
            row.names = T, col.names = T, quote = F)

## save the contribution of each bioclim variable
gap_env <- cbind(gap$vectors, eigen = gap$eigenvalues)
colnames(gap_env)[1:19] <- colnames(new_clim)[1:19]
rownames(gap_env) <- colnames(new_clim)[1:19]
write.table(gap_env, file.path(outdir, "genetic_gap_env.tsv"), 
            row.names = T, col.names = T, quote = F)
