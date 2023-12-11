#!/usr/bin/env Rscript --vanilla

## estimate genetic offset using all 18M significant kmers
## when moving the african inds to Vietnam in the present/future
## input: 1- file containing matrix of latent factors, 2- directory containing kmer matrices, 3- file containing original bioclim, 
## 4- file containing new bioclim, 5- occurrence in new bioclim data (index), 6- output dir
## output: 1- genetic offet, rona, climate distance for each individual, 2- contribution of each bioclim variable


# .libPaths("/home/baotram/R/x86_64-pc-linux-gnu-library/4.3")
.libPaths("/shared/home/baotram/R/x86_64-pc-linux-gnu-library/4.3")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(LEA))


arg <- commandArgs(trailingOnly = T)
latfac_file <- arg[1]
matrix_dir <- arg[2]
clim_file <- arg[3]
new_env_file <- arg[4]
outdir <- arg[5]



# latfac_file <- "/shared/projects/most_kmer/afterkiss/gea/lfmm_rda/latent_factors_15_9.txt"
# matrix_dir <- "/shared/projects/most_kmer/afterkiss/gea/lfmm_rda"
# clim_file <- "/shared/projects/most_kmer/afterkiss/gea/explanatory_wc2.1.csv"
# l <- 1
# new_env_file <- "/shared/projects/most_kmer/afterkiss/bioclim/wc2.1_30s_bioc_EC-Earth3-Veg_ssp585_2041-2060_640occ_VN.csv"

go_ex <- function(Y, X, K = 0, X.new) {
  mod.lm = lm(Y ~ X)
  B = t(mod.lm$coefficients[-1,])
  # B[is.na(B)] <- 0
  genetic.gap = rowSums(((X.new - X) %*% t(B))^2)/nrow(B) 
  rona = rowSums(abs((X.new - X) %*% t(B)))/nrow(B)
  return(list(offset = genetic.gap, distance = rona))
}


## use only common inds with Remi
com_inds <- fread("/shared/projects/most_kmer/afterkiss/offset/IPSL-CM6A-LR_ssp245_2041-2060/Genomic_Vulnerability_rcp4.5.txt",
                  data.table = F)
com_inds <- na.omit(com_inds)

## Y
cat("read kmer matrix\n")
mat_files <- list.files(matrix_dir, "matrix_significant", full.names = T)
signif_mat <- do.call(cbind, sapply(mat_files, function(f) {
  m <- fread(f, data.table = F)
  m <- m %>% column_to_rownames(colnames(m)[1]) ## get ind label for rownames
  m <- m[rownames(m) %in% com_inds$Label,] ## remove uncommon inds
}, simplify = T, USE.NAMES = F))
signif_mat <- as.matrix(signif_mat)



## retained bioclim
bc <- c("bio_18", "bio_15", "bio_12", "bio_10", "bio_2")

## X
cat("read present bioclim\n")
clim <- fread(clim_file, data.table = F)
clim <- clim %>% column_to_rownames(colnames(clim)[1]) ## get ind label for rownames
colnames(clim) <- gsub("bio", "bio_", colnames(clim)) ## rename bioclim variables
clim <- clim[rownames(clim) %in% com_inds$Label,] ## remove uncommon inds
clim <- clim[, bc]
pres_clim <- matrix(unlist(clim), nrow = nrow(clim))


# ## U
# cat("read latent factors\n")
# latfac <- read.table(latfac_file)
# latfac <- as.matrix(latfac)
# 
# 
# ## K
# k <- ncol(latfac)
k <- 0

## X.new
cat("read new bioclim\n")
new_clim <- fread(new_env_file, data.table = F)
new_clim <- new_clim %>% column_to_rownames(colnames(new_clim)[1]) ## get ind label for rownames
# colnames(new_clim) <- gsub("c_.+_", "_", colnames(new_clim)) ## get ind label for rownames
colnames(new_clim) <- gsub("^.+bi50", "bio_", colnames(new_clim)) ## rename bioclim variables
new_clim <- new_clim[rownames(new_clim) %in% com_inds$Label,] ## remove uncommon inds
new_clim <- new_clim %>% 
  relocate(str_sort(colnames(new_clim), numeric = T)) 
# colnames(fut_clim) <- gsub("bio", "bio_", colnames(fut_clim))
new_clim <- as.matrix(new_clim[,bc])


## genetic offset using external latent factors
cat("run genetic offset\n")
gap <- go_ex(Y = as.matrix(signif_mat),
             X = pres_clim, 
             # U = as.matrix(latfac),
             K = k,
             X.new = new_clim)



## estimate climate euclidean distance
env_dist <- sqrt(rowSums((new_clim - pres_clim)^2))

cat("output results\n")
## save genetic offet, rona, climate distance for each individual
gap_ind <- data.frame(offset = gap$offset, 
                      rona = gap$distance,
                      env_dist = env_dist)
rownames(gap_ind) <- rownames(clim)
write.table(gap_ind, file.path(outdir, "genetic_gap_ind.tsv"), 
            row.names = T, col.names = T, quote = F)