#!/usr/bin/env Rscript --vanilla

## estimate genetic offset in each kmer subset, then take the average
## when moving the african inds to Vietnam in the present/future
## input: 1- new bioclim data file containing 640 points, 2- kmer batch, 3- output dir, 4- path to present clim matrix
## output: 1- genetic offet, rona, climate distance for each individual, 2- contribution of each bioclim variable


# .libPaths("/home/baotram/R/x86_64-pc-linux-gnu-library/4.3")
.libPaths("/shared/home/baotram/R/x86_64-pc-linux-gnu-library/4.3")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(LEA))
suppressPackageStartupMessages(library(BEDMatrix))

# setwd("/scratch/most_kmer")
setwd("/shared/projects/most_kmer/afterkiss")

## input
arg <- commandArgs(trailingOnly = T)
new_env_file <- arg[1]
b <- arg[2]
outdir <- arg[3]
# matrix_dir <- arg[4]
clim_file <- arg[4]
# k <- as.numeric(arg[6])
dir.create(outdir, showWarnings = F)


## read present bioclim data
cat("read present bioclim data\n")
# clim_file <- "./bioclim/explanatory_wc2.1.csv"
# clim_file <- "./gea/explanatory_wc2.1.csv"
clim <- fread(clim_file, data.table = F)
clim <- clim %>% column_to_rownames(colnames(clim)[1])
pres_clim <- matrix(unlist(clim[,1:19]), nrow = nrow(clim))

## read new bioclim data
cat("read new bioclim data\n")
# new_env_file <- "./bioclim/wc2.1_30s_bioc_present_1970-2000_640occ_VN.csv"
new_clim <- fread(new_env_file, data.table = F)
new_clim <- new_clim %>% 
  relocate(str_sort(colnames(new_clim), numeric = T)) 
# colnames(fut_clim) <- gsub("bio", "bio_", colnames(fut_clim))
# pred_clim <- matrix(rep(as.numeric(new_clim[l,1:19]), nrow(clim)), nrow = nrow(clim), byrow = T)


## read genotype
cat("read genotype\n")
lfmm_out <- "./gea/lfmm_rda/out"
submat_file <- file.path(lfmm_out, b, "kmer_submat.tsv")
candfile <- file.path(lfmm_out, b, "candidates_rda_lfmm.txt")

submat <- fread(submat_file)
submat <- submat %>% 
  column_to_rownames("V1")

candkmer <- readLines(candfile)

## genetic offset
cat("genetic offset\n")
for (l in 1:nrow(new_clim)) {
  output <- file.path(outdir, paste0("genetic_gap_env_", l, ".tsv"))
  if (file.exists(output)) next
  cat(l, "| ")
  pred_clim <- matrix(rep(as.numeric(new_clim[l,1:19]), nrow(clim)), nrow = nrow(clim), byrow = T)
  gap <- genetic.gap(input = submat,
                     env = pres_clim,
                     pred.env = pred_clim,
                     scale = TRUE,
                     K = 5,
                     candidate.loci = which(colnames(submat) %in% candkmer))
  gap_batch <- data.frame(ind = rownames(clim),
                          offset = gap$offset,
                          rona = gap$distance,
                          batch = b,
                          env_dist = sqrt(rowSums((pred_clim - pres_clim)^2)))
  write.table(gap_batch, output,
              row.names = F, col.names = T, quote = F)
}


# gap_ind <- data.frame()


# batch <- list.files(lfmm_out)

# for (b in c(0:18)) {
#   cat("## bed file:", b, "\n")
#   bedfile <- paste0("/shared/projects/most_kmer/kiss_af/output/3.TABLE2BED/output_file.", b, ".bed")
#   mat <- BEDMatrix(bedfile, simple_names = T)
#   nb_sets <- length(list.files(paste0("/shared/projects/most_kmer/kiss_af/output/5.RANGES/output_file.", b), ".txt"))
#   cat("subset:")
#   for (s in c(1:nb_sets)) {
#     ## kmer subset
#     cat(s, "| ")
#     rangefile <- paste0("/shared/projects/most_kmer/kiss_af/output/5.RANGES/output_file.", b, "/", s, ".txt")
#     kmer <- scan(rangefile, quiet = T)
#     submat <- mat[which(rownames(mat) %in% rownames(clim)), sort(kmer)]
#     kmer_count <- colSums(submat)
#     submat <- submat[, which(kmer_count < nrow(submat)*2 & kmer_count > 0)]
#     colnames(submat) <- gsub("_\\d+", "", colnames(submat))
#     ## significant kmer
#     candfile <- paste0("/shared/projects/most_kmer/afterkiss/gea/lfmm_rda/out/", b, "_", s, "/candidates_rda_lfmm.txt")
#     candkmer <- readLines(candfile)
#     ## genetic offset
#     gap <- genetic.gap(input = submat,
#                        env = pres_clim,
#                        pred.env = pred_clim,
#                        scale = TRUE,
#                        K = 5,
#                        candidate.loci = which(colnames(submat) %in% candkmer))
#     gap_batch <- data.frame(ind = rownames(clim),
#                             offset = gap$offset, 
#                             rona = gap$distance,
#                             batch = paste0(b, "_", s),
#                             env_dist = sqrt(rowSums((pred_clim - pres_clim)^2)))
#     gap_ind <- rbind(gap_ind, gap_batch)
#   }
# }

# gap_ind <- gap_ind %>% 
#   group_by(ind) %>% 
#   summarise(mean_go = mean(offset),
#             mean_rona = mean(rona))
# gap_ind$env_dist <- sqrt(rowSums((pred_clim - pres_clim)^2))
