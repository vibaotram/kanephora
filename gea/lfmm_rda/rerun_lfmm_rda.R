#!/usr/bin/env Rscript --vanilla
#
#########
######### input (from kiss): bed/bim files (in 3.TABLE2BED), a list of tested kmers (in 5.RANGES), a climatic data file (id in first column), K (number of groups), an output dir
######### output: a list of candidate kmers, a histogram of calibrated p-values, GIF, saved them to the outdir
######### 1, run pca on the climatic data (optional)
######### 2, run lfmm with the genetic matrix and clim variables in 2 ways: using PCA of clim variables, or all the variables

.libPaths("/shared/home/baotram/R/x86_64-pc-linux-gnu-library/4.0")


## load packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(lfmm))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(BEDMatrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
# suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pcadapt))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(robust))
suppressPackageStartupMessages(library(vegan))

option_list <- list(
  make_option(c("-b", "--bedFile"),
              type = "character",
              default = NULL,
              help = "Path of bed file"),
  make_option(c("-o", "--outdir"),
              type = "character",
              default = NULL,
              help = "Path where output sequence file will be written."),
  make_option(c("-n", "--kmerList"),
              type = "character",
              default = NULL,
              help = "path to kmer list"),
  make_option(c("-c", "--climFile"),
              type = "character",
              default = NULL,
              help = "climatic variable file"),
  # make_option(c("-p", "--pca"),
  #             type = "numeric",
  #             default = NULL,
  #             help = "number of pc of climatic variables to use in lfmm"),
  # make_option(c("-v", "--variable"),
  #             type = "numeric",
  #             default = NULL,
  #             help = "variable indexes of climatic variables to use in lfmm"),
  make_option(c("-K", "--numgroup"),
              type = "numeric",
              default = NULL,
              help = "number of ancestral groups used as latent factor in lfmm"),
  make_option(c("-r", "--retain"),
              type = "numeric",
              default = NULL,
              help = "number of rda axis to retain"), 
  make_option(c("-f", "--fdr"),
              type = "numeric",
              default = NULL,
              help = "false discovery rate to select candidate kmers")
)



myArgs <- parse_args(
  OptionParser(usage = "%prog [options]", option_list = option_list,
               description = "run lfmm using genotype in bed file and climatic variables in the clim file"))

cat("Reading input bed file\n")
bedfile <- myArgs$bedFile
bed <- BEDMatrix(bedfile, simple_names = T)
outdir <- myArgs$outdir
dir.create(outdir, showWarnings = F, recursive = T)
cat("Reading input subset kmers\n")
kmerfile <- myArgs$kmerList
kmer <- scan(kmerfile)

r <- as.numeric(myArgs$retain)
f <- as.numeric(myArgs$fdr)



cat("Reading input clim file\n")
clim <- read.csv(myArgs$climFile)
clim <- clim %>% 
  select(Label, bio_1:bio_19) %>% 
  column_to_rownames(colnames(clim)[1])
clim <- clim[, !colnames(clim) %in% c("bio_10", "bio_1", "bio_3", "bio_14", "bio_9", "bio_16")]
# clim_pca <- rda(clim, scale = T)
# clim_pc <- scores(clim_pca, choices=1:2, display="sites", scaling=0)
# pc <- myArgs$pca
# vi <- myArgs$variable
K <- myArgs$numgroup


#######   #######   #######

# # detect outlier k-mers
# kmermat <- bed[, sort(kmer)]
# kmerpca <- read.pcadapt(t(kmermat), type = "pcadapt")
# pcadapt <- pcadapt(input = kmerpca, K = 2)
# 
# padj <- p.adjust(pcadapt$pvalues,method="BH")
# alpha <- 0.05
# outliersBH0.05 <- which(padj < alpha)
# snp_pc <- get.pc(pcadapt, outliersBH0.05)


## remove individuals missing in the climfile
## keep kmers in the kmerlist
## 18.10.2022 - keep outlier kmers
# submat <- bed[which(rownames(bed) %in% rownames(clim)), outliersBH0.05]
cat("Extracting bed matrix for subset kmers\n")
submat <- bed[which(rownames(bed) %in% rownames(clim)), sort(kmer)]
submat <- submat/2
kmer_count <- colSums(submat)
submat <- submat[, which(kmer_count < nrow(submat) & kmer_count > 0)]

# explanatory variables: pca of climatic variables if pc is set numeric, else selected variables, else all the variables

# if (is.numeric(pc)) {
#   print(paste0("run lfmm using ", pc, " PCs from PCA of climatic variables"))
#   clim_pca <- rda(clim, scale = T)
#   expl <- scores(clim_pca, choices = c(1:pc), display = "sites", scaling = 0)
# } else if (is.numeric(vi)) {
#   print("run lfmm with climatic variables: ")
#   print(colnames(clim)[vi])
#   expl <- clim[,vi]
# } else {
#   print("run lfmm with all the provided climatic data\n")
#   expl <- clim
# }


# lfmm
cat("Running LFMM\n")
sub_lfmm <- lfmm_ridge(Y = submat, X = clim, K = K)
pv_lfmm <- lfmm_test(Y = submat, X = clim, 
                lfmm = sub_lfmm, 
                calibrate = "gif")

signif_kmers_lfmm <- unique(unlist(lapply(1:ncol(pv_lfmm$calibrated.pvalue), function(c) {
  # filter out NaN pvalues
  pv_filt <- na.omit(as.data.frame(pv_lfmm$calibrated.pvalue[,c]), cols = 1)
  colnames(pv_filt) <- "pvalue"
  
  # fdr control
  qv <- qvalue(pv_filt$pvalue, fdr.level = f, pfdr = T)
  fdr_cands <- which(qv$qvalues < f)
  fdr_cands
})))

top_kmers_lfmm <- unique(unlist(lapply(1:ncol(pv_lfmm$calibrated.pvalue), function(c) {
  # filter out NaN pvalues
  pv_filt <- na.omit(as.data.frame(pv_lfmm$calibrated.pvalue[,c]), cols = 1)
  colnames(pv_filt) <- "pvalue"
  
  # fdr control
  qv <- qvalue(pv_filt$pvalue, fdr.level = f, pfdr = T)
  fdr_cands <- which(qv$qvalues < f)
  q10 <- quantile(qv$qvalues[fdr_cands], probs = 0.1)
  top_cands <- which(!qv$qvalues > q10)
  top_cands
})))

# pvalues_lfmm <- qvalue(pv_lfmm$calibrated.pvalue, fdr.level = 0.01, pfdr = T)
# signif_kmers_lfmm <- which(pvalues_lfmm$qvalues < 0.01)


cat("Writing output of LFMM\n")
# lfmm_res <- data.frame(kmer = colnames(submat),
#                        pvalue = as.numeric(pvalues_lfmm$pvalues),
#                        qvalue = as.numeric(pvalues_lfmm$qvalues))
write.table(pv_lfmm$calibrated.pvalue, file.path(outdir, "lfmm_results.tsv"), 
            quote = FALSE, sep = "\t", row.names = T, col.names = T)


## rda
cat("Running RDA\n")
# conf_fac <- sub_lfmm$U
# rownames(conf_fac) <- rownames(submat)
# colnames(conf_fac) <- paste("V", 1:K, sep = "")
# env <- cbind(clim, conf_fac)
# sub_rda <- rda(submat ~ bio_1:bio_19 + Condition(V1:V5), data = env, scale = T) #  + Condition(PC1 + PC2)
# pdf(file.path(outdir, "rda_results.pdf"))
# screeplot(sub_rda)
# plot(sub_rda, scaling = 1)

sub_rda <- rda(submat ~ ., clim, scale = F)
## select candidate kmers
rdadapt <- function(rda, K) {
  zscores <- rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action = na.omit, estim = "pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df = K)
  reschi2test <- pchisq(resmaha/lambda,K, lower.tail = FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt <- qval$qvalues
  return(data.frame(p.values = reschi2test, q.values = q.values_rdadapt))
}

signif_rdadapt <- rdadapt(sub_rda, r)
# hist(signif_rdadapt$p.values)
stopifnot(nrow(signif_rdadapt) == ncol(submat))
signif_kmers_rda <- which(signif_rdadapt$q.values < f)

rda_q10 <- quantile(signif_rdadapt$q.values[signif_kmers_rda], probs = 0.1)
top_kmers_rda <- which(!signif_rdadapt$q.values > rda_q10)

# load.rda <- scores(sub_rda, choices=c(1:3), display="species")
# 
# outliers <- function(x,z){
#   lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
#   which(x < lims[1] | x > lims[2])               # locus names in these tails
# }
# rda_cand1 <- outliers(load.rda[,1],3) # 38
# rda_cand2 <- outliers(load.rda[,2],3) # 69
# rda_cand3 <- outliers(load.rda[,3],3)
# 
# signif_kmers_rda <- unique(c(rda_cand1, rda_cand2, rda_cand3))

cat("Writing output of RDA\n")
rda_res <- data.frame(kmer = colnames(submat),
                      pvalue = signif_rdadapt$p.values,
                      qvalue = signif_rdadapt$q.values)
write.table(rda_res, file.path(outdir, "rda_results.tsv"),
            quote = FALSE, sep = "\t", row.names = F, col.names = T)

cat("Number of candidate k-mers detected by methods:\n")
cat("RDA", length(signif_kmers_rda), "\n")

# signif_lfmm <- fread(myArgs$lfmm)
cat("LFMM", length(signif_kmers_lfmm), "\n")

signif_kmers <- intersect(signif_kmers_rda, signif_kmers_lfmm)
# signif_data <- data.frame(kmer = colnames(submat)[signif_kmers],
#                           pvalue_lfmm = pvalues_lfmm$pvalues[signif_kmers],
#                           pvalue_rda = signif_rdadapt$p.values[signif_kmers])

top_kmers <- intersect(top_kmers_lfmm, top_kmers_rda)

cat("RDA+LFMM", length(signif_kmers), sprintf("%.2f", length(signif_kmers)/ncol(submat)),"\n")

cat("Writing final list of significant kmers\n")
# write.table(signif_data, file.path(outdir, "lfmm_rda.tsv"), 
 #            quote = FALSE, sep = "\t", row.names = F, col.names = T)
writeLines(colnames(submat)[signif_kmers], file.path(outdir, "candidates_rda_lfmm.txt"))
writeLines(colnames(submat)[signif_kmers_lfmm], file.path(outdir, "candidates_lfmm.txt"))
writeLines(colnames(submat)[signif_kmers_rda], file.path(outdir, "candidates_rda.txt"))
writeLines(colnames(submat)[top_kmers], file.path(outdir, "top_kmer_rda_lfmm.txt"))

cand_mat <- submat[, signif_kmers]
write.table(cand_mat, file.path(outdir, "candidates_matrix.txt"),
            col.names = T, row.names = T, sep = "\t", quote = F)

# tot_cands <- pv$calibrated.pvalue[sapply(1:nrow(pv$calibrated.pvalue), function(i) any(pv$calibrated.pvalue[i,] < 0.05)),]
# cat("Number of total candidate k-mers:", nrow(tot_cands), "\n")


# bedfile <- "./gea/tmp/3.TABLE2BED/output_file.0.bed"
# kmerfile <- "./gea/tmp/5.RANGES/output_file.0/1.txt"

# bedfile <- "/shared/projects/most_kmer/kiss_af/output/3.TABLE2BED/output_file.0.bed"
# kmerfile <- "/shared/projects/most_kmer/kiss_af/output/5.RANGES/output_file.0/1.txt"
# climFile <- "/shared/projects/most_kmer/afterkiss/bioclim/updated_bioclim_AF.csv"
# clim <- fread(climFile)
# selectivefile <- "/shared/projects/most_kmer/kiss_af/output/6.PCADAPT/output_file.0_1_BH0.05.pcadapt_outliers.csv"
