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
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pcadapt))
suppressPackageStartupMessages(library(qvalue))

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
  make_option(c("-p", "--pca"),
              type = "numeric",
              default = NULL,
              help = "number of pc of climatic variables to use in lfmm"),
  make_option(c("-v", "--variable"),
              type = "numeric",
              default = NULL,
              help = "variable indexes of climatic variables to use in lfmm"),
  make_option(c("-K", "--numgroup"),
              type = "numeric",
              default = NULL,
              help = "number of ancestral groups used as latent factor in lfmm"),
  make_option(c("-s", "--selectiveKmers"),
              type = "character",
              default = NULL,
              help = "file containing selective/outlier kmers"),
  make_option(c("-g", "--backgroundKmers"),
              type = "character",
              action = NULL,
              help = "file containing background kmers in lfmm")
)



myArgs <- parse_args(
  OptionParser(usage = "%prog [options]", option_list = option_list,
               description = "run lfmm using genotype in bed file and climatic variables in the clim file"))

bedfile <- myArgs$bedFile
bed <- BEDMatrix(bedfile, simple_names = T)
outdir <- myArgs$outdir
dir.create(outdir, showWarnings = F, recursive = T)
kmerfile <- myArgs$kmerList
selectivefile <- myArgs$selectiveKmers
outliers <- read.table(selectivefile, sep = "\t")
selective_kmers <- outliers[,1]
selective_kmers <- gsub("_0", "", selective_kmers)

bg_kmer_file <- myArgs$backgroundKmers

if (kmerfile != selectivefile) {
  cat("run lfmm on all k-mers\n")
  kmer <- scan(kmerfile) ### NOTE: the object `kmer` here contains the index of the kmers out of the bed file
} else {
  cat("run lfmm on outlier k-mers\n")
  outkmer <- outliers[,4] ### NOTE: the object `outkmer` here contains the kmer indexes out of the bed file
  
  if (!is.null(bg_kmer_file)) {
    cat("add", length(outkmer), "random background k-mers to lfmm analysis\n")
    allkmer <- scan(bg_kmer_file) ### NOTE: the object `allkmer` here contains the kmer indexes out of the bed file
    bgkmer <- sample(allkmer[!allkmer %in% outkmer], length(outkmer), replace = F)
    kmer <- c(outkmer, bgkmer)
  } else {
    kmer <- outkmer
  }
  
}


clim <- fread(myArgs$climFile)
clim <- clim %>% column_to_rownames(colnames(clim)[1])
pc <- myArgs$pca
vi <- myArgs$variable
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

submat <- bed[which(rownames(bed) %in% rownames(clim)), sort(kmer)]


# explanatory variables: pca of climatic variables if pc is set numeric, else selected variables, else all the variables

if (is.numeric(pc)) {
  print(paste0("run lfmm using ", pc, " PCs from PCA of climatic variables"))
  clim_pca <- rda(clim, scale = T)
  expl <- scores(clim_pca, choices = c(1:pc), display = "sites", scaling = 0)
} else if (is.numeric(vi)) {
  print("run lfmm with climatic variables: ")
  print(colnames(clim)[vi])
  expl <- clim[,vi]
} else {
  print("run lfmm with all the provided climatic data\n")
  expl <- clim
}


# lfmm
sub_lfmm <- lfmm_ridge(Y = submat, X = expl, K = K)
pv <- lfmm_test(Y = submat, X = expl, 
                     lfmm = sub_lfmm, 
                     calibrate = "gif")



hist_plot <- file.path(outdir, "pvalues_histogram.pdf")
pdf(file = hist_plot)
for (i in 1:ncol(expl)) {
  pvalues <- pv$calibrated.pvalue[,i]
  hist(pvalues,
       main = paste0("Calibrated p-values for ", colnames(expl)[i],
                     " (GIF = ", sprintf("%.2f", pv$gif[i]), ")"),
       xlab = "p-values")
}
dev.off()


qq_plot <- file.path(outdir, "pvalues_qqplot.pdf")
pdf(file = qq_plot)
for (i in 1:ncol(expl)) {
  pvalues <- pv$calibrated.pvalue[,i]
  qqplot(rexp(length(pvalues), rate = log(10)),
         -log10(pvalues), xlab = "Expected quantile",
         main = paste0("Calibrated p-values for ", colnames(expl)[i],
                       " (GIF = ", sprintf("%.2f", pv$gif[i]), ")"),
         pch = 19, cex = .4)
  abline(0,1)
}

dev.off()

## output all candidates together (to be deprecated)
# cands <- pv$calibrated.pvalue[sapply(1:nrow(pv$calibrated.pvalue), function(i) any(pv$calibrated.pvalue[i,] < 0.05)),]
# cands <- cands %>% 
#   as.data.frame() %>% 
#   mutate(sequence = rownames(cands),
#          bedname = gsub(".bed", "", basename(bedfile)),
#          segment_nb = gsub(".txt", "", basename(kmerfile))) %>% 
#   relocate(bedname, segment_nb, sequence)
# 
# cat("Number of candidate k-mers:", nrow(cands))
# 
# candfile <-  file.path(outdir, "candidates.csv")
# write.table(cands, candfile, quote = FALSE, sep="\t", row.names = F, col.names = F)

## output candidates separately for each explanatory variable (under construction)

cat("variable\tgif_calibration_kmers\tgif_calibration_outliers\tfdr_control_kmer\tfdr_control_outliers\n")
for (c in 1:ncol(pv$calibrated.pvalue)) {
  # select kmers with pvalue < 0.05
  pv_cands <- rownames(pv$calibrated.pvalue)[pv$calibrated.pvalue[,c] < 0.05]
  pv_cands_df <- data.frame(bedname = gsub(".bed", "", basename(bedfile)),
                         segment_nb = gsub(".csv", "", basename(kmerfile)),
                         sequence = ifelse(length(pv_cands) == 0, NA, pv_cands),
                         outlier = ifelse(length(pv_cands) == 0, NA, pv_cands %in% selective_kmers))
  
  # try with FDR control of 0.05
  qv <- qvalue(pv$calibrated.pvalue[,c], fdr.level = 0.01, pfdr = T)
  fdr_cands <- rownames(pv$calibrated.pvalue)[qv$significant]
  fdr_cands_df <- data.frame(bedname = gsub(".bed", "", basename(bedfile)),
                             segment_nb = gsub("(.txt|.csv)", "", basename(kmerfile)),
                             sequence = ifelse(length(fdr_cands) == 0, NA, fdr_cands),
                             outlier = ifelse(length(fdr_cands) == 0, NA, fdr_cands %in% selective_kmers))
  
  
  # cat("Number of candidate k-mers associated with", colnames(pv$calibrated.pvalue)[c], "is", length(cands), "\n")
  # cat("Number of outlier k-mers associated with", colnames(pv$calibrated.pvalue)[c], "is", length(cands[cands %in% selective_kmers]), "\n")
  cat(colnames(pv$calibrated.pvalue)[c], 
      length(pv_cands), 
      length(pv_cands[pv_cands %in% selective_kmers]),
      length(fdr_cands), 
      length(fdr_cands[fdr_cands %in% selective_kmers]),
      # kmerfile,
      sep = "\t"
      )
  cat("\n")
  
  pvcandfile <-  file.path(outdir, paste0("candidates_gif_", colnames(pv$calibrated.pvalue)[c],".csv"))
  write.table(pv_cands_df, pvcandfile, quote = FALSE, sep = "\t", row.names = F, col.names = F)

  fdrcandfile <- file.path(outdir, paste0("candidates_fdr_", colnames(pv$calibrated.pvalue)[c],".csv"))
  write.table(fdr_cands_df, fdrcandfile, quote = FALSE, sep = "\t", row.names = F, col.names = F)
}

# tot_cands <- pv$calibrated.pvalue[sapply(1:nrow(pv$calibrated.pvalue), function(i) any(pv$calibrated.pvalue[i,] < 0.05)),]
# cat("Number of total candidate k-mers:", nrow(tot_cands), "\n")


# bedfile <- "./gea/tmp/3.TABLE2BED/output_file.0.bed"
# kmerfile <- "./gea/tmp/5.RANGES/output_file.0/1.txt"

# bedfile <- "/shared/projects/most_kmer/kiss_af/output/3.TABLE2BED/output_file.0.bed"
# kmerfile <- "/shared/projects/most_kmer/kiss_af/output/5.RANGES/output_file.0/1.txt"
# climFile <- "/shared/projects/most_kmer/kiss_af/samples/climatic_variables.csv"
# clim <- fread(climFile)
# selectivefile <- "/shared/projects/most_kmer/kiss_af/output/6.PCADAPT/output_file.0_1_BH0.05.pcadapt_outliers.csv"
