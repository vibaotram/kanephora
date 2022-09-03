
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
              help = "number of ancestral groups used as latent factor in lfmm")
)



myArgs <- parse_args(
  OptionParser(usage = "%prog [options]", option_list = option_list,
               description = "run lfmm using genotype in bed file and climatic variables in the clim file"))

bedfile <- myArgs$bedFile
bed <- BEDMatrix(bedfile, simple_names = T)
outdir <- myArgs$outdir
dir.create(outdir, showWarnings = F, recursive = T)
kmerfile <- myArgs$kmerList
kmer <- scan(kmerfile)
clim <- fread(myArgs$climFile)
clim <- clim %>% column_to_rownames(colnames(clim)[1])
pc <- myArgs$pca
vi <- myArgs$variable
K <- myArgs$numgroup

## remove individuals missing in the climfile
## keep kmers in the kmerlist
submat <- bed[which(rownames(bed) %in% rownames(clim)), sort(kmer)]

#######   #######   #######

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
  print("run lfmm with all the provided climatic data")
  expl <- clim
}


# lfmm
sub_lfmm <- lfmm_ridge(Y = submat, X = expl, K = K)
pv <- lfmm_test(Y = submat, X = expl, 
                     lfmm = sub_lfmm, 
                     calibrate = "gif")

lfmmplot <- file.path(outdir, "lfmm_results.pdf")
pdf(file = lfmmplot)
for (i in 1:ncol(expl)) {
  pvalues <- pv$calibrated.pvalue[,i]
  hist(pvalues, 
       main = paste0("Calibrated p-values for ", colnames(expl)[i], 
                     " (GIF = ", sprintf("%.2f", pv$gif[i]), ")"),
       xlab = "p-values")
  # qqplot(rexp(length(pvalues), rate = log(10)),
  #        -log10(pvalues), xlab = "Expected quantile",
  #        pch = 19, cex = .4)
}

# plot(-log10(pv_pca2$calibrated.pvalue[,1]), 
#      pch = 19, 
#      cex = .2, 
#      xlab = "k-mer", ylab = "-Log P",
#      col = "grey")
# points(-log10(pv_pca2$calibrated.pvalue[,2]), 
#        pch = 19, 
#        cex = .2, 
#        col = "grey")
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
for (c in 1:ncol(pv$calibrated.pvalue)) {
  cands <- rownames(pv$calibrated.pvalue)[pv$calibrated.pvalue[,c] < 0.05]
  cands_df <- data.frame(bedname = gsub(".bed", "", basename(bedfile)),
                         segment_nb = gsub(".txt", "", basename(kmerfile)),
                         sequence = cands)
  
  cat("Number of candidate k-mers associated with", colnames(pv$calibrated.pvalue)[c], "is", length(cands), "\n")
  
  candfile <-  file.path(outdir, paste0("candidates_", colnames(pv$calibrated.pvalue)[c],".csv"))
  write.table(cands_df, candfile, quote = FALSE, sep = "\t", row.names = F, col.names = F)
}

tot_cands <- pv$calibrated.pvalue[sapply(1:nrow(pv$calibrated.pvalue), function(i) any(pv$calibrated.pvalue[i,] < 0.05)),]
cat("Number of total candidate k-mers:", nrow(tot_cands), "\n")

# bedfile <- "./gea/tmp/3.TABLE2BED/output_file.0.bed"
# kmerfile <- "./gea/tmp/5.RANGES/output_file.0/1.txt"
