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
  make_option(c("-s", "--selectiveKmers"),
              type = "character",
              default = NULL,
              help = "file containing selective/outlier kmers"),
  make_option(c("-g", "--backgroundKmers"),
              type = "character",
              action = NULL,
              help = "file containing background kmers in lfmm"),  
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

r <- as.numeric(myArgs$retain)
f <- as.numeric(myArgs$fdr)


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

submat <- bed[which(rownames(bed) %in% rownames(clim)), sort(kmer)]
kmer_count <- colSums(submat)
submat <- submat[, which(kmer_count < nrow(submat)*2 & kmer_count > 0)]

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

expl <- clim[,-20]

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

signif_lfmm <- data.frame()

cat("variable\tsignificant_kmers\tsignificant_outliers\n")
for (c in 1:ncol(pv$calibrated.pvalue)) {
  # select kmers with pvalue < 0.05
  # pv_cands <- rownames(pv$calibrated.pvalue)[pv$calibrated.pvalue[,c] < 0.05]
  # if (length(pv_cands) == 0) {
  #   pv_cands_df <- NULL
  #   } else {
  #     pv_cands_df <- data.frame(bedname = gsub(".bed", "", basename(bedfile)),
  #                               segment_nb = gsub(".csv", "", basename(kmerfile)),
  #                               sequence = pv_cands,
  #                               outlier = pv_cands %in% selective_kmers)
  #   }
  
  
  # filter out NaN pvalues
  pv_filt <- na.omit(as.data.frame(pv$calibrated.pvalue[,c]), cols = 1)
  colnames(pv_filt) <- "pvalue"
  
  # fdr control
  qv <- qvalue(pv_filt$pvalue, fdr.level = 0.01, pfdr = T)
  fdr_cands <- rownames(pv_filt)[qv$significant]
  
  # qv <- qvalue(pv$calibrated.pvalue[,c], fdr.level = 0.01, pfdr = T)
  # fdr_cands <- rownames(pv$calibrated.pvalue)[qv$significant]
  if (length(fdr_cands) == 0) {
    fdr_cands_var <- data.frame(bedname = gsub(".bed", "", basename(bedfile)),
                                segment_nb = gsub("(.txt|.csv)", "", basename(kmerfile)),
                                sequence = NA,
                                outlier = NA,
                                variable = FALSE)
  } else {
    fdr_cands_var <- data.frame(bedname = gsub(".bed", "", basename(bedfile)),
                                segment_nb = gsub("(.txt|.csv)", "", basename(kmerfile)),
                                sequence = fdr_cands,
                                outlier = fdr_cands %in% selective_kmers,
                                variable = TRUE)
  }
  colnames(fdr_cands_var)[5] <- colnames(pv$calibrated.pvalue)[c]
  signif_lfmm <- merge(signif_lfmm, fdr_cands_var, all = T)
  
  cat(
    colnames(pv$calibrated.pvalue)[c],
    length(fdr_cands),
    length(fdr_cands[fdr_cands %in% selective_kmers]),
    sep = "\t"
  )
  cat("\n")
  
  pv_filt$significant <- qv$significant
  pvfile <-  file.path(outdir, paste0("pvalue_", colnames(pv$calibrated.pvalue)[c],".csv"))
  write.table(pv_filt, pvfile, quote = FALSE, sep = "\t", row.names = T, col.names = F)
}

signif_lfmm <- signif_lfmm %>% 
  dplyr::filter(!is.na(outlier))

signif_lfmm[is.na(signif_lfmm)] <- FALSE

fdrcandfile <- file.path(outdir, "candidates_kmers_lfmm.csv")
write.table(signif_lfmm, fdrcandfile, quote = FALSE, sep = "\t", row.names = F, col.names = T)


## rda
sub_rda <- rda(submat ~ ., data = expl, scale = F) #  + Condition(PC1 + PC2)
pdf(file.path(outdir, "rda_results.pdf"))
screeplot(sub_rda)
plot(sub_rda, scaling = 1)

## select candidate kmers
rdadapt <- function(rda,K) {
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
hist(signif_rdadapt$p.values)
stopifnot(nrow(signif_rdadapt) == ncol(submat))
signif_rda <- rownames(sub_rda$CCA$v)[signif_rdadapt$q.values < f]

ggplot() +
  geom_point(aes(x=sub_rda$CCA$v[which(signif_rdadapt[,2] >= 0.01),1]*100, 
                 y=sub_rda$CCA$v[which(signif_rdadapt[,2] >= 0.01),2]*100), col = "gray50", size = .1) +
  geom_point(aes(x=sub_rda$CCA$v[which(signif_rdadapt[,2] < 0.01),1]*100,
                 y=sub_rda$CCA$v[which(signif_rdadapt[,2] < 0.01),2]*100), col = "red", size = .1) +
  geom_segment(aes(xend=sub_rda$CCA$biplot[,1], yend=sub_rda$CCA$biplot[,2], x=0, y=0),
               colour="black", size=0.5, linetype=1,
               arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(aes(x=1.1*sub_rda$CCA$biplot[,1], y=1.1*sub_rda$CCA$biplot[,2],
                label = rownames(sub_rda$CCA$biplot))) +
  xlab("RDA 1") + ylab("RDA 2") +
  theme_minimal() +
  theme(legend.position="none")

dev.off()

cat("Number of candidate k-mers detected by methods:\n")
cat("RDA", length(signif_rda), "\n")

# signif_lfmm <- fread(myArgs$lfmm)
cat("LFMM", nrow(signif_lfmm), "\n")
signif_com <- signif_rda[signif_rda %in% signif_lfmm$sequence]
cat("RDA+LFMM", length(signif_com), sprintf("%.2f", length(signif_com)/nrow(signif_lfmm)),"\n")

signif_lfmm$RDA <- signif_lfmm$sequence %in% signif_rda
write.table(signif_lfmm, file.path(outdir, "lfmm_rda.tsv"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
pval_rdadapt <- data.frame(sequence = rownames(sub_rda$CCA$v),
                           pvalue = signif_rdadapt$p.values,
                           significant = signif_rdadapt$q.values < f,
                           lfmm = rownames(sub_rda$CCA$v) %in% signif_lfmm$sequence)
write.table(pval_rdadapt, file.path(outdir, "candidates_rda.tsv"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
writeLines(signif_com, file.path(outdir, "candidates_rda_lfmm.txt"))

# tot_cands <- pv$calibrated.pvalue[sapply(1:nrow(pv$calibrated.pvalue), function(i) any(pv$calibrated.pvalue[i,] < 0.05)),]
# cat("Number of total candidate k-mers:", nrow(tot_cands), "\n")


# bedfile <- "./gea/tmp/3.TABLE2BED/output_file.0.bed"
# kmerfile <- "./gea/tmp/5.RANGES/output_file.0/1.txt"

# bedfile <- "/shared/projects/most_kmer/kiss_af/output/3.TABLE2BED/output_file.0.bed"
# kmerfile <- "/shared/projects/most_kmer/kiss_af/output/5.RANGES/output_file.0/1.txt"
# climFile <- "/shared/projects/most_kmer/kiss_af/samples/climatic_variables.csv"
# clim <- fread(climFile)
# selectivefile <- "/shared/projects/most_kmer/kiss_af/output/6.PCADAPT/output_file.0_1_BH0.05.pcadapt_outliers.csv"
