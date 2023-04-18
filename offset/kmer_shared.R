## get number kmers shared between groups (all combination of groups)
## input:
#### 1- kmer p/a table
#### 2- sample info
#### 3- number of total kmers (103800513)
#### 4- number of kmers each batch (100000)
#### 5- batch index ($SLURM_ARRAY_TASK_ID)
#### 6- outdir
#### 7- number of cores($SLURM_CPUS_PER_TASK)

library(optparse)
library(data.table)
library(tidyverse)
# library(eulerr)
# library(UpSetR)

option_list <- list(
  make_option(c("-t", "--table"),
              type = "character",
              default = NULL,
              help = "Path to kmer p/a table"),
  make_option(c("-i", "--info"),
              type = "character",
              default = NULL,
              help = "Path to sample info file"),
  make_option(c("-n", "--nk"),
              type = "integer",
              default = NULL,
              help = "number of total kmers"),
  make_option(c("-s", "--sk"),
              type = "integer",
              default = NULL,
              help = "number of kmers each batch"),
  make_option(c("-b", "--batch"),
              type = "integer",
              default = NULL,
              help = "batch index"),
  make_option(c("-c", "--cores"),
              type = "integer",
              default = NULL,
              help = "number of cores"),
  make_option(c("-o", "--outdir"),
              type = "character",
              default = NULL,
              help = "outdir")
)



myArgs <- parse_args(
  OptionParser(usage = "%prog [options]", option_list = option_list,
               description = "get number kmers shared between groups (all combination of groups)"))

print(myArgs)
table <- myArgs$table
samples <- fread(myArgs$info)
nk <- as.numeric(myArgs$nk)
sk <- as.numeric(myArgs$sk)
batch <- as.numeric(myArgs$batch)
cores <- as.numeric(myArgs$cores)
outdir <- myArgs$outdir
# dir.create(outdir, showWarnings = F)

all_comb <- function(n) {
  cb <- lapply(1:length(n), function(i) combn(n, i))
  names(cb) <- 1:length(n)
  return(cb)
}

kmer_inds <- fread(table, drop = 1, data.table = F, nThread = cores, nrows = 0)
ls <- seq((batch - 1)*sk + 1, min(batch*sk, nk))


cat("process kmers from line", ls[1], "to", max(ls), "\n")

kmer_count_gr <- sapply(unique(samples$phenotype_value), function(gr) {
  cat(gr, "- ")
  kmer_gr_table <- fread(table, nThread = cores, 
                         # select = samples$accession_id[samples$phenotype_value == gr],
                         drop = 1, skip = ls[1], nrows = length(ls),
                         data.table = F)
  colnames(kmer_gr_table) <- colnames(kmer_inds)
  kmer_lfmm_ind <- rowSums(kmer_gr_table[,colnames(kmer_gr_table) %in% samples$accession_id[1:53]])
  kmer_gr_table <- kmer_gr_table[which(kmer_lfmm_ind > 0),colnames(kmer_gr_table) %in% samples$accession_id[samples$phenotype_value == gr]]
  cat(ncol(kmer_gr_table), "individuals\n")
  kmer_gr_count <- rowSums(kmer_gr_table)
  # kmer_gr_count[kmer_gr_count > 0] <- 1
  return(rownames(kmer_gr_table)[kmer_gr_count > 0])
}, simplify = F)

# nkmer_grp <- sapply(kmer_count_gr, length)
# saveRDS(nkmer_grp, file.path(outdir, paste0("nb_kmer_group_", batch, ".RDS")))

combo <- all_comb(names(kmer_count_gr))


kmer_shared <- sapply(names(combo), function(c) {
  sapply(1:ncol(combo[[c]]), function(i) {
    g <- combo[[c]][,i]
    # if (length(g) == 1) {
    #   gk <- kmer_count_gr[[g]]
    #   ogk <- unlist(kmer_count_gr[names(kmer_count_gr) != g])
    #   kc <- length(gk[!gk %in% ogk])
    # } else {
    #   gk <- kmer_count_gr[g]
    #   kc <- length(Reduce(intersect, gk))
    # }
    # gk <- unique(unlist(kmer_count_gr[g]))
    gk <- Reduce(intersect, kmer_count_gr[g])
    ogk <- unique(unlist(kmer_count_gr[!names(kmer_count_gr) %in% c(g, "Vietnam")]))
    kc <- length(gk[!gk %in% ogk])
    names(kc) <- paste(g, collapse = "&")
    return(kc)
  })
}, simplify = T, USE.NAMES = F)

saveRDS(kmer_shared, file.path(outdir, paste0("kmer_shared_", batch, ".RDS")))

# kmer_comb <- unlist(kmer_comb)

# kmer_fit <- euler(unlist(kmer_comb), input = c("disjoint", "union"),
#                   shape = "ellipse", loss = "abs")
# 
# plot(kmer_fit, quantities = TRUE,
#      fill = group_col[names(kmer_count_gr)])
# 
# 
# upset(fromExpression(unlist(kmer_comb)), 
#       nsets = 6, nintersects = NA, intersections = unlist(combo[c(1,2,6)]),
#       keep.order = T, 
#       sets.bar.color = group_col[-c(6,8)])


