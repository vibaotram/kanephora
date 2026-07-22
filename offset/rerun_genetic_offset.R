library(data.table)
library(tidyverse)
library(sf)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(aplot)

## genetic offset function
go_ex <- function(Y, X, U, K, X.new) {
  mod.lm = lm(Y ~ X + U)
  B = t(mod.lm$coefficients[-c(1, (2+ncol(X)):(1+ncol(X)+K)),])
  B[is.na(B)] <- 0
  genetic.gap = rowSums(((X.new - X) %*% t(B))^2)/nrow(B) 
  rona = rowSums(abs((X.new - X) %*% t(B)))/nrow(B)
  return(list(offset = genetic.gap, distance = rona))
}


## input
#### significant kmers
mat_list <- list.files("/shared/projects/most_kmer/afterkiss/gea/lfmm_rda/out", "candidates_matrix.txt",
                       full.names = T, recursive = T)

can_mat_list <- lapply(mat_list, function(mf) {
  print(which(mat_list == mf))
  can <- read.delim(mf, row.names = 1)
  # print(nrow(can))
  if (nrow(can) == 56) return(can)
})

can_mat_list <- can_mat_list[!unlist(lapply(can_mat_list, is.null))]
can_mat <- do.call(cbind, can_mat_list)

#### current climate
clim_file <- "/shared/projects/most_kmer/afterkiss/bioclim/updated_bioclim_AF.csv"
clim <- read.csv(clim_file)
clim <- clim %>% 
  column_to_rownames(colnames(clim)[1])
clim <- clim[, !colnames(clim) %in% c("bio_10", "bio_1", "bio_3", "bio_14", "bio_9", "bio_16")]
pres_clim <- matrix(unlist(clim[, grep("bio", colnames(clim))]), nrow = nrow(clim))

#### new climate
env_files <- list.files("/shared/projects/most_kmer/afterkiss/bioclim", "_640occ_VN.csv",
                        full.names = T)

#### latent factor
latfac_file <- "/shared/projects/most_kmer/afterkiss/gea/lfmm_rda/latent_factors_15_9_updated.txt"
latfac <- read.table(latfac_file)
latfac <- as.matrix(latfac)
k <- ncol(latfac)

## GO

gap_list <- lapply(env_files, function(e) {
  print(e)
  new_clim <- read.csv(e)
  colnames(new_clim) <- gsub("bioc.+_", "bio_", colnames(new_clim))
  new_clim <- new_clim[, !colnames(new_clim) %in% c("bio_10", "bio_1", "bio_3", "bio_14", "bio_9", "bio_16")]
  new_clim <- new_clim %>% 
    relocate(str_sort(colnames(new_clim), numeric = T)) 
  # new_clim <- matrix(rep(as.numeric(new_clim[l, grep("bio", colnames(new_clim))]), nrow(pres_clim)), nrow = nrow(pres_clim), byrow = T)
  
  gap <- do.call(rbind, lapply(1:nrow(new_clim), function(l) {
    print(l)
    pred_clim <- as.numeric(matrix(rep(new_clim[l, grep("bio", colnames(new_clim))], 
                                       nrow(clim)), nrow = nrow(clim), byrow = T))
    gap <- go_ex(Y = as.matrix(can_mat),
                 X = pres_clim, 
                 U = as.matrix(latfac),
                 K = k,
                 X.new = pred_clim)
    gap_occ <- data.frame(ind = rownames(clim),
                          group = clim$snp_group,
                          occurrence = new_clim[l, "Label"],
                          province = new_clim[l, "province"],
                          district = new_clim[l, "district"],
                          lat = new_clim[l, "Lat"],
                          long = new_clim[l, "Long"],
                          model = gsub("wc2.1_30s_bioc_|_640occ_VN.csv", "", basename(e)),
                          offset = gap$offset,
                          rona = gap$distance,
                          env_dist = sqrt(rowSums((pred_clim - pres_clim)^2)))
    
  }))
  return(gap)
})


all_gap <- do.call(rbind, gap_list)
all_gap %>% group_by(model, location) %>% 
  slice_min(offset) %>% 
  group_by(group) %>% 
  summarise(n = n())

snmf_assign <- read.csv("/shared/projects/most_kmer/afterkiss/structure/af_samples_ancestry_kmers.csv")

all_gap <- left_join(all_gap %>% select(-group), 
          snmf_assign %>% rename(ind = Label),
          by = "ind")

# colnames(all_gap)[colnames(all_gap) == "scenario"] <- "model"

write.csv(all_gap, "offset/all_gap_rerun.csv", quote = F, row.names = F)

## plot
group_col <- c("#30123BFF", "#28BBECFF", "#94e637", "#FB8022FF", "#7A0403FF", "gray", "gray50", "black")
names(group_col) <- c("ER", "OB", "C", "AG", "D", "admixed", "VN present (1970-2000)", "VN future (2041-2060)")

group_col_af <- c("#30123BFF", "#28BBECFF", "#A2FC3CFF", "#FB8022FF", "#7A0403FF", "gray")
names(group_col_af) <- c("ER", "OB", "C", "AG", "D", "hybrid")

best_group <- all_gap %>% 
  filter(group != "hybrid") %>% 
  group_by(model, location) %>% 
  slice_min(offset)

ggplot(best_group) + 
  facet_wrap(vars(model)) + 
  geom_point(aes(x = long, y = lat, col = group)) +
  scale_color_manual(values = group_col_af)

## central higland map
vn <- getData("GADM", country = "VN", level = 1)
vn_pol <- st_as_sf(vn)

district <- getData("GADM", country = "VN", level = 2)
dt_pol <- st_as_sf(district)
dt_pol <- dt_pol[dt_pol$NAME_1 %in% c(rob_ch_clim$province, "Lâm Đồng"),]

best_plots <- lapply(unique(best_group$model), function(s) {
  best_model <- best_group %>% filter(model == s)
  pl <- ggplot() +
    geom_sf(data = dt_pol, fill = "white", color = "#77a37b", size = 0.3) +
    geom_sf(data = dt_pol, fill = NA, color = "#255429", linewidth = .6) + 
    geom_point(data = best_model, aes(x = long, y = lat, color = group), size = 1) +
    # show.legend = (m == "present, 1970-2000")) +
    # geom_line(data = mig_present, 
    #           aes(x = Long, y = Lat, 
    #               color = group,
    #               group = occurrence),
    #               alpha = 0.2) +
    scale_color_manual(values = c(group_col[1:5], "hybrid" = "gray")) +
    theme_void() +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.margin = unit(c(0,2,0,2), "cm")) +
    ggtitle(ifelse(grepl("present", s), "Present", gsub("_", ", ", s))) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
  return(pl)
})

plot_list(best_plots[[7]], best_plots[[1]], best_plots[[2]],
          ggplot() + theme_void(), best_plots[[3]], best_plots[[4]],
          ggplot() + theme_void(), best_plots[[5]], best_plots[[6]],
          labels = c(LETTERS[1:3], "", LETTERS[4:5], "", LETTERS[6:7]))
