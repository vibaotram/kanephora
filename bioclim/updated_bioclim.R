library(raster)
library(xlsx)
library(dplyr)
library(vegan)

## bioclim data
climfile <- "./wc2.1_30s_bio.zip" # destination path to download the elevation data,  change if needed
download.file("https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_bio.zip", 
              destfile = climfile,
              method = "curl")
unzip(climfile)
climdata_dir <- "." # path to the worldclim dataset folder
climdata <- list.files(climdata_dir, pattern = "wc2.1_30s_bio_\\d+.tif", full.names = T)
bioclim <- raster::stack(climdata)

## download elevation data
elevfile <- "./wc2.1_30s_elev.zip" # destination path to download the elevation data,  change if needed
download.file("https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_elev.zip", 
              destfile = elevfile,
              method = "curl")
unzip(elevfile)
elev <- raster(gsub(".zip", ".tif", elevfile)) 

## samples
af_samples <- read.xlsx("bioclim/@passeport-climatic_variables_af_viet.xlsx", sheetIndex = 1)

upd_sample <- af_samples %>% filter(is.na(bio_1)) %>% filter(!is.na(Long))

upd_point <- SpatialPoints(data.frame(long = upd_sample$Long, lat = upd_sample$Lat), proj4string = bioclim@crs)
upd_clim_values <- raster::extract(bioclim, upd_point)

af_samples[is.na(af_samples$bio_1) & !is.na(af_samples$Long),3:21] <- upd_clim_values[,stringr::str_sort(colnames(upd_clim_values), numeric = T)]

upd_elev <- raster::extract(elev, upd_point)
af_samples[is.na(af_samples$elevation) & !is.na(af_samples$Long),22] <-upd_elev

write.csv(af_samples, "bioclim/updated_bioclim_AF.csv", row.names = F, quote = F)

## pca
af_samples <- read.csv("bioclim/updated_bioclim_AF.csv")
af_env <- af_samples %>% select(Label, bio_1:bio_19)
env_pca <- rda(af_env[,-1], scale = T)
summary(env_pca)$cont
round(scores(env_pca, choices=1:10, display="species", scaling=0), digits=3)
screeplot(env_pca)
biplot(env_pca)

## correlation between variables
library(corrplot)
library(caret)
clim_cor <- cor(clim)

pdf("/shared/projects/most_kmer/afterkiss/bioclim/updated_bioclim_correlation.pdf")
corrplot(clim_cor, method = "number", type = "upper", number.cex = 0.7, tl.col = "black")
dev.off()

## remove variables with cor > 0.9
high_cor <- findCorrelation(clim_cor, cutoff = 0.9)
colnames(clim)[high_cor]
# "bio_10", "bio_1", "bio_3", "bio_14", "bio_9", "bio_16"