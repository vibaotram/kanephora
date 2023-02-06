.libPaths("/shared/home/baotram/R/x86_64-pc-linux-gnu-library/4.1")

setwd("/shared/projects/most_kmer/afterkiss/bioclim")

if (!require("raster")) install.packages('raster')

library(raster)
library(data.table)
library(dplyr)
library(stringr)

samples <- fread("./lfmm_samples.csv")
samples <- samples %>% filter(!is.na(Lat))


## bioclim data
# climfile <- "./wc2.1_30s_bio.zip" # destination path to download the elevation data,  change if needed
# download.file("https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_bio.zip",
#               destfile = climfile,
#               method = "curl")
# unzip(climfile)
climdata_dir <- "./" # path to the worldclim dataset folder
climdata <- list.files(climdata_dir, pattern = "wc2.1_30s_bio_\\d+.tif", full.names = T)
clim <- raster::stack(climdata)

## download elevation data
elevfile <- "./wc2.1_30s_elev.zip" # destination path to download the elevation data,  change if needed
# download.file("https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_elev.zip", 
#               destfile = elevfile,
#               method = "curl")
# unzip(elevfile)
elev <- raster(gsub(".zip", ".tif", elevfile)) 

## extract values for points
envirvar <- do.call(rbind, lapply(1:nrow(samples), function(i) {
  # clim <- getData(name = "worldclim", var = "bio", res = 0.5, lat = na.omit(samples$Lat)[i], lon = na.omit(samples$Long)[i])
  # elev <- getData(name = "SRTM", lat = na.omit(samples$Lat)[i], lon = na.omit(samples$Long)[i])
  
  point <- SpatialPoints(data.frame(long = na.omit(samples$Long)[i], lat = na.omit(samples$Lat)[i]), proj4string = clim@crs)
  clim_values <- raster::extract(clim, point)
  elev_value <- raster::extract(elev, point)
  colnames(clim_values) <- gsub("wc2.1_30s_", "", colnames(clim_values))
  clim_values <- clim_values[, str_sort(colnames(clim_values), numeric = T)]
  df <- data.frame(Label = samples$Label[i], 
                   t(clim_values), 
                   elevation = elev_value,
                   Long = samples$Long[i], Lat = samples$Lat[i],
                   snp_group = samples$snp_group[i])
  return(df)
}))

write.csv(envirvar, "./climatic_variables_af_wc2.1.csv", row.names = F)
write.csv(envirvar[,1:21], "../gea/lfmm/explanatory_wc2.1.csv", row.names = F)