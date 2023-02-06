.libPaths("/shared/home/baotram/R/x86_64-pc-linux-gnu-library/4.1")

setwd("/shared/projects/most_kmer/afterkiss/bioclim")

if (!require("raster")) install.packages('raster')

library(raster)
library(data.table)
library(dplyr)


samples <- fread("./lfmm_samples.csv")
samples <- samples %>% filter(!is.na(Lat))


## bioclim data
climdata_chelsa <- list.files(".", pattern = "CHELSA_bio\\d+_1981-2010_V.2.1.tif", full.names = T)
bioclim_chelsa <- raster::stack(climdata_chelsa)
elevfile <- "./wc2.1_30s_elev.tif"
elev <- raster(elevfile)

envirvar_chelsa <- do.call(rbind, lapply(1:nrow(samples), function(i) {
  # clim <- getData(name = "worldclim", var = "bio", res = 0.5, lat = na.omit(samples$Lat)[i], lon = na.omit(samples$Long)[i])
  # elev <- getData(name = "SRTM", lat = na.omit(samples$Lat)[i], lon = na.omit(samples$Long)[i])
  print(i)
  point <- SpatialPoints(data.frame(long = na.omit(samples$Long)[i], lat = na.omit(samples$Lat)[i]), proj4string = bioclim_chelsa@crs)
  clim_values <- raster::extract(bioclim_chelsa, point)
  elev_value <- raster::extract(elev, point)
  colnames(clim_values) <- gsub("(CHELSA_|_1981.2010_V.2.1)", "", colnames(clim_values))
  df <- data.frame(Label = samples$Label[i], 
                   clim_values, 
                   elevation = elev_value,
                   Long = samples$Long[i], Lat = samples$Lat[i],
                   snp_group = samples$snp_group[i])
  return(df)
}))

write.csv(envirvar_chelsa, "./climatic_variables_af_chelsa.csv", row.names = F)
write.csv(envirvar_chelsa[,1:21], "./lfmm/explanatory_af_chelsa.csv", row.names = F)