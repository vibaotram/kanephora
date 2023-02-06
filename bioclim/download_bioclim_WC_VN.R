.libPaths("/shared/home/baotram/R/x86_64-pc-linux-gnu-library/4.1")

setwd("/shared/projects/most_kmer/afterkiss/bioclim")

if (!require("raster")) install.packages('raster')

library(raster)
library(data.table)
library(dplyr)
library(stringr)


occ_file <- "occ_coffee_vnm.csv"
occurence <- read.csv(occ_file)

vn <- getData("GADM", country = "VN", level = 2)

## bioclim data
climdata_dir <- "./" # path to the worldclim dataset folder
climdata <- list.files(climdata_dir, pattern = "wc2.1_30s_bio_\\d+.tif", full.names = T)
clim <- raster::stack(climdata)

elevfile <- "./wc2.1_30s_elev.tif"
elev <- raster(elevfile)

envirvar <- do.call(rbind, lapply(1:nrow(occurence), function(i) {
  print(i)
  point <- SpatialPoints(data.frame(long = na.omit(occurence$Lon)[i], lat = na.omit(occurence$Lat)[i]), proj4string = clim@crs)
  clim_values <- raster::extract(clim, point)
  elev_value <- raster::extract(elev, point)
  colnames(clim_values) <- gsub("wc2.1_30s_", "", colnames(clim_values))
  clim_values <- clim_values[, str_sort(colnames(clim_values), numeric = T)]
  geo <- raster::extract(vn, 
                         SpatialPoints(data.frame(
                           long = na.omit(occurence$Lon)[i], 
                           lat = na.omit(occurence$Lat)[i]
                         ), proj4string = vn@proj4string))
  df <- data.frame(Label = paste0("occ_", i), 
                   t(clim_values), 
                   elevation = elev_value,
                   Long = occurence$Lon[i], Lat = occurence$Lat[i],
                   snp_group = geo$NAME_0,
                   district = geo$VARNAME_2,
                   province = geo$NAME_1)
  return(df)
}))

write.csv(envirvar, "./climatic_variables_vn_wc2.1.csv", row.names = F)
# write.csv(envirvar_chelsa[,1:21], "./lfmm/explanatory_af_chelsa.csv", row.names = F)