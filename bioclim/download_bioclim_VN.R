.libPaths("/shared/home/baotram/R/x86_64-pc-linux-gnu-library/4.1")

setwd("/shared/projects/most_kmer/afterkiss/bioclim")

if (!require("raster")) install.packages('raster')

library(raster)
library(data.table)
library(dplyr)
library(stringr)

arg <- commandArgs(trailingOnly = T)

# occ_file <- "occ_coffee_vnm.csv"
# occurence <- read.csv(occ_file)
occurrence <- read.csv("./wc2.1_30s_bioc_present_1970-2000_640occ_VN.csv")

vn <- getData("GADM", country = "VN", level = 2)

## bioclim data
# climdata_dir <- "./" # path to the worldclim dataset folder
climdata <- arg[1] #"wc2.1_30s_bioc_CNRM-CM6-1_ssp126_2041-2060.tif"
clim <- raster::stack(climdata, native = F)

envirvar <- do.call(rbind, lapply(1:nrow(occurrence), function(i) {
  print(i)
  point <- SpatialPoints(data.frame(long = na.omit(occurrence$Long)[i], lat = na.omit(occurrence$Lat)[i]), proj4string = clim@crs)
  clim_values <- raster::extract(clim, point)
  colnames(clim_values) <- gsub(".+_", "bio_", colnames(clim_values))
  clim_values <- clim_values[, str_sort(colnames(clim_values), numeric = T)]
  geo <- raster::extract(vn, 
                         SpatialPoints(data.frame(
                           long = na.omit(occurrence$Long)[i], 
                           lat = na.omit(occurrence$Lat)[i]
                         ), proj4string = vn@proj4string))
  df <- data.frame(Label = occurrence$Label[i], 
                   t(clim_values), 
                   Long = occurrence$Lon[i], Lat = occurrence$Lat[i],
                   snp_group = geo$NAME_0,
                   district = geo$VARNAME_2,
                   province = geo$NAME_1)
  return(df)
}))


write.csv(envirvar, gsub(".tif", "_640occ_VN.csv", climdata), row.names = F, quote = F)
# write.csv(envirvar_chelsa[,1:21], "./lfmm/explanatory_af_chelsa.csv", row.names = F)
