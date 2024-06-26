.libPaths("/shared/home/baotram/R/x86_64-pc-linux-gnu-library/4.1")

setwd("/shared/projects/most_kmer/afterkiss/bioclim")

if (!require("raster")) install.packages('raster')

library(raster)
library(data.table)
library(dplyr)


occ_file <- "occ_coffee_vnm.csv"
occurence <- read.csv(occ_file)

vn <- getData("GADM", country = "VN", level = 2)

## bioclim data
climdata_chelsa <- list.files(".", pattern = "CHELSA_bio\\d+_2041-2070_gfdl-esm4_ssp585_V.2.1.tif", full.names = T)
bioclim_chelsa <- raster::stack(climdata_chelsa)

envirvar_chelsa <- do.call(rbind, lapply(1:nrow(occurence), function(i) {
  print(i)
  point <- SpatialPoints(data.frame(long = na.omit(occurence$Lon)[i], lat = na.omit(occurence$Lat)[i]), proj4string = bioclim_chelsa@crs)
  clim_values <- raster::extract(bioclim_chelsa, point)
  geo <- raster::extract(vn, 
                         SpatialPoints(data.frame(
                           long = na.omit(occurence$Lon)[i], 
                           lat = na.omit(occurence$Lat)[i]
                           ), proj4string = vn@proj4string))
  colnames(clim_values) <- gsub("(CHELSA_|_2041.2070_ipsl.cm6a.lr_ssp585_V.2.1)", "", colnames(clim_values))
  df <- data.frame(Label = paste0("occ_", i), 
                   clim_values, 
                   Long = occurence$Lon[i], Lat = occurence$Lat[i],
                   snp_group = geo$NAME_0,
                   district = geo$VARNAME_2,
                   province = geo$NAME_1)
  return(df)
}))

write.csv(envirvar_chelsa, "./climatic_variables_vn_future_gfdl.csv", row.names = F)
# write.csv(envirvar_chelsa[,1:21], "./lfmm/explanatory_af_chelsa.csv", row.names = F)