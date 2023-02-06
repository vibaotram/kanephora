.libPaths("/shared/home/baotram/R/x86_64-pc-linux-gnu-library/4.1")

setwd("/shared/projects/most_kmer/afterkiss")

if (!require("raster")) install.packages('raster')
if (!require("tibble")) install.packages('tibble')

library(raster)
library(data.table)
library(dplyr)
library(tibble)


## coordinates in vn
vn <- getData("GADM", country = "VN", level = 2)

ch_provinces <- c("Lâm Đồng", "Đắk Lắk", "Đắk Nông", "Kon Tum", "Gia Lai")
# ch_provinces <- "Đắk Lắk"
ch_map <- vn[vn$NAME_1 %in% ch_provinces, ]
length(ch_map$NAME_2)

chdist_coord <- do.call(rbind, lapply(ch_map$NAME_2, function(d) {
  dist_map <- ch_map[ch_map$NAME_2 == d,]
  coord <- data.frame(province = dist_map$NAME_1,
                      district = dist_map$NAME_2,
                      coordinates(dist_map)
  )
  colnames(coord)[3:4] <- c("long", "lat")
  return(coord)
}))


## get bioclim
bioclim <- readRDS("./gea/bioclim_2.1.RDS")
elev_data <- readRDS("./gea/elevation0.5_v2.1.RDS")
chvar <- do.call(rbind, lapply(1:nrow(chdist_coord), function(i) {
  # clim <- getData(name = "worldclim", var = "bio", res = 0.5, lat = na.omit(chdist_coord$lat)[i], lon = na.omit(chdist_coord$long)[i])
  # elev <- getData(name = "SRTM", lat = na.omit(chdist_coord$lat)[i], lon = na.omit(chdist_coord$long)[i])
  point <- SpatialPoints(data.frame(long = na.omit(chdist_coord$long)[i], lat = na.omit(chdist_coord$lat)[i]), proj4string = bioclim@crs)
  clim_values <- raster::extract(bioclim, point)
  elev_value <- raster::extract(elev_data, point)
  colnames(clim_values) <- gsub("_.+", "", colnames(clim_values))
  df <- data.frame(district = chdist_coord$district[i], 
                   clim_values, 
                   elevation = elev_value,
                   Long = chdist_coord$long[i], Lat = chdist_coord$lat[i],
                   group = "VN")
  return(df)
}))
# chvar <- chvar %>% 
#   column_to_rownames("district")

write.table(chvar, "./offset/climatic_variables_2.1_VN.csv", sep = ",",
            row.names = F, quote = F)