library(data.table)
library(raster)

rcp45 <- fread("Genomic_Vulnerability_rcp4.5.txt",
               data.table = F)

source("krige_function_4share.R")
test_raster <- getData('worldclim', var='tmin', res=10)
test_raster <- setValues(test_raster, 1)
plot.krige(long = rcp45$Longitude, lat = rcp45$Latitude, z = rcp45$Genomic.Vulnerability,
           raster = test_raster)