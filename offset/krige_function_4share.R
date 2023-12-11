


plot.krige <- function(long, lat, z, dot.col.factor = NULL, extremal.value=0, raster, kriging.mode = "spherical", pixels=100, range.values=NULL, cols = viridis(20), size_points = 2, map) {
  
  require(raster); require(kriging); require(modeest); require(ggplot2); require(maptools); require(viridis)
  
  coords <- apply(cbind(long, lat), c(1,2), function(x) as.numeric(as.character(x)))
  colnames(coords) <- c("x", "y")
  coords.ini <- as.data.frame(coords)
  
  # add extremal coords
  r <- raster
  ext <- extent(r)
  coords <- rbind(coords, c(ext@xmin, ext@ymax), c(ext@xmax, ext@ymin))
  z <- c(z, extremal.value, extremal.value)
  
  # OBS
  kr <- kriging(coords[,1], coords[,2], z, 
                    model=kriging.mode, 
                    lags=2, 
                    pixels=pixels, 
                    polygons=NULL)
  krmap <- kr$map
  krpoints <- SpatialPoints(krmap[,1:2])
  krmap <- krmap[!is.na(raster::extract(r, krpoints)),]
  krmap <- na.omit(krmap)
  
  names(krmap) <- c("x", "y", "RONA")
  
  if (is.null(range.values)) {
    range.values <- range(krmap$RONA)
  }

  g <- ggplot() +
    geom_tile(data=krmap, aes(x=x, y=y, fill=RONA)) +
    scale_fill_gradientn(
      colours=cols, 
      limits=range.values) +
    theme_void() +
    geom_point(data=coords.ini, aes(x,y), size=size_points, shape=21, fill=dot.col.factor, color = "white", alpha=1) +
    geom_sf(data = af,
            fill = NA,
            # aes(fill = yield_change/present_yield*100),
            color = "gray", size = 0.3) +
    scale_y_continuous(limits = c(-15, 14), expand = c(0, 0)) +
    scale_x_continuous(limits = c(-17, 36), expand = c(0, 0)) +
    # coord_equal() + 
    xlab("Longitude (N)") + ylab("Latitude (E)")

if (F) {
	if (is.null(dot.col.factor)) {
		g <- g + geom_point(data=coords.ini, aes(x, y, ), size=3, shape=19, col="black")
	} else {
		group <- dot.col.factor
		g <- g + geom_point(data=coords.ini, aes(x, y, col = group), size=3, shape=19)
	}
}
  
  print(g)
  return(krmap)
}

if (F) {
  # the background is set as a raster
  # could improve function to deal with polygons
  
  H <- read.table("H.txt")
  load("Tb"); H <- Tb
  
  r <- getData('worldclim', var='alt', res=10)
  plot(r)
  
  tmean <- getData('worldclim', var='tmean', res=10)
  tmax <- getData('worldclim', var='tmax', res=10)
  prec <- getData('worldclim', var='prec', res=10)
  
  tme <- crop(tmean[[1]], extent(-20, 40, -7, 10)); plot(tme); points(H$long, H$lat, pch=19)
  tma <- crop(tmax[[1]], extent(-20, 40, -7, 10)); plot(tma); points(H$long, H$lat, pch=19)
  pr <- crop(prec[[1]], extent(-20, 40, -7, 10)); plot(pr); points(H$long, H$lat, pch=19)
  
  
  coord <- H
  coordinates(coord) <- ~long+lat
  
  coord$tmean <- raster::extract(tme, coord)
  coord$tmax <- raster::extract(tma, coord)
  coord$prec <- raster::extract(pr, coord)
  
  boxplot(coord$tmean)
  boxplot(coord$tmax)
  boxplot(coord$prec*12)
  
  bios <- getData('worldclim', var='bio', res=2.5)
  bi <- crop(bios, extent(-20, 40, -7, 10))
  plot(bi); points(H$long, H$lat, pch=19)
  
  summary(raster::extract(bi[[12]], coord))
  
  for ( i in 1:19) {
    cat(i,".")
    png(paste("bio",i,".png",sep=""),1500,1000)
      #plot(bi[[i]], main=i); 
    plot(crop(bi[[i]], extent(28, 35, -3, 3)), main=i); 
    points(H$long, H$lat, pch=19)
    
    dev.off()
  
  }
  
  af <- crop(r, extent(-20, 60, -40, 50))
  af <- crop(r, extent(-20, 40, -7, 10))
  af[!is.na(af)] <- 1
  plot(af)
  points(H$long, H$lat, pch=19)
  
  # z is the vector of values for each point
  # extremal.value is the value to set at the most extreme point
  
  df <- H
  mh <- aggregate(df$Ho, by=list(long=df$long, lat=df$lat), data=df, FUN=mean)
  
  png(paste("krige.png",sep=""),1500,500)
  plot.krige(long = mh$long, 
             lat = mh$lat, 
             z = mh$x, 
             extremal.value = 0,
             raster = af,
             kriging.mode = "spherical", # spherical | exponential | gaussian
             pixels = 500, # resolution
             range.values = NULL, # either NULL (do not constraint range) or a range to bound mapped z values 
             cols = terrain.colors(20) )
  dev.off()

}
