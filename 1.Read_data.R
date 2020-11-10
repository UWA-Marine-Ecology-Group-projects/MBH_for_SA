### Read data ####

library( rgdal)
library( sp)
library( raster)


# Set working directory ####
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
s.dir <- paste(w.dir, "spatial_data", sep ='/')


# read poly of marine parks ----
list.files(s.dir)
mp <- readOGR(paste(s.dir, "SA-mp.shp", sep='/'))
mp$ORIG_NAME@data
plot(mp[1,])
plot(mp[2,])
plot(mp[3,])

# read polys of state and commonwealth mps ---
smp <- readOGR(paste(s.dir, "StateMP_SA.shp", sep='/'))
amp <- readOGR(paste(s.dir, "AustralianMP_SA.shp", sep='/'))



#########################
#read in the country boundary -- for plotting mostly
#coastLine <- readOGR( dsn="C:/Users/00093391/Dropbox/UWA/Research Associate/MBHpackage/Ningaloo19_Data/nsaasr9nnd_02211a04es_geo___ (1)")
#proj4string( coastLine) <- CRS("+init=epsg:4283")
#coastLine <- spTransform( coastLine, "+init=epsg:4326")

############################
#read in survey areas data


# zones is a list of polygons

zones <- list()
zones$state <- smp
zones$commonw <- amp
zones$both <- mp


#intial look to see area
plot( zones$both, border='black')
plot( zones$state, add=TRUE, col='orange')
plot( zones$commonw, add=TRUE, col='green')



#### Read bathy ----
dir(s.dir)
b <- raster(paste(s.dir, "bathy-for-mbh-SA.tif", sep='/'))
plot(b)
plot(zones$both, add=T)
# crop to area --
b2 <- crop(b, zones$both)
plot(b2)
plot(zones$both, add=T) # 100-130m

# save new bathy --
writeRaster(b2, paste(s.dir, "bathy-for-mbh-SA2.tif", sep='/'))

###R read TPI ----
t <- raster(paste(s.dir, "tpi-for-mbh-SA.tif", sep='/'))
plot(t)

s <- raster(paste(s.dir, "slope-for-mbh-SA.tif", sep='/'))
plot(s)
s2 <- crop(s, zones$both)
plot(s2)

a <- raster(paste(s.dir, "aspect-for-mbh-SA.tif", sep='/'))
plot(a)


r<- raster(paste(s.dir, "roughness-for-mbh-SA.tif", sep='/'))
plot(r)


# save new sole --
writeRaster(s2, paste(s.dir, "slope-for-mbh-SA2.tif", sep='/'))



###################################################
#### converting polygons to a common raster.

state_raster <- rasterize( x=zones$state, y=b2, field=zones$state@data[,1], bkg.value=-999, fun="first")
plot(state_raster)
amp_raster <- rasterize( zones$commonw, y=b2, field=zones$commonw@data[,1], bkg.value=-999, fun="first")
plot(amp_raster)

###################################
#convert and combine
tmp1 <- as.data.frame( state_raster, xy=TRUE)
tmp2 <- as.data.frame( amp_raster, xy=TRUE)
tmp3 <- as.data.frame( b2, xy=TRUE)
tmp4 <- as.data.frame( s2, xy=TRUE)



SWDat <- cbind( tmp1, tmp2[,3])
SWDat <- cbind( SWDat, tmp3[,3])
SWDat <- cbind( SWDat, tmp4[,3])

head(SWDat)

colnames( SWDat) <- c("Eastern", "Northing", "StateMP", "AMP", "Bathy", "Slope")


## New directory ----
d.dir <- paste(w.dir, "data", sep='/')
setwd(d.dir)

saveRDS( SWDat, file="SWData_forSA.RDS")

sa_rasters <- list()
sa_rasters$bathy <- b2
sa_rasters$slope <- s2

saveRDS( sa_rasters, file="SARasters.RDS")
saveRDS(zones, file="SAZones.RDS")





