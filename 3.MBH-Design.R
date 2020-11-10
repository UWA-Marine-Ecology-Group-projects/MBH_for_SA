### Create inclusion probabilities ####

## for multibeam and lidar #######
#install.packages("MBHdesign")
library( MBHdesign)
library( parallel)
library( class)
library( fields)
#install.packages("pdist")
library( pdist)
library( raster)
library( rgdal)
library( sp)

# clear environment ----
rm(list = ls())

# Set working directory ####
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
s.dir <- paste(w.dir, "spatial_data", sep ='/')
d.dir <- paste(w.dir, "data", sep ='/')


#read in data -----

# read in the inclusion probs

inclProbs <- raster(paste(d.dir, "inclProbs_SA.tif", sep='/'))
plot(inclProbs)
inclProbs <- setValues( inclProbs, values( inclProbs) / sum( values( inclProbs), na.rm=TRUE))
plot(inclProbs)
rootInclProbs <- inclProbs
rootInclProbs <- setValues( rootInclProbs, sqrt( values( rootInclProbs)))
zones <- readRDS( "SAZones.RDS") # this one in different folder
rast <- readRDS("SARasters.RDS")
straw.nums <- readRDS( "StrawmanNumbers_Zones.RDS")
names(straw.nums) <- c("state", "commonw")




############################
####  Spatial sample of new sites
####  from altered incl. probs.
############################

### Here use quasiSamp to get random points ####
## these points will be the center of buffer for transects ###

####  Set the seed for reproducability
#set.seed( 777)
#### HAVE NOT BEEN ABLE TO MAKE THIS FUNCTION WORK ----
newSites <- list(SMP=NULL,AMP=NULL)
for( zz in c("state","commonw")){
  print( zz)
  numby <- floor( (straw.nums[zz]))
  #the number of samples to take (specified minus the legacy number)
  #numby <- floor( (straw.nums[zz] - numRef[zz])/2)
  #numby <- floor( (straw.nums[zz] - numRef[zz]))
  #set up spatial domain
  myZone <- zones[[zz]]
  #if( zz == "AMP"){
  # myZone = zones$AMP - zones$IUCN2
  #set.seed( 747)
  #}
  #tmpIP <- mask( rootInclProbs_agg_100m, myZone)
  tmpIP <- mask( inclProbs, myZone)
  tmpIP <- crop( tmpIP, myZone)
  #take the sample of clusters based on root incl probs
  newSites[[zz]] <- quasiSamp( n=numby, potential.sites=coordinates( tmpIP), inclusion.probs=values(tmpIP), nSampsToConsider=5000)
  
  #plotting (maybe remove at a later date?)
  tmpIPFull <- mask( inclProbs, myZone)
  tmpIPFull <- crop( tmpIPFull, myZone)
  plot( tmpIPFull)
  #plot( legacySites, add=TRUE, pch=1, col='red')
  points( newSites[[zz]][,c("x","y")], pch=20, col='black')
}
newSites <- do.call( "rbind", newSites)
newSites <- SpatialPointsDataFrame( coords=newSites[,c("x","y")], data=newSites, proj4string=CRS(proj4string(inclProbs)))
#some of the spatial balance is not great...  Presumably because the balance of the reference sites is also not great...


# check clusters
plot(inclProbs)
points(newSites, pch=20, cex =2, col="blue")
plot(zones$both, add=T)



###############################
####  Choose new points within clusters
####  Here I need to choose transects not points
##############################

getlocal <- function(ii){
  point <- newSites[ii,c("x","y")]
  r2 <- rasterize( point, inclProbs, field=1)
  pbuf <- buffer( r2, width=1500) ## units are in metres
  buf <- mask( inclProbs, pbuf)
  buffer <- trim(buf, pad=0)
  return( buffer)
}


sampWithOver <- 12 # how many points to use : use many, since then they need to be filtered with a 400 m radius


fullSample <- list()
fullZones <- list()



## Get bruvs accoring to clusters ---
for( ii in 1:nrow(newSites)){
  tmp <- getlocal(ii)
  fullZones[[ii]] <- rownames( newSites@data)[ii]
  tmpm <- raster::as.matrix(tmp)
  tmpm <- t(tmpm)
  tmpdf <- as.data.frame (
    cbind (coordinates (tmp), as.numeric (tmpm)))
  colnames(tmpdf) <- c("x", "y", "inclProbs_design1")
  tmpdf <- tmpdf[ order(tmpdf$y, tmpdf$x),]  # order ascending first by northing and then by easting
  fullSample[[ii]] <- quasiSamp( n=sampWithOver, potential.sites=coordinates(tmp), inclusion.probs=values(tmp), nSampsToConsider=10000)
  #fullSample[[ii]] <- transectSamp( n=sampWithOver, potential.sites = tmpdf[,c(1,2)], 
                                    #potential.sites= tmpdf[,c("x","y")],
                                    #inclusion.probs= incprobdf[,3],
                                    #inclusion.probs= tmpdf[,3],
                                    #control=gb.control
                                    #constrainedSet=gb.constraints.bool
  
  plot( tmp)
  points( fullSample[[ii]]$points[,c("x","y")], pch=20, col='red')
  #plot( legacySites, add=TRUE, pch=4, col='blue')
}


fullSample <- do.call( "rbind", fullSample[1:length(fullSample)])
bruvs <- SpatialPointsDataFrame( coords=fullSample[,c("x","y")], data=fullSample, proj4string=CRS(proj4string(inclProbs)))



plot(inclProbs)
points(bruvs, pch =20)
points(newSites, pch=20, cex =2, col="blue")
plot(zones$both, add=T)

head(newSites) # cluster centres
head(bruvs) # bruvs


# To save ----


### plot nicely w bathy #####
library(dichromat)
library(RColorBrewer)

o.dir <- paste(w.dir, "outputs", sep='/')

### CHANGE THIS NO. EVERY TIME ----
nodesign <- "design7"

namebathy <- paste("SA-bathymetry", nodesign, sep='-')
nameslope <- paste("SA-slope", nodesign, sep='-')


#pal <- colorRampPalette(c("red","blue"))
pdf(paste(o.dir,  paste(namebathy, "pdf", sep='.'), sep='/'), height=7, width=8)
plot(rast$bathy, main =namebathy, col = rev(brewer.pal(20, "RdYlBu")))
plot(zones$both, add=T)
points(bruvs, pch=21, cex = 0.9, bg='white')
points(newSites, pch=21, cex=1.5, bg = 'black')
#plot(swnp, add=T)
#points(wp, add = T, pch=20, cex = 0.7, col='green')
dev.off()

### plot nicely w slope #####

#pal <- colorRampPalette(c("green","orange"))
pdf(paste(o.dir,  paste(nameslope, "pdf", sep='.'), sep='/'), height=7, width=8)
plot(rast$slope, main =nameslope, col = rev(brewer.pal(20, "RdYlBu")))
plot(zones$both, add=T)
points(bruvs, pch=21, cex = 0.9, bg='white')
points(newSites, pch=21, cex=1.5, bg = 'black')
#plot(swnp, add=T)
#points(wp, add = T, pch=20, cex = 0.7, col='green')
dev.off()



####  Write the shape files

writeOGR(bruvs, dsn=o.dir, layer=paste("SA-bruvs", nodesign, sep='-'), driver="ESRI Shapefile", overwrite_layer=TRUE)
#writeOGR(fullSample3, dsn=d.dir, layer=paste( "Bruvs4.8", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
writeOGR(newSites, dsn=o.dir, layer=paste("SA-clusters", nodesign, sep='-'), driver="ESRI Shapefile", overwrite_layer=TRUE)
