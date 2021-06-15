# GEt inclusion probabilities ---

## Create inclusion probabilities ####

library( rgdal)
library( sp)
library( raster)

# clear environment ----
rm(list = ls())

# Set working directory ####
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
s.dir <- paste(w.dir, "spatial_data", sep ='/')
d.dir <- paste(w.dir, "data", sep ='/')



#read in data

SWDat <- readRDS( "SWData_forSA.RDS")
sa_rasters <- readRDS(paste(d.dir, "SARasters.RDS", sep='/'))
zones <- readRDS(paste(d.dir, "SAZones.RDS", sep='/'))

# SW marine park polygon ----
mp <- readOGR(paste(s.dir, "SA-mp.shp", sep='/'))
plot(mp)

####  Straw man for numbers of samples in each region ----

straw.nums <- c(4, 4)  #numbers of drops in and out
straw.props <- straw.nums / sum( straw.nums) # 0.5 / 0.5
names( straw.nums) <- names( straw.props) <- c("smp", "amp")
saveRDS( straw.nums, file="StrawmanNumbers_Zones.RDS")



####  Deciding Bathy cut points for strata : Decided to use slope for this design ----
####  And their numbers of drops


#Bathy.quant <- c(0,0.5,0.9,1)
#Bathy.cuts <- quantile(sa_rasters$bathy, Bathy.quant)#c( -Inf,0.02,0.04,0.08,0.16,Inf)
#Bathy.cuts

Slope.quant <- c(0,0.8,0.95,0.99,1) # this proportions can be changed according to the amount of slope or bathy in the area
Slope.cuts <- quantile(sa_rasters$slope, Slope.quant)
Slope.cuts


#trying to make it so there is no hand-picking (except for the hand-picked function)
tmp <- cumsum(Slope.quant)
#Bathy.targetNums <- rep( floor( 18/8), 4)#floor( ( tmp / sum( tmp))[-1] * 200)#rep( 40, 5)#c( 20,20,30,65,65)
Bathy.targetNums <- rep( floor( 20/5), 2) #  4 4 # 20 BRUVs, in 5 clusters, this is handy when you have several zones w different size
Bathy.targetProps <-  Bathy.targetNums / sum( Bathy.targetNums) # 0.5 0.5


# Proportion of potential sites in each zone ----

# SWDat_small <- SWDat[!is.na( SWDat$Slope),]
# tmp <- colSums( SWDat_small[,c("StateMP", "AMP")], na.rm=TRUE)  # so similar amount of coordinates in each zone     
# tmp[1] # 1710 
# tmp[2] # 1950 
# props <- tmp / nrow( SWDat_small) # state 0.4672131 - common w 0.5327869 
# props <- props / sum( props) # state 0.4672131 - common w 0.5327869


###################################
####  To get cut points
###################################

catS <- cut( sa_rasters$slope, breaks=Slope.cuts, na.rm=TRUE)
plot(catS)

plot( zones$state); plot( catS, add=TRUE); plot( zones$commonw, add=TRUE)
#plot(mp, add=T)

writeRaster(catS, file='Slope_cuts_SA.tif', overwrite=TRUE)

plot(catS)


##################################
####  Within each zone (incl probs)
####  Weight according to straw.props
##################################



inclProbs <- catS
for( zz in c( "state", "commonw")){
  print( zz)
  #if( zz == "MUZ")
  #zoneID <- extract( x=catB, y=zones$MUZ, cellnumbers=TRUE)
  #zoneID <- extract( x=catB, y=zones$MUZ-zones$NPZ, cellnumbers=TRUE)
  #else
  zoneID <- extract( x=catS, y=zones[[zz]], cellnumbers=TRUE)
  propsOfbathy <- table( catS@data@values[zoneID[[1]][,"cell"]]) # it says bathy, but it is slope
  propsOfbathy <- propsOfbathy / sum( propsOfbathy)
  tmp <- Bathy.targetProps / propsOfbathy #the desired inclusion probs (unstandardised)
  for( ii in 1:length( propsOfbathy)){
    inclProbs[zoneID[[1]][,"cell"]][zoneID[[1]][,"value"]==ii] <- tmp[ii]
  }
  inclProbs[zoneID[[1]][,"cell"]][is.na( inclProbs[zoneID[[1]][,"cell"]])] <- 0
  inclProbs[zoneID[[1]][,"cell"]] <- inclProbs[zoneID[[1]][,"cell"]] / sum( inclProbs[zoneID[[1]][,"cell"]])
}
inclProbs@data@values[inclProbs@data@values %in% c(0,1,2,3,4,5,6,7,8)] <- NA  #cheats way to crop
plot( inclProbs)


#standardising so that the zone totals are correct according to straw.props | straw.nums
SMP <- extract( x=catS, y=zones$state, cellnumbers=TRUE)
AMP <- extract( x=catS, y=zones$commonw, cellnumbers=TRUE)


inclProbs@data@values[SMP[[1]][,'cell']] <- inclProbs@data@values[SMP[[1]][,'cell']] * straw.props["smp"]
inclProbs@data@values[AMP[[1]][,'cell']] <- inclProbs@data@values[AMP[[1]][,'cell']] * straw.props["amp"]


plot(inclProbs)


writeRaster(inclProbs, file='inclProbs_SA.tif', overwrite=TRUE)
