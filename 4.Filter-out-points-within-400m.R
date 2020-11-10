#### Get the clustered points design and filterout the points that are within 250 m of eachother ####

library( MBHdesign)
library( parallel)
library( class)
library( fields)
library( pdist)
library( raster)
library(rgeos)
library( rgdal)
library( sp)


# clear environment ----
rm(list = ls())


# Set working directory ####
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#.dir <- paste(w.dir, "plots", sep ='/')
s.dir <- paste(w.dir, "spatial_data", sep='/')
#d.dir <- paste(w.dir, "Fish-HWY", sep='/')
d.dir <- paste(w.dir, "data", sep='/')
o.dir <- paste(w.dir, "outputs", sep='/')

## load raster of incl probs ----
ip <- raster(paste(d.dir, "inclProbs_SA.tif", sep='/'))
plot(ip)
llcrs <- proj4string(ip)

b <- raster(paste(s.dir, "bathy-for-mbh-SA2.tif", sep='/'))
plot(b)

s <- raster(paste(s.dir, "slope-for-mbh-SA2.tif", sep='/'))
plot(s)

## read the spatial points created with MBH ----

p1 <- readOGR(paste(o.dir, "SA-bruvs-design-1.shp", sep='/')) # 96 features

proj4string(p1) <- proj4string(ip)

p2 <- readOGR(paste(o.dir, "SA-bruvs-design-2.shp", sep='/')) # 96 features

proj4string(p2) <- proj4string(ip)

p3 <- readOGR(paste(o.dir, "SA-bruvs-design3.shp", sep='/')) # 96 features

proj4string(p3) <- proj4string(ip)

p4 <- readOGR(paste(o.dir, "SA-bruvs-design4.shp", sep='/')) # 96 features

proj4string(p5) <- proj4string(ip)

p5 <- readOGR(paste(o.dir, "SA-bruvs-design5.shp", sep='/')) # 96 features

proj4string(p5) <- proj4string(ip)

p6 <- readOGR(paste(o.dir, "SA-bruvs-design6.shp", sep='/')) # 96 features

proj4string(p6) <- proj4string(ip)

p7 <- readOGR(paste(o.dir, "SA-bruvs-design7.shp", sep='/')) # 96 features

proj4string(p7) <- proj4string(ip)



## read cluster centres ----
c1 <- readOGR(paste(o.dir, "SA-clusters-design-1.shp", sep='/'))

c2 <- readOGR(paste(o.dir, "SA-clusters-design-2.shp", sep='/'))

c3 <- readOGR(paste(o.dir, "SA-clusters-design3.shp", sep='/'))

c4 <- readOGR(paste(o.dir, "SA-clusters-design4.shp", sep='/'))

c5 <- readOGR(paste(o.dir, "SA-clusters-design5.shp", sep='/'))

c6 <- readOGR(paste(o.dir, "SA-clusters-design6.shp", sep='/'))

c7 <- readOGR(paste(o.dir, "SA-clusters-design7.shp", sep='/'))




## read object in utm ----
ga <- raster(paste(s.dir, "ga4858_grid1_MSL.tiff", sep='/'))
gacrs <-  proj4string(ga)

## transform the points into UTM --
p1u <- spTransform(p1, gacrs)
p2u <- spTransform(p2, gacrs)
p3u <- spTransform(p3, gacrs)
p4u <- spTransform(p4, gacrs)
p5u <- spTransform(p5, gacrs)
p6u <- spTransform(p6, gacrs)
p7u <- spTransform(p7, gacrs)

# read zones ---
mp <- readOGR(paste(s.dir, "SA-mp.shp", sep='/'))

## calculate if 2 points fall within 250m of eachother ----
# https://gis.stackexchange.com/questions/102796/remove-points-within-x-distance

## p1 ----
p1_matrix <- gWithinDistance(p1u, dist = 500, byid = TRUE)
diag(p1_matrix) <- NA
p1_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p1_matrix[lower.tri(p1_matrix, diag=TRUE)] <- NA
p1_matrix

colSums(p1_matrix, na.rm=TRUE) == 0
v1 <- colSums(p1_matrix, na.rm=TRUE) == 0
p1u[v1, ] # 42 features left


## Save points of design ---
swbruvs <- p1u[v1, ]
swbruvs

swbruvs <- spTransform(swbruvs, llcrs)
swbruvs

# plot --
plot(ip)
plot(swbruvs, pch=20, col="black", add=T)
plot(c1, pch=20, cex = 2, col='blue', add=T)
plot(mp, add=T)
#plot(c1, pch=20, col="blue", add=T)



# save points in  latlong ----

### CHANGE THIS NO. EVERY TIME ----
nodesign <- "design1"
named <- paste("SA-bruvs", nodesign, sep='-')

writeOGR(swbruvs, dsn=paste(o.dir, 'filtered', sep='/'), layer=named, driver="ESRI Shapefile", overwrite_layer=TRUE)

#writeOGR( swbruvsu, dsn=paste(o.dir, 'filtered', sep='/'), layer=paste("InNOutMP-Bruvs-d3_utm", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
#writeOGR( swbruvs, dsn=paste(o.dir, 'filtered', sep='/'), layer=paste("InNOutMP-Bruvs-d3_latlong", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)

## PLOT---

f.dir <- paste(o.dir,'filtered', sep='/')
pdf(paste(f.dir,  paste(named, "pdf", sep='.'), sep='/'), height=7, width=8)
plot(b, main =named, col = rev(brewer.pal(20, "RdYlBu")))
plot(mp, add=T)
points(swbruvs, pch=21, cex = 0.9, bg='white')
points(c1, pch=21, cex=1.5, bg = 'black')
#plot(swnp, add=T)
#points(wp, add = T, pch=20, cex = 0.7, col='green')
dev.off()



## p2 ----
p2_matrix <- gWithinDistance(p2u, dist = 400, byid = TRUE)
diag(p2_matrix) <- NA
p2_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p2_matrix[lower.tri(p2_matrix, diag=TRUE)] <- NA
p2_matrix

colSums(p2_matrix, na.rm=TRUE) == 0
v2 <- colSums(p2_matrix, na.rm=TRUE) == 0
p2u[v2, ] # 14 features left

## Save points of design ---
swbruvs <- p2u[v2, ]
swbruvs

swbruvs <- spTransform(swbruvs, llcrs)
swbruvs


# save ----

### CHANGE THIS NO. EVERY TIME ----
nodesign <- "design2"
named <- paste("SA-bruvs", nodesign, sep='-')

writeOGR(swbruvs, dsn=paste(o.dir, 'filtered', sep='/'), layer=named, driver="ESRI Shapefile", overwrite_layer=TRUE)

#writeOGR( swbruvsu, dsn=paste(o.dir, 'filtered', sep='/'), layer=paste("InNOutMP-Bruvs-d3_utm", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
#writeOGR( swbruvs, dsn=paste(o.dir, 'filtered', sep='/'), layer=paste("InNOutMP-Bruvs-d3_latlong", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)

## PLOT---

f.dir <- paste(o.dir,'filtered', sep='/')
pdf(paste(f.dir,  paste(named, "pdf", sep='.'), sep='/'), height=7, width=8)
plot(b, main =named, col = rev(brewer.pal(20, "RdYlBu")))
plot(mp, add=T)
points(swbruvs, pch=21, cex = 0.9, bg='white')
points(c2, pch=21, cex=1.5, bg = 'black')
#plot(swnp, add=T)
#points(wp, add = T, pch=20, cex = 0.7, col='green')
dev.off()





## p3 ----
p3_matrix <- gWithinDistance(p3u, dist = 400, byid = TRUE)
diag(p3_matrix) <- NA
p3_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p3_matrix[lower.tri(p3_matrix, diag=TRUE)] <- NA
p3_matrix

colSums(p3_matrix, na.rm=TRUE) == 0
v3 <- colSums(p3_matrix, na.rm=TRUE) == 0
p3[v3, ] # 15 features left

## Save points of design ---
swbruvs <- p3u[v3, ]
swbruvs

swbruvs <- spTransform(swbruvs, llcrs)
swbruvs

# save ----

### CHANGE THIS NO. EVERY TIME ----
nodesign <- "design3"
named <- paste("SA-bruvs", nodesign, sep='-')

writeOGR(swbruvs, dsn=paste(o.dir, 'filtered', sep='/'), layer=named, driver="ESRI Shapefile", overwrite_layer=TRUE)

#writeOGR( swbruvsu, dsn=paste(o.dir, 'filtered', sep='/'), layer=paste("InNOutMP-Bruvs-d3_utm", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
#writeOGR( swbruvs, dsn=paste(o.dir, 'filtered', sep='/'), layer=paste("InNOutMP-Bruvs-d3_latlong", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)

## PLOT---

f.dir <- paste(o.dir,'filtered', sep='/')
pdf(paste(f.dir,  paste(named, "pdf", sep='.'), sep='/'), height=7, width=8)
plot(b, main =named, col = rev(brewer.pal(20, "RdYlBu")))
plot(mp, add=T)
points(swbruvs, pch=21, cex = 0.9, bg='white')
points(c3, pch=21, cex=1.5, bg = 'black')
#plot(swnp, add=T)
#points(wp, add = T, pch=20, cex = 0.7, col='green')
dev.off()



## p4 ----
p4_matrix <- gWithinDistance(p4u, dist = 400, byid = TRUE)
diag(p4_matrix) <- NA
p4_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p4_matrix[lower.tri(p4_matrix, diag=TRUE)] <- NA
p4_matrix

colSums(p4_matrix, na.rm=TRUE) == 0
v4 <- colSums(p4_matrix, na.rm=TRUE) == 0
p4[v4, ] # 12 features left

## Save points of design ---
swbruvs <- p4u[v4, ]
swbruvs

swbruvs <- spTransform(swbruvs, llcrs)
swbruvs

# save ----

### CHANGE THIS NO. EVERY TIME ----
nodesign <- "design4"
named <- paste("SA-bruvs", nodesign, sep='-')

writeOGR(swbruvs, dsn=paste(o.dir, 'filtered', sep='/'), layer=named, driver="ESRI Shapefile", overwrite_layer=TRUE)

#writeOGR( swbruvsu, dsn=paste(o.dir, 'filtered', sep='/'), layer=paste("InNOutMP-Bruvs-d3_utm", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
#writeOGR( swbruvs, dsn=paste(o.dir, 'filtered', sep='/'), layer=paste("InNOutMP-Bruvs-d3_latlong", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)

## PLOT---

f.dir <- paste(o.dir,'filtered', sep='/')
pdf(paste(f.dir,  paste(named, "pdf", sep='.'), sep='/'), height=7, width=8)
plot(b, main =named, col = rev(brewer.pal(20, "RdYlBu")))
plot(mp, add=T)
points(swbruvs, pch=21, cex = 0.9, bg='white')
points(c4, pch=21, cex=1.5, bg = 'black')
#plot(swnp, add=T)
#points(wp, add = T, pch=20, cex = 0.7, col='green')
dev.off()


## p5 ----
p5_matrix <- gWithinDistance(p5u, dist = 400, byid = TRUE)
diag(p5_matrix) <- NA
p5_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p5_matrix[lower.tri(p5_matrix, diag=TRUE)] <- NA
p5_matrix

colSums(p5_matrix, na.rm=TRUE) == 0
v5 <- colSums(p5_matrix, na.rm=TRUE) == 0
p5u[v5, ] # 15 features left

## Save points of design ---
swbruvs <- p5u[v5, ]
swbruvs

swbruvs <- spTransform(swbruvs, llcrs)
swbruvs

# save ----

### CHANGE THIS NO. EVERY TIME ----
nodesign <- "design5"
named <- paste("SA-bruvs", nodesign, sep='-')

writeOGR(swbruvs, dsn=paste(o.dir, 'filtered', sep='/'), layer=named, driver="ESRI Shapefile", overwrite_layer=TRUE)

#writeOGR( swbruvsu, dsn=paste(o.dir, 'filtered', sep='/'), layer=paste("InNOutMP-Bruvs-d3_utm", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
#writeOGR( swbruvs, dsn=paste(o.dir, 'filtered', sep='/'), layer=paste("InNOutMP-Bruvs-d3_latlong", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)

## PLOT---

f.dir <- paste(o.dir,'filtered', sep='/')
pdf(paste(f.dir,  paste(named, "pdf", sep='.'), sep='/'), height=7, width=8)
plot(b, main =named, col = rev(brewer.pal(20, "RdYlBu")))
plot(mp, add=T)
points(swbruvs, pch=21, cex = 0.9, bg='white')
points(c5, pch=21, cex=1.5, bg = 'black')
#plot(swnp, add=T)
#points(wp, add = T, pch=20, cex = 0.7, col='green')
dev.off()


## p6 ----
p6_matrix <- gWithinDistance(p6u, dist = 400, byid = TRUE)
diag(p6_matrix) <- NA
p6_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p6_matrix[lower.tri(p6_matrix, diag=TRUE)] <- NA
p6_matrix

colSums(p6_matrix, na.rm=TRUE) == 0
v6 <- colSums(p6_matrix, na.rm=TRUE) == 0
p6u[v6, ] # 13 features left

## Save points of design ---
swbruvs <- p6u[v6, ]
swbruvs

swbruvs <- spTransform(swbruvs, llcrs)
swbruvs

# save ----

### CHANGE THIS NO. EVERY TIME ----
nodesign <- "design6"
named <- paste("SA-bruvs", nodesign, sep='-')

writeOGR(swbruvs, dsn=paste(o.dir, 'filtered', sep='/'), layer=named, driver="ESRI Shapefile", overwrite_layer=TRUE)

#writeOGR( swbruvsu, dsn=paste(o.dir, 'filtered', sep='/'), layer=paste("InNOutMP-Bruvs-d3_utm", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
#writeOGR( swbruvs, dsn=paste(o.dir, 'filtered', sep='/'), layer=paste("InNOutMP-Bruvs-d3_latlong", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)

## PLOT---

f.dir <- paste(o.dir,'filtered', sep='/')
pdf(paste(f.dir,  paste(named, "pdf", sep='.'), sep='/'), height=7, width=8)
plot(b, main =named, col = rev(brewer.pal(20, "RdYlBu")))
plot(mp, add=T)
points(swbruvs, pch=21, cex = 0.9, bg='white')
points(c6, pch=21, cex=1.5, bg = 'black')
#plot(swnp, add=T)
#points(wp, add = T, pch=20, cex = 0.7, col='green')
dev.off()


## p7 ----
p7_matrix <- gWithinDistance(p7u, dist = 400, byid = TRUE)
diag(p7_matrix) <- NA
p7_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p7_matrix[lower.tri(p7_matrix, diag=TRUE)] <- NA
p7_matrix

colSums(p7_matrix, na.rm=TRUE) == 0
v7 <- colSums(p7_matrix, na.rm=TRUE) == 0
p7u[v7, ] # 12 features left

## Save points of design ---
swbruvs <- p7u[v7, ]
swbruvs

swbruvs <- spTransform(swbruvs, llcrs)
swbruvs

# save ----

### CHANGE THIS NO. EVERY TIME ----
nodesign <- "design7"
named <- paste("SA-bruvs", nodesign, sep='-')

writeOGR(swbruvs, dsn=paste(o.dir, 'filtered', sep='/'), layer=named, driver="ESRI Shapefile", overwrite_layer=TRUE)

#writeOGR( swbruvsu, dsn=paste(o.dir, 'filtered', sep='/'), layer=paste("InNOutMP-Bruvs-d3_utm", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
#writeOGR( swbruvs, dsn=paste(o.dir, 'filtered', sep='/'), layer=paste("InNOutMP-Bruvs-d3_latlong", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)

## PLOT---

f.dir <- paste(o.dir,'filtered', sep='/')
pdf(paste(f.dir,  paste(named, "pdf", sep='.'), sep='/'), height=7, width=8)
plot(b, main =named, col = rev(brewer.pal(20, "RdYlBu")))
plot(mp, add=T)
points(swbruvs, pch=21, cex = 0.9, bg='white')
points(c7, pch=21, cex=1.5, bg = 'black')
#plot(swnp, add=T)
#points(wp, add = T, pch=20, cex = 0.7, col='green')
dev.off()



## Plot slope ----

f.dir <- paste(o.dir,'filtered', sep='/')
pdf(paste(f.dir,  paste("Slope", "pdf", sep='.'), sep='/'), height=7, width=8)
plot(s, main ="Slope", col = rev(brewer.pal(20, "RdYlBu")))
plot(mp, add=T)
#points(swbruvs, pch=21, cex = 0.9, bg='white')
#points(c7, pch=21, cex=1.5, bg = 'black')
#plot(swnp, add=T)
#points(wp, add = T, pch=20, cex = 0.7, col='green')
dev.off()

## Plot bathy ----

f.dir <- paste(o.dir,'filtered', sep='/')
pdf(paste(f.dir,  paste("Bathy", "pdf", sep='.'), sep='/'), height=7, width=8)
plot(b, main ="Bathymetry", col = rev(brewer.pal(20, "RdYlBu")))
plot(mp, add=T)
#points(swbruvs, pch=21, cex = 0.9, bg='white')
#points(c7, pch=21, cex=1.5, bg = 'black')
#plot(swnp, add=T)
#points(wp, add = T, pch=20, cex = 0.7, col='green')
dev.off()




