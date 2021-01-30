library(raster)
library(ncdf4)
library(doParallel)
library(velox)

setwd("E:/pop_hr/ssp1")

#make area raster
r1.2010 <- raster("D:/human_population/population_projection/ssp1/SSP1_1km/ssp1_total_2010.nc4")
area <- area(r1.2010, filename="area_1km.tif", na.rm=FALSE)

#crop to only the coastal areas near reef sites. then calculate pop density for all years
path.ssp1 <- "D:/human_population/population_projection/ssp1/SSP1_1km/"
ssp1files <- list.files(path.ssp1, pattern="_total_", full.names=T)
nfiles <- length(ssp1files)
ssp1_list  <- paste("ssp1.20",1:nfiles,"0",sep="")
path_out <- "D:/omega_arag/historic/co2sys_input/"
ssp1_dens_list <- paste("ssp1.20",1:nfiles,"0d",sep="")
reefbuff <- raster("E:/pop_hr/reeflandbuffer.tif")
ssp1densfiles <- paste("ssp1_20",1:nfiles,"0_dens.tif",sep="")

for (i in 1:nfiles){
  
  x <- raster(ssp1files[i]) #open each raster
  assign(ssp1_list[i],x)
  
  y <- extend(x,extent(reefbuff),value=NA) #change extent so matches reef buffer 
  extent(y) <- extent(reefbuff)
  z <- overlay(y,reefbuff,fun=function(x,y){(x*y)}) #get rid of inland areas we aren't interested in
  
  d <- overlay(z,area,fun=function(x,y,na.rm=T){return(x/y)},filename=ssp1densfiles[i]) #calculate density
  assign(ssp1_dens_list[i],d)
}

#calculate 5 year averages
ssp1dfiles <- list.files("E:/pop_hr/ssp1", pattern="_dens.tif", full.names=T)
ssp15yrnames <- paste("ssp1_20",1:nfiles,"5_dens.tif",sep="")

for (i in 1:nfiles){
  x <- raster(ssp1dfiles[i])
  n <- i+1
  y <- raster(ssp1dfiles[n])
  xy5 <- overlay(x,y,fun=function(x){ mean(x,na.rm=T)},filename=ssp15yrnames[i])
}


#make brick
# maybe this isn't necessary? It takes forever.
ssp1allfiles <- list.files("E:/pop_hr/ssp1", pattern="_dens.tif",full.names=T)
nallfiles <- length(ssp1allfiles)
ssp1alllist <- paste("ssp1_",1:nallfiles,"d",sep="")
for (i in 1:nallfiles){
  x <- raster(ssp1allfiles[i]) #open each raster
  assign(ssp1alllist[i],x)
}
ssp1fiveyr <- brick(ssp1_1d,ssp1_2d,ssp1_3d,ssp1_4d,ssp1_5d,ssp1_6d,ssp1_7d,ssp1_8d,ssp1_9d,ssp1_10d,ssp1_11d,ssp1_12d,ssp1_13d,ssp1_14d,ssp1_15d,ssp1_16d,
                   ssp1_17d,ssp1_18d,ssp1_19d)
writeRaster(ssp1fiveyr, filename="ssp1_fiveyr_dens.nc", format="CDF", overwrite=TRUE)
