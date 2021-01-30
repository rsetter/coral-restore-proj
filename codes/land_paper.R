library(raster)
library(ncdf4)
library(doParallel)
library(velox)
library(sf)
library(foreach)
library(rgis)
library(dplyr)
library(tibble)

###
#SSP1 - RCP26
###

setwd("E:/land_hr/ssp1-rcp26")

###
#calculate percent of unsuitable land per cell
##
reefbuff <- raster("E:/pop_hr/reeflandbuffer.tif")
ssp1unsuitfiles <- paste("ssp1_",2000+1:10*10,"_us.tif",sep="") 

#test for 2005 (run 2005 separately bc this data annoyingly has an uneven temporal scale from 2005, 2010, 2020, etc)
c <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP1-Baseline-v.1.0.1.nc",level=1,band=1) #cropland
cb <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP1-Baseline-v.1.0.1.nc",level=2,band=1) #cropland-bioenergy
g <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP1-Baseline-v.1.0.1.nc",level=3,band=1) #grassland
b <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP1-Baseline-v.1.0.1.nc",level=7,band=1) #built-up
tu <- overlay(c,cb,g,b,fun=function(w,x,y,z,na.rm=T){return(w+x+y+z)},filename="ssp1_2005_us.tif") #sum of unsuitable land

for (i in 2:11){
  
  c <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP1-Baseline-v.1.0.1.nc",level=1,band=i) #cropland
  cb <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP1-Baseline-v.1.0.1.nc",level=2,band=i) #cropland-bioenergy
  g <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP1-Baseline-v.1.0.1.nc",level=3,band=i) #grassland
  b <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP1-Baseline-v.1.0.1.nc",level=7,band=i) #built-up
  tu <- overlay(c,cb,g,b,fun=function(w,x,y,z,na.rm=T){return(w+x+y+z)},filename=ssp1unsuitfiles[i-1]) #add up all unsuitable land
  
}




#calculate 5 year averages
ssp1files <- list.files("E:/land_hr/ssp1-rcp26", pattern="_us.tif", full.names=T)
ssp15yrnames <- paste("ssp1_20",0:9,"5_us.tif",sep="")

for (i in 2:10){
  x <- raster(ssp1files[i])
  n <- i+1
  y <- raster(ssp1files[n])
  xy5 <- overlay(x,y,fun=function(x){ mean(x,na.rm=T)},filename=ssp15yrnames[i])
}


####
#extract reef suit vals
#find maximum unsuitable land % within 50km radius of reef site - save max % unsuitability as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("E:/reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reef <- st_read("D:/distribution/all_reef_pt.shp")
ssp1files <- list.files("E:/land_hr/ssp1-rcp26",pattern="_us", full.names=T)
nssp1files <- length(ssp1files)
ssp1rout <- paste("ssp1_",1995+(1:nssp1files+1)*5,"_maxunsuit.tif",sep="")
ssp1suitout <- paste("ssp1_",1995+(1:nssp1files+1)*5,"_suit.tif",sep="")
ssp1suitcsvout <- paste("ssp1_",1995+(1:nssp1files+1)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:20){
  land <- raster(ssp1files[i])
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxunsuit")
  result <- transform(dfext,pointid=as.numeric(pointid),maxunsuit=as.numeric(maxunsuit))
  result$suit <- with(result,ifelse(maxunsuit <0.5, 1, 0)) #if the percentage of unsuitable land is less than 0.5, it is suitable (1)
  result$suit <- with(result, ifelse(is.na(maxunsuit),1,suit)) #if land is NA, then land is too far away to return value so is suitable (1)
  landpt <- inner_join(reef,result,by="pointid")
  landdf <- as.data.frame(landpt)
  write.csv(landdf,ssp1suitcsvout[i], row.names=FALSE)
  rasterize(landpt,r008333,field="maxunsuit",fun=max,na.rm=T,filename=ssp1rout[i],overwrite=T)
  rasterize(landpt,r008333,field="suit",fun=max,na.rm=T,filename=ssp1suitout[i],overwrite=T)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/land_hr/ssp1-rcp26",pattern="_suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  names(a) <- c("pointid","grid_code","zone","maxunsuit","suit","geometryx","geometryy")
  b<-merge(a,filter,by="pointid")
  b$land_coral <- b$suit * b$iscoral
  c <- count(b, vars="land_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 
 






###
##SSP2 - RCP45
###

setwd("E:/land_hr/ssp2-rcp45")

###
#calculate percent of unsuitable land per cell
##
reefbuff <- raster("E:/pop_hr/reeflandbuffer.tif")
ssp2unsuitfiles <- paste("ssp2_",2000+1:10*10,"_us.tif",sep="") 

#test for 2005 (run 2005 separately bc this data annoyingly has an uneven temporal scale from 2005, 2010, 2020, etc)
c <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP2-Baseline-v.1.0.1.nc",level=1,band=1) #cropland
cb <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP2-Baseline-v.1.0.1.nc",level=2,band=1) #cropland-bioenergy
g <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP2-Baseline-v.1.0.1.nc",level=3,band=1) #grassland
b <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP2-Baseline-v.1.0.1.nc",level=7,band=1) #built-up
tu <- overlay(c,cb,g,b,fun=function(w,x,y,z,na.rm=T){return(w+x+y+z)},filename="ssp2_2005_us.tif") #sum of unsuitable land

for (i in 2:11){
  
  c <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP2-Baseline-v.1.0.1.nc",level=1,band=i) #cropland
  cb <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP2-Baseline-v.1.0.1.nc",level=2,band=i) #cropland-bioenergy
  g <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP2-Baseline-v.1.0.1.nc",level=3,band=i) #grassland
  b <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP2-Baseline-v.1.0.1.nc",level=7,band=i) #built-up
  tu <- overlay(c,cb,g,b,fun=function(w,x,y,z,na.rm=T){return(w+x+y+z)},filename=ssp2unsuitfiles[i-1]) #add up all unsuitable land
  
}




#calculate 5 year averages
ssp2files <- list.files("E:/land_hr/ssp2-rcp45", pattern="_us.tif", full.names=T)
ssp25yrnames <- paste("ssp2_20",0:9,"5_us.tif",sep="")

for (i in 2:10){
  x <- raster(ssp2files[i])
  n <- i+1
  y <- raster(ssp2files[n])
  xy5 <- overlay(x,y,fun=function(x){ mean(x,na.rm=T)},filename=ssp25yrnames[i])
}


####
#extract reef suit vals
#find maximum unsuitable land % within 50km radius of reef site - save max % unsuitability as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("E:/reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reef <- st_read("D:/distribution/all_reef_pt.shp")
ssp2files <- list.files("E:/land_hr/ssp2-rcp45",pattern="_us", full.names=T)
nssp2files <- length(ssp2files)
ssp2rout <- paste("ssp2_",1995+(1:nssp2files+1)*5,"_maxunsuit.tif",sep="")
ssp2suitout <- paste("ssp2_",1995+(1:nssp2files+1)*5,"_suit.tif",sep="")
ssp2suitcsvout <- paste("ssp2_",1995+(1:nssp2files+1)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:20){
  land <- raster(ssp2files[i])
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxunsuit")
  result <- transform(dfext,pointid=as.numeric(pointid),maxunsuit=as.numeric(maxunsuit))
  result$suit <- with(result,ifelse(maxunsuit <0.5, 1, 0)) #if the percentage of unsuitable land is less than 0.5, it is suitable (1)
  result$suit <- with(result, ifelse(is.na(maxunsuit),1,suit)) #if land is NA, then land is too far away to return value so is suitable (1)
  landpt <- inner_join(reef,result,by="pointid")
  landdf <- as.data.frame(landpt)
  write.csv(landdf,ssp2suitcsvout[i], row.names=FALSE)
  rasterize(landpt,r008333,field="maxunsuit",fun=max,na.rm=T,filename=ssp2rout[i],overwrite=T)
  rasterize(landpt,r008333,field="suit",fun=max,na.rm=T,filename=ssp2suitout[i],overwrite=T)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/land_hr/ssp2-rcp45",pattern="_suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  names(a) <- c("pointid","grid_code","zone","maxunsuit","suit","geometryx","geometryy")
  b<-merge(a,filter,by="pointid")
  b$land_coral <- b$suit * b$iscoral
  c <- count(b, vars="land_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 











###
##SSP5 - RCP85 (baseline)
###

setwd("E:/land_hr/ssp5-rcp85")

###
#calculate percent of unsuitable land per cell
##
reefbuff <- raster("E:/pop_hr/reeflandbuffer.tif")
ssp5unsuitfiles <- paste("ssp5_",2000+1:10*10,"_us.tif",sep="")

#test for 2005
c <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP5-Baseline-v.1.0.1.nc",level=1,band=1) #cropland
cb <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP5-Baseline-v.1.0.1.nc",level=2,band=1) #cropland-bioenergy
g <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP5-Baseline-v.1.0.1.nc",level=3,band=1) #grassland
b <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP5-Baseline-v.1.0.1.nc",level=7,band=1) #built-up
tu <- overlay(c,cb,g,b,fun=function(w,x,y,z,na.rm=T){return(w+x+y+z)},filename="ssp5_2005_us.tif") #sum of unsuitable land

for (i in 2:11){
  
  c <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP5-Baseline-v.1.0.1.nc",level=1,band=i) #cropland
  cb <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP5-Baseline-v.1.0.1.nc",level=2,band=i) #cropland-bioenergy
  g <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP5-Baseline-v.1.0.1.nc",level=3,band=i) #grassland
  b <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP5-Baseline-v.1.0.1.nc",level=7,band=i) #built-up
  tu <- overlay(c,cb,g,b,fun=function(w,x,y,z,na.rm=T){return(w+x+y+z)},filename=ssp5unsuitfiles[i-1]) #add up all unsuitable land
  
}




#calculate 5 year averages
ssp5files <- list.files("E:/land_hr/ssp5-rcp85", pattern="_us.tif", full.names=T)
ssp55yrnames <- paste("ssp5_20",0:9,"5_us.tif",sep="")

for (i in 2:10){
  x <- raster(ssp5files[i])
  n <- i+1
  y <- raster(ssp5files[n])
  xy5 <- overlay(x,y,fun=function(x){ mean(x,na.rm=T)},filename=ssp55yrnames[i])
}


####
#extract reef suit vals
#find maximum unsuitable land % within 50km radius of reef site - save max % unsuitability as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("E:/reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reef <- st_read("D:/distribution/all_reef_pt.shp")
ssp5files <- list.files("E:/land_hr/ssp5-rcp85",pattern="_us", full.names=T)
nssp5files <- length(ssp5files)
ssp5rout <- paste("ssp5_",1995+(1:nssp5files+1)*5,"_maxunsuit.tif",sep="")
ssp5suitout <- paste("ssp5_",1995+(1:nssp5files+1)*5,"_suit.tif",sep="")
ssp5suitcsvout <- paste("ssp5_",1995+(1:nssp5files+1)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:20){
  land <- raster(ssp5files[i])
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxunsuit")
  result <- transform(dfext,pointid=as.numeric(pointid),maxunsuit=as.numeric(maxunsuit))
  result$suit <- with(result,ifelse(maxunsuit <0.5, 1, 0)) #if the percentage of unsuitable land is less than 0.5, it is suitable (1)
  result$suit <- with(result, ifelse(is.na(maxunsuit),1,suit)) #if land is NA, then land is too far away to return value so is suitable (1)
  landpt <- inner_join(reef,result,by="pointid")
  landdf <- as.data.frame(landpt)
  write.csv(landdf,ssp5suitcsvout[i], row.names=FALSE)
  rasterize(landpt,r008333,field="maxunsuit",fun=max,na.rm=T,filename=ssp5rout[i],overwrite=T)
  rasterize(landpt,r008333,field="suit",fun=max,na.rm=T,filename=ssp5suitout[i],overwrite=T)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/land_hr/ssp5-rcp85",pattern="_suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  names(a) <- c("pointid","grid_code","zone","maxunsuit","suit","geometryx","geometryy")
  b<-merge(a,filter,by="pointid")
  b$land_coral <- b$suit * b$iscoral
  c <- count(b, vars="land_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 






###
##SSP3 - RCP85
###

setwd("E:/land_hr/ssp3-rcp85")


###
#calculate percent of unsuitable land per cell
##
reefbuff <- raster("E:/pop_hr/reeflandbuffer.tif")
ssp3unsuitfiles <- paste("ssp3_",2000+1:10*10,"_us.tif",sep="") 

#test for 2005 (run 2005 separately bc this data annoyingly has an uneven temporal scale from 2005, 2010, 2020, etc)
c <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP3-Baseline-v.1.0.1.nc",level=1,band=1) #cropland
cb <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP3-Baseline-v.1.0.1.nc",level=2,band=1) #cropland-bioenergy
g <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP3-Baseline-v.1.0.1.nc",level=3,band=1) #grassland
b <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP3-Baseline-v.1.0.1.nc",level=7,band=1) #built-up
tu <- overlay(c,cb,g,b,fun=function(w,x,y,z,na.rm=T){return(w+x+y+z)},filename="ssp3_2005_us.tif") #sum of unsuitable land

for (i in 2:11){
  
  c <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP3-Baseline-v.1.0.1.nc",level=1,band=i) #cropland
  cb <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP3-Baseline-v.1.0.1.nc",level=2,band=i) #cropland-bioenergy
  g <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP3-Baseline-v.1.0.1.nc",level=3,band=i) #grassland
  b <- raster("D:/landuse/aim-projected/AIM-SSPRCP-LUmap-SSP3-Baseline-v.1.0.1.nc",level=7,band=i) #built-up
  tu <- overlay(c,cb,g,b,fun=function(w,x,y,z,na.rm=T){return(w+x+y+z)},filename=ssp3unsuitfiles[i-1]) #add up all unsuitable land
  
}




#calculate 5 year averages
ssp3files <- list.files("E:/land_hr/ssp3-rcp85", pattern="_us.tif", full.names=T)
ssp35yrnames <- paste("ssp3_20",0:9,"5_us.tif",sep="")

for (i in 2:10){
  x <- raster(ssp3files[i])
  n <- i+1
  y <- raster(ssp3files[n])
  xy5 <- overlay(x,y,fun=function(x){ mean(x,na.rm=T)},filename=ssp35yrnames[i])
}


####
#extract reef suit vals
#find maximum unsuitable land % within 50km radius of reef site - save max % unsuitability as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("E:/reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reef <- st_read("D:/distribution/all_reef_pt.shp")
ssp3files <- list.files("E:/land_hr/ssp3-rcp85",pattern="_us", full.names=T)
nssp3files <- length(ssp3files)
ssp3rout <- paste("ssp3_",1995+(1:nssp3files+1)*5,"_maxunsuit.tif",sep="")
ssp3suitout <- paste("ssp3_",1995+(1:nssp3files+1)*5,"_suit.tif",sep="")
ssp3suitcsvout <- paste("ssp3_",1995+(1:nssp3files+1)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:20){
  land <- raster(ssp3files[i])
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxunsuit")
  result <- transform(dfext,pointid=as.numeric(pointid),maxunsuit=as.numeric(maxunsuit))
  result$suit <- with(result,ifelse(maxunsuit <0.5, 1, 0)) #if the percentage of unsuitable land is less than 0.5, it is suitable (1)
  result$suit <- with(result, ifelse(is.na(maxunsuit),1,suit)) #if land is NA, then land is too far away to return value so is suitable (1)
  landpt <- inner_join(reef,result,by="pointid")
  landdf <- as.data.frame(landpt)
  write.csv(landdf,ssp3suitcsvout[i], row.names=FALSE)
  rasterize(landpt,r008333,field="maxunsuit",fun=max,na.rm=T,filename=ssp3rout[i],overwrite=T)
  rasterize(landpt,r008333,field="suit",fun=max,na.rm=T,filename=ssp3suitout[i],overwrite=T)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/land_hr/ssp3-rcp85",pattern="_suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  names(a) <- c("pointid","grid_code","zone","maxunsuit","suit","geometryx","geometryy")
  b<-merge(a,filter,by="pointid")
  b$land_coral <- b$suit * b$iscoral
  c <- count(b, vars="land_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 







###
##historic
###

setwd("E:/land_hr/historic")

###
#calculate percent of unsuitable land per cell. Note: the historic data uses scale 1-100 while projected data used 0-1
#unsuitable land = c3 cropland, c4 cropland, c3 pastureland, c4 pastureland, urban land
#varnames = C3crop, C4crop, C3past, C4past, Urban
##
reefbuff <- raster("E:/pop_hr/reeflandbuffer.tif")
histfiles <- list.files("D:/landuse/historic",pattern="\\.nc$",full.names=T)
histunsuitfiles <- paste("hist_",1975+1:7*5,"_us.tif",sep="")


for (i in 1:7){
  c3 <- raster(histfiles[i], varname="C3crop") #cropland c3
  c4 <- raster(histfiles[i], varname="C4crop") #cropland c4
  p3 <- raster(histfiles[i], varname="C3past") #pastureland c3
  p4 <- raster(histfiles[i], varname="C4past") #pastureland c4
  u <- raster(histfiles[i], varname="Urban") #urban
  tu <- overlay(c3,c4,p3,p4,u,fun=function(v,w,x,y,z,na.rm=T){return(v+w+x+y+z)}) #sum of unsuitable land
  tur <- rotate(tu) #rotate to extent (-180,180,-90,90)
  writeRaster(tur, filename=histunsuitfiles[i])
}



####
#extract reef suit vals
#find maximum unsuitable land % within 50km radius of reef site - save max % unsuitability as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("E:/reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reef <- st_read("D:/distribution/all_reef_pt.shp")
histfiles <- list.files("E:/land_hr/historic",pattern="_us", full.names=T)
nhistfiles <- length(histfiles)
histrout <- paste("hist_",1975+(1:nhistfiles)*5,"_maxunsuit.tif",sep="")
histsuitout <- paste("hist_",1975+(1:nhistfiles)*5,"_suit.tif",sep="")
histsuitcsvout <- paste("hist_",1975+(1:nhistfiles)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:7){
  land <- raster(histfiles[i])
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=land,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxunsuit")
  result <- transform(dfext,pointid=as.numeric(pointid),maxunsuit=as.numeric(maxunsuit))
  result$suit <- with(result,ifelse(maxunsuit <50, 1, 0)) #if the percentage of unsuitable land is less than 50%, it is suitable (1)
  result$suit <- with(result, ifelse(is.na(maxunsuit),1,suit)) #if land is NA, then land is too far away to return value so is suitable (1)
  landpt <- inner_join(reef,result,by="pointid")
  landdf <- as.data.frame(landpt)
  write.csv(landdf,histsuitcsvout[i], row.names=FALSE)
  rasterize(landpt,r008333,field="maxunsuit",fun=max,na.rm=T,filename=histrout[i],overwrite=T)
  rasterize(landpt,r008333,field="suit",fun=max,na.rm=T,filename=histsuitout[i],overwrite=T)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/land_hr/historic",pattern="_suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  names(a) <- c("pointid","grid_code","zone","maxunsuit","suit","geometryx","geometryy")
  b<-merge(a,filter,by="pointid")
  b$land_coral <- b$suit * b$iscoral
  c <- count(b, vars="land_coral")
  y <- add_column(c,year=paste0(i*5+1975),.before=T)
  count <- rbind(y,count)
} 

