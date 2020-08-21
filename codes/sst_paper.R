library(raster)
library(ncdf4)
library(doParallel)
library(velox)
library(sf)
library(foreach)
library(rgis)
library(dplyr)
library(tibble)


###RCP26

####
#extract reef suit vals
#find sst for each point
#classify as suitable or not - save as .csv
###

setwd("E:/sst/rcp26")


reef <- st_read("D:/distribution/all_reef_pt.shp")
sst <- brick("D:/SST/CMIP5/RCP26/sst26_mean_fiveyr.nc")
sstr <- rotate(sst)
rcp26rout <- paste("rcp26_",2000+(1:19+1)*5,"_sst.tif",sep="")
rcp26suitout <- paste("rcp26_",2000+(1:19+1)*5,"_suit.tif",sep="")
rcp26suitcsvout <- paste("rcp26_",2000+(1:19+1)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

  reefext <- rgis::fast_extract(sf=reef,ras=sstr,funct=function(x){max(x,na.rm=T)},parallel=T)
  dfext <- as.data.frame(reefext)
  for (i in 1:19){
    x <- paste0("sstr_X",i)
    dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
    names(dfi) <- c("pointid","grid_code", "geometry","sst")
    dfi$suit <- with(dfi, ifelse(sst >= 297, ifelse(sst <= 305,1,0), 0)) #in kelvin. 24 to 32C. 297 to 305K
    sstpt <- inner_join(reef,dfi,by="pointid")
    sstdf <- as.data.frame(sstpt)
    write.csv(sstdf,rcp26suitcsvout[i], row.names=FALSE)
    rasterize(sstpt,r008333,field="sst",fun=max,na.rm=T,filename=rcp26rout[i],overwrite=T)
    rasterize(sstpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp26suitout[i],overwrite=T)
  }

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/sst/rcp26",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_codex","grid_codey","sst","suit")
  b<-merge(a,filter,by="pointid")
  b$sst_coral <- b$suit * b$iscoral
  c <- count(b, vars="sst_coral")
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  count <- rbind(y,count)
} 




###RCP45

####
#extract reef suit vals
#find sst for each point
#classify as suitable or not - save as .csv
###

setwd("E:/sst/rcp45")


reef <- st_read("D:/distribution/all_reef_pt.shp")
sst <- brick("D:/SST/CMIP5/RCP45/sst45_mean_fiveyr.nc")
sstr <- rotate(sst)
rcp45rout <- paste("rcp45_",2000+(1:19+1)*5,"_sst.tif",sep="")
rcp45suitout <- paste("rcp45_",2000+(1:19+1)*5,"_suit.tif",sep="")
rcp45suitcsvout <- paste("rcp45_",2000+(1:19+1)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reef,ras=sstr,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:19){
  x <- paste0("sstr_X",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
  names(dfi) <- c("pointid","grid_code", "geometry","sst")
  dfi$suit <- with(dfi, ifelse(sst >= 297, ifelse(sst <= 305,1,0), 0)) #in kelvin. 24 to 32C. 297 to 305K
  sstpt <- inner_join(reef,dfi,by="pointid")
  sstdf <- as.data.frame(sstpt)
  write.csv(sstdf,rcp45suitcsvout[i], row.names=FALSE)
  rasterize(sstpt,r008333,field="sst",fun=max,na.rm=T,filename=rcp45rout[i],overwrite=T)
  rasterize(sstpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp45suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/sst/rcp45",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_codex","grid_codey","sst","suit")
  b<-merge(a,filter,by="pointid")
  b$sst_coral <- b$suit * b$iscoral
  c <- count(b, vars="sst_coral")
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  count <- rbind(y,count)
} 






###RCP85

####
#extract reef suit vals
#find sst for each point
#classify as suitable or not - save as .csv
###

setwd("E:/sst/rcp85")


reef <- st_read("D:/distribution/all_reef_pt.shp")
sst <- brick("D:/SST/CMIP5/rcp85/sst85_mean_fiveyr.nc")
sstr <- rotate(sst)
rcp85rout <- paste("rcp85_",2000+(1:19+1)*5,"_sst.tif",sep="")
rcp85suitout <- paste("rcp85_",2000+(1:19+1)*5,"_suit.tif",sep="")
rcp85suitcsvout <- paste("rcp85_",2000+(1:19+1)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reef,ras=sstr,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:19){
  x <- paste0("sstr_X",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
  names(dfi) <- c("pointid","grid_code", "geometry","sst")
  dfi$suit <- with(dfi, ifelse(sst >= 297, ifelse(sst <= 305,1,0), 0)) #in kelvin. 24 to 32C. 297 to 305K
  sstpt <- inner_join(reef,dfi,by="pointid")
  sstdf <- as.data.frame(sstpt)
  write.csv(sstdf,rcp85suitcsvout[i], row.names=FALSE)
  rasterize(sstpt,r008333,field="sst",fun=max,na.rm=T,filename=rcp85rout[i],overwrite=T)
  rasterize(sstpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp85suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/sst/rcp85",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_codex","grid_codey","sst","suit")
  b<-merge(a,filter,by="pointid")
  b$sst_coral <- b$suit * b$iscoral
  c <- count(b, vars="sst_coral")
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  count <- rbind(y,count)
} 




###historic ### 

####
#extract reef suit vals
#find sst for each point
#classify as suitable or not - save as .csv
###

setwd("E:/sst/historic")

sst <- brick("D:/SST/OISST/oisstmasked.nc")
sstr <- rotate(sst)


#make brick with 5 year average. 
#raw file has monthly dates 1981/12 to 2019/6 (450 layers). make averages for 1985, 1990, 1995, 2000, 2005, 2010, 2015

sstd <- dropLayer(sstr, c(1:13,433:450)) #drop layers. we don't need the year 1982 or 2018, 2019
mvlist <- unstack(sstd) # now a list of rasters 
grpsize <- 60 # desired layers per stack - 60 months or five years
# uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
names(mvlist) <- rep(1:(ceiling(length(mvlist) / grpsize )), 
                     each = grpsize, length.out = length(mvlist))
# make list of rasters into a list of stacks. basically each year is now one separate stack
stacks <- lapply(unique(names(mvlist)), function(y) {
  b <- mvlist[names(mvlist) == y]
  stack(b)
})
# find the mean temp for 5 years - w is each stack/5year. 
stacksout <- lapply(stacks, function(w) {
  g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
})
sst_five <- stack(stacksout)
writeRaster(sst_five, filename="oisst_fiveyr.tif", overwrite=TRUE)



reef <- st_read("D:/distribution/all_reef_pt.shp")
sst <- sst_five
histrout <- paste("hist_",1980+(1:7)*5,"_sst.tif",sep="")
histsuitout <- paste("hist_",1980+(1:7)*5,"_suit.tif",sep="")
histsuitcsvout <- paste("hist_",1980+(1:7)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reef,ras=sst,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:7){
  x <- paste0("sst_layer.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
  names(dfi) <- c("pointid","grid_code", "geometry","sst")
  dfi$suit <- with(dfi, ifelse(sst >= 24, ifelse(sst <= 32,1,0), 0)) #in celsius. 24 to 32C. 297 to 305K
  sstpt <- inner_join(reef,dfi,by="pointid")
  sstdf <- as.data.frame(sstpt)
  write.csv(sstdf,histsuitcsvout[i], row.names=FALSE)
  rasterize(sstpt,r008333,field="sst",fun=max,na.rm=T,filename=histrout[i],overwrite=T)
  rasterize(sstpt,r008333,field="suit",fun=max,na.rm=T,filename=histsuitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/sst/historic",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  names(a) <- c("pointid","grid_code","grid_codex","grid_codey","sst","suit","geometryx","geometryy")
  b<-merge(a,filter,by="pointid")
  b$sst_coral <- b$suit * b$iscoral
  c <- count(b, vars="sst_coral")
  y <- add_column(c,year=paste0(i*5+1980),.before=T)
  count <- rbind(y,count)
} 
