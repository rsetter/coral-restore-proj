library(raster)
library(ncdf4)
library(doParallel)
library(velox)
library(sf)
library(foreach)
library(rgis)
library(dplyr)
library(tibble)



###hist

####
#extract reef suit vals
#find storm wind for each point
#classify as suitable or not - save as .csv
###

setwd("E:/storms/historic")


polyfiles <- list.files("E:/reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reef <- st_read("D:/distribution/all_reef_pt.shp")
histfiles <- list.files("E:/storms/historic",pattern="\\.tif$", full.names=T)
nhistfiles <- length(histfiles)
histrout <- paste("hist_",1:nhistfiles*5+1965,"_storm.tif",sep="")
histsuitout <- paste("hist_",1:nhistfiles*5+1965,"_suit.tif",sep="")
histsuitcsvout <- paste("hist_",1:nhistfiles*5+1965,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")



start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:nhistfiles){
  storm <- raster(histfiles[i])
  reefext <- rgis::fast_extract(sf=reef,ras=storm,funct=function(x){max(x,na.rm=T)},parallel=T)
  gc()
  dfext <- as.data.frame(reefext)
  names(dfext) <- c("pointid","grid_code", "geometry","maxwind")
  result <- transform(dfext,pointid=as.numeric(pointid),maxwind=as.numeric(maxwind))
  result$suit <- with(result, ifelse(maxwind >113, 0, 1)) #if windspeed is >113 knots, it's unsuitable. otherwise suitable
  result$suit <- with(result, ifelse(is.na(maxwind),1,suit)) #if no windspeed is recorded for the pixel, then it's suitable
  stormpt <- inner_join(reef,result,by="pointid")
  stormdf <- as.data.frame(stormpt)
  write.csv(stormdf,histsuitcsvout[i], row.names=FALSE)
  rasterize(stormpt,r008333,field="maxwind",fun=max,na.rm=T,filename=histrout[i],overwrite=T)
  rasterize(stormpt,r008333,field="suit",fun=max,na.rm=T,filename=histsuitout[i],overwrite=T)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/storms/historic",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_code", "geometry","maxwind","suit")
  b<-merge(a,filter,by="pointid")
  b$storms_coral <- b$suit * b$iscoral
  c <- count(b, vars="storms_coral")
  y <- add_column(c,year=paste0(i*5+1965),.before=T)
  count <- rbind(y,count)
} 




#####
##test: define threshold for PDI (power dissipation index)
#####
setwd("E:/storms")


reef <- st_read("D:/distribution/all_reef_pt.shp")
storm <- brick("D:/storms/projected/storm45_5yr.nc")
stormr <- rotate(storm)
rcp45rout <- paste("rcp45_",2000+(1:19+1)*5,"_storm.tif",sep="")
rcp45suitout <- paste("rcp45_",2000+(1:19+1)*5,"_suit.tif",sep="")
rcp45suitcsvout <- paste("rcp45_",2000+(1:19+1)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")



start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
reefext <- rgis::fast_extract(sf=reef,ras=stormr,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:2){
  x <- paste0("stormr_X",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
  names(dfi) <- c("pointid","grid_code", "geometry","maxwind")
  dfi$suit <- with(dfi, ifelse(maxwind >274628329, 0, 1)) #if windspeed is >x m3/s2, it's unsuitable. otherwise suitable
                                                          #test1: 93346429; test2: 159220840; test3: 178317062; test4: 274628329 (mean+sd)
  dfi$suit <- with(dfi, ifelse(is.na(maxwind),1,suit)) #if no windspeed is recorded for the pixel, then it's suitable  
  stormpt <- inner_join(reef,dfi,by="pointid")
  stormdf <- as.data.frame(stormpt)
  write.csv(stormdf,rcp45suitcsvout[i], row.names=FALSE)
  rasterize(stormpt,r008333,field="maxwind",fun=max,na.rm=T,filename=rcp45rout[i],overwrite=T)
  rasterize(stormpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp45suitout[i],overwrite=T)
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/storms",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_code", "geometry","maxwind","suit")
  b<-merge(a,filter,by="pointid")
  b$storms_coral <- b$suit * b$iscoral
  c <- count(b, vars=storms_coral)
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  count <- rbind(y,count)
} 






##########################################################


###RCP26

####
#extract reef suit vals
#find storm wind for each point
#classify as suitable or not - save as .csv
###

setwd("E:/storms/rcp26")


reef <- st_read("D:/distribution/all_reef_pt.shp")
storm <- brick("D:/storms/projected/storm26_5yr.nc")
stormr <- rotate(storm)
rcp26rout <- paste("rcp26_",2000+(1:19+1)*5,"_storm.tif",sep="")
rcp26suitout <- paste("rcp26_",2000+(1:19+1)*5,"_suit.tif",sep="")
rcp26suitcsvout <- paste("rcp26_",2000+(1:19+1)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")



start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
reefext <- rgis::fast_extract(sf=reef,ras=stormr,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:19){
  x <- paste0("stormr_X",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
  names(dfi) <- c("pointid","grid_code", "geometry","maxwind")
  dfi$suit <- with(dfi, ifelse(maxwind >274628329, 0, 1)) #if windspeed is >274628329 m3/s2, it's unsuitable. otherwise suitable 
                                                          #old threshold 93,346,429 m3/s2 too low. 
  dfi$suit <- with(dfi, ifelse(is.na(maxwind),1,suit)) #if no windspeed is recorded for the pixel, then it's suitable  
  stormpt <- inner_join(reef,dfi,by="pointid")
  stormdf <- as.data.frame(stormpt)
  write.csv(stormdf,rcp26suitcsvout[i], row.names=FALSE)
  rasterize(stormpt,r008333,field="maxwind",fun=max,na.rm=T,filename=rcp26rout[i],overwrite=T)
  rasterize(stormpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp26suitout[i],overwrite=T)
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/storms/rcp26",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_code", "geometry","maxwind","suit")
  b<-merge(a,filter,by="pointid")
  b$storms_coral <- b$suit * b$iscoral
  c <- count(b, vars=storms_coral)
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  count <- rbind(y,count)
} 








###RCP45

####
#extract reef suit vals
#find storm wind for each point
#classify as suitable or not - save as .csv
###

setwd("E:/storms/rcp45")


reef <- st_read("D:/distribution/all_reef_pt.shp")
storm <- brick("D:/storms/projected/storm45_5yr.nc")
stormr <- rotate(storm)
rcp45rout <- paste("rcp45_",2000+(1:19+1)*5,"_storm.tif",sep="")
rcp45suitout <- paste("rcp45_",2000+(1:19+1)*5,"_suit.tif",sep="")
rcp45suitcsvout <- paste("rcp45_",2000+(1:19+1)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")



start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
reefext <- rgis::fast_extract(sf=reef,ras=stormr,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:19){
  x <- paste0("stormr_X",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
  names(dfi) <- c("pointid","grid_code", "geometry","maxwind")
  dfi$suit <- with(dfi, ifelse(maxwind >274628329, 0, 1)) #if windspeed is >274628329 m3/s2, it's unsuitable. otherwise suitable
  dfi$suit <- with(dfi, ifelse(is.na(maxwind),1,suit)) #if no windspeed is recorded for the pixel, then it's suitable  
  stormpt <- inner_join(reef,dfi,by="pointid")
  stormdf <- as.data.frame(stormpt)
  write.csv(stormdf,rcp45suitcsvout[i], row.names=FALSE)
  rasterize(stormpt,r008333,field="maxwind",fun=max,na.rm=T,filename=rcp45rout[i],overwrite=T)
  rasterize(stormpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp45suitout[i],overwrite=T)
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/storms/rcp45",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_code", "geometry","maxwind","suit")
  b<-merge(a,filter,by="pointid")
  b$storms_coral <- b$suit * b$iscoral
  c <- count(b, vars=storms_coral)
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  count <- rbind(y,count)
} 




###RCP85

####
#extract reef suit vals
#find storm wind for each point
#classify as suitable or not - save as .csv
###

setwd("E:/storms/rcp85")


reef <- st_read("D:/distribution/all_reef_pt.shp")
storm <- brick("D:/storms/projected/storm85_5yr.nc")
stormr <- rotate(storm)
rcp85rout <- paste("rcp85_",2000+(1:19+1)*5,"_storm.tif",sep="")
rcp85suitout <- paste("rcp85_",2000+(1:19+1)*5,"_suit.tif",sep="")
rcp85suitcsvout <- paste("rcp85_",2000+(1:19+1)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")



start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
reefext <- rgis::fast_extract(sf=reef,ras=stormr,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:19){
  x <- paste0("stormr_X",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
  names(dfi) <- c("pointid","grid_code", "geometry","maxwind")
  dfi$suit <- with(dfi, ifelse(maxwind >274628329, 0, 1)) #if windspeed is >274628329 m3/s2, it's unsuitable. otherwise suitable
  dfi$suit <- with(dfi, ifelse(is.na(maxwind),1,suit)) #if no windspeed is recorded for the pixel, then it's suitable  
  stormpt <- inner_join(reef,dfi,by="pointid")
  stormdf <- as.data.frame(stormpt)
  write.csv(stormdf,rcp85suitcsvout[i], row.names=FALSE)
  rasterize(stormpt,r008333,field="maxwind",fun=max,na.rm=T,filename=rcp85rout[i],overwrite=T)
  rasterize(stormpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp85suitout[i],overwrite=T)
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/storms/rcp85",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_code", "geometry","maxwind","suit")
  b<-merge(a,filter,by="pointid")
  b$storms_coral <- b$suit * b$iscoral
  c <- count(b, vars=storms_coral)
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  count <- rbind(y,count)
} 

