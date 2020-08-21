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
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("E:/dhw/rcp26")


reef <- st_read("D:/distribution/all_reef_pt.shp")
dhw <- brick("dhw26_mean_fiveyr.nc")
rcp26rout <- paste("rcp26_",2000+(1:19+1)*5,"_dhw.tif",sep="")
rcp26suitout <- paste("rcp26_",2000+(1:19+1)*5,"_suit.tif",sep="")
rcp26suitcsvout <- paste("rcp26_",2000+(1:19+1)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reef,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:19){
  x <- paste0("dhw_X",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
  names(dfi) <- c("pointid","grid_code", "geometry","dhw")
  dfi$suit <- with(dfi, ifelse(dhw <= 8,1, 0)) #suitable number of dhw is <8
  dhwpt <- inner_join(reef,dfi,by="pointid")
  dhwdf <- as.data.frame(dhwpt)
  write.csv(dhwdf,rcp26suitcsvout[i], row.names=FALSE)
  rasterize(dhwpt,r008333,field="dhw",fun=max,na.rm=T,filename=rcp26rout[i],overwrite=T)
  rasterize(dhwpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp26suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/dhw/rcp26",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_codex","grid_codey","dhw","suit")
  b<-merge(a,filter,by="pointid")
  b$dhw_coral <- b$suit * b$iscoral
  c <- count(b, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  count <- rbind(y,count)
} 




###RCP45

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("E:/dhw/rcp45")


reef <- st_read("D:/distribution/all_reef_pt.shp")
dhw <- brick("dhw45_mean_fiveyr.nc")
rcp45rout <- paste("rcp45_",2000+(1:19+1)*5,"_dhw.tif",sep="")
rcp45suitout <- paste("rcp45_",2000+(1:19+1)*5,"_suit.tif",sep="")
rcp45suitcsvout <- paste("rcp45_",2000+(1:19+1)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reef,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:19){
  x <- paste0("dhw_X",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
  names(dfi) <- c("pointid","grid_code", "geometry","dhw")
  dfi$suit <- with(dfi, ifelse(dhw <= 8,1, 0)) #suitable number of dhw is <8
  dhwpt <- inner_join(reef,dfi,by="pointid")
  dhwdf <- as.data.frame(dhwpt)
  write.csv(dhwdf,rcp45suitcsvout[i], row.names=FALSE)
  rasterize(dhwpt,r008333,field="dhw",fun=max,na.rm=T,filename=rcp45rout[i],overwrite=T)
  rasterize(dhwpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp45suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/dhw/rcp45",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_codex","grid_codey","dhw","suit")
  b<-merge(a,filter,by="pointid")
  b$dhw_coral <- b$suit * b$iscoral
  c <- count(b, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  count <- rbind(y,count)
} 



###RCP85

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("E:/dhw/rcp85")


reef <- st_read("D:/distribution/all_reef_pt.shp")
dhw <- brick("dhw85_mean_fiveyr.nc")
rcp85rout <- paste("rcp85_",2000+(1:19+1)*5,"_dhw.tif",sep="")
rcp85suitout <- paste("rcp85_",2000+(1:19+1)*5,"_suit.tif",sep="")
rcp85suitcsvout <- paste("rcp85_",2000+(1:19+1)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reef,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:19){
  x <- paste0("dhw_X",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
  names(dfi) <- c("pointid","grid_code", "geometry","dhw")
  dfi$suit <- with(dfi, ifelse(dhw <= 8,1, 0)) #suitable number of dhw is <8
  dhwpt <- inner_join(reef,dfi,by="pointid")
  dhwdf <- as.data.frame(dhwpt)
  write.csv(dhwdf,rcp85suitcsvout[i], row.names=FALSE)
  rasterize(dhwpt,r008333,field="dhw",fun=max,na.rm=T,filename=rcp85rout[i],overwrite=T)
  rasterize(dhwpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp85suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/dhw/rcp85",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_codex","grid_codey","dhw","suit")
  b<-merge(a,filter,by="pointid")
  b$dhw_coral <- b$suit * b$iscoral
  c <- count(b, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  count <- rbind(y,count)
} 



###historic

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("E:/dhw/historic")


reef <- st_read("D:/distribution/all_reef_pt.shp")
dhw <- brick("hist_mean_fiveyr.nc")
histrout <- paste("hist_",1980+(1:7)*5,"_dhw.tif",sep="")
histsuitout <- paste("hist_",1980+(1:7)*5,"_suit.tif",sep="")
histsuitcsvout <- paste("hist_",1980+(1:7)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reef,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:7){
  x <- paste0("dhw_X",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
  names(dfi) <- c("pointid","grid_code", "geometry","dhw")
  dfi$suit <- with(dfi, ifelse(dhw <= 8,1, 0)) #suitable number of dhw is <8
  dhwpt <- inner_join(reef,dfi,by="pointid")
  dhwdf <- as.data.frame(dhwpt)
  write.csv(dhwdf,histsuitcsvout[i], row.names=FALSE)
  rasterize(dhwpt,r008333,field="dhw",fun=max,na.rm=T,filename=histrout[i],overwrite=T)
  rasterize(dhwpt,r008333,field="suit",fun=max,na.rm=T,filename=histsuitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/dhw/historic",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_codex","grid_codey","dhw","suit")
  b<-merge(a,filter,by="pointid")
  b$dhw_coral <- b$suit * b$iscoral
  c <- count(b, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+1980),.before=T)
  count <- rbind(y,count)
} 

