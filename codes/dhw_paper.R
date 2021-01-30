library(raster)
library(ncdf4)
library(doParallel)
library(velox)
library(sf)
library(foreach)
library(rgis)
library(dplyr)
library(tibble)
library(plyr)


###RCP26

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("E:/dhw/rcp26")


reefs <- st_read("E:/currentcoral/reefs.shp")
dhw <- brick("dhw26_mean_fiveyr.tif")
#rcp26rout <- paste("rcp26_",2000+(1:19+1)*5,"_dhw.tif",sep="")
#rcp26suitout <- paste("rcp26_",2000+(1:19+1)*5,"_suit.tif",sep="")
rcp26suitcsvout <- paste("rcp26_",2000+(1:20)*5,"_suit.csv",sep="")
#r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:20){
  x <- paste0("dhw_dhw26_mean_fiveyr.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","dhw")
  dfi$suit <- with(dfi, ifelse(dhw <= 8,1, 0)) #suitable number of dhw is <8
  #dhwpt <- inner_join(reef,dfi,by="pointid")
  #dhwdf <- as.data.frame(dhwpt)
  write.csv(dfi,rcp26suitcsvout[i], row.names=FALSE)
  #rasterize(dhwpt,r008333,field="dhw",fun=max,na.rm=T,filename=rcp26rout[i],overwrite=T)
  #rasterize(dhwpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp26suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/dhw/rcp26",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$dhw_coral <- a$suit * a$iscoral
  c <- count(a, vars=dhw_coral)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 




###RCP45

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("E:/dhw/rcp45")


reefs <- st_read("E:/currentcoral/reefs.shp")
dhw <- brick("dhw45_mean_fiveyr.tif")
#rcp45rout <- paste("rcp45_",2000+(1:19+1)*5,"_dhw.tif",sep="")
#rcp45suitout <- paste("rcp45_",2000+(1:19+1)*5,"_suit.tif",sep="")
rcp45suitcsvout <- paste("rcp45_",2000+(1:20)*5,"_suit.csv",sep="")
#r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:20){
  x <- paste0("dhw_dhw45_mean_fiveyr.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","dhw")
  dfi$suit <- with(dfi, ifelse(dhw <= 8,1, 0)) #suitable number of dhw is <8
  #dhwpt <- inner_join(reef,dfi,by="pointid")
  #dhwdf <- as.data.frame(dhwpt)
  write.csv(dfi,rcp45suitcsvout[i], row.names=FALSE)
  #rasterize(dhwpt,r008333,field="dhw",fun=max,na.rm=T,filename=rcp45rout[i],overwrite=T)
  #rasterize(dhwpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp45suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/dhw/rcp45",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$dhw_coral <- a$suit * a$iscoral
  c <- count(a, vars=dhw_coral)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 



###RCP85

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("E:/dhw/rcp85")


reefs <- st_read("E:/currentcoral/reefs.shp")
dhw <- brick("dhw85_mean_fiveyr.nc")
#rcp85rout <- paste("rcp85_",2000+(1:19+1)*5,"_dhw.tif",sep="")
#rcp85suitout <- paste("rcp85_",2000+(1:19+1)*5,"_suit.tif",sep="")
rcp85suitcsvout <- paste("rcp85_",2000+(1:20)*5,"_suit.csv",sep="")
#r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:20){
  x <- paste0("dhw_X",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","dhw")
  dfi$suit <- with(dfi, ifelse(dhw <= 8,1, 0)) #suitable number of dhw is <8
  #dhwpt <- inner_join(reef,dfi,by="pointid")
  #dhwdf <- as.data.frame(dhwpt)
  write.csv(dfi,rcp85suitcsvout[i], row.names=FALSE)
  #rasterize(dhwpt,r008333,field="dhw",fun=max,na.rm=T,filename=rcp85rout[i],overwrite=T)
  #rasterize(dhwpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp85suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/dhw/rcp85",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$dhw_coral <- a$suit * a$iscoral
  c <- count(a, vars=dhw_coral)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 



###historic oisst

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("E:/dhw/historic/empirical")


reefs <- st_read("E:/currentcoral/reefs.shp")
dhw <- brick("hist_mean_fiveyr.tif")
#histrout <- paste("hist_",1980+(1:7)*5,"_dhw.tif",sep="")
#histsuitout <- paste("hist_",1980+(1:7)*5,"_suit.tif",sep="")
histsuitcsvout <- paste("hist_",1980+(1:7)*5,"_suit.csv",sep="")
#r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 3:7){
  x <- paste0("dhw_hist_mean_fiveyr.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","dhw")
  dfi$suit <- with(dfi, ifelse(dhw <= 8,1, 0)) #suitable number of dhw is <8
  #dhwpt <- inner_join(reef,dfi,by="pointid")
  #dhwdf <- as.data.frame(dhwpt)
  write.csv(dfi,histsuitcsvout[i], row.names=FALSE)
  #rasterize(dhwpt,r008333,field="dhw",fun=max,na.rm=T,filename=histrout[i],overwrite=T)
  #rasterize(dhwpt,r008333,field="suit",fun=max,na.rm=T,filename=histsuitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/dhw/historic/empirical",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$dhw_coral <- a$suit * a$iscoral
  c <- count(a, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+1980),.before=T)
  count <- rbind(y,count)
} 





#historic model

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("E:/dhw/historic")


reefs <- st_read("E:/currentcoral/reefs.shp")
dhw <- brick("E:/dhw/historic/dhwhist_mean_tenyr.tif")
#histrout <- paste("hist_",1980+(1:7)*5,"_dhw.tif",sep="")
#histsuitout <- paste("hist_",1980+(1:7)*5,"_suit.tif",sep="")
histsuitcsvout <- paste("hist_",1845+(1:16)*10,"_suit.csv",sep="")
#r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:16){
  x <- paste0("dhw_dhwhist_mean_tenyr.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","dhw")
  dfi$suit <- with(dfi, ifelse(dhw <= 8,1, 0)) #suitable number of dhw is <8
  #dhwpt <- inner_join(reef,dfi,by="pointid")
  #dhwdf <- as.data.frame(dhwpt)
  write.csv(dfi,histsuitcsvout[i], row.names=FALSE)
  #rasterize(dhwpt,r008333,field="dhw",fun=max,na.rm=T,filename=histrout[i],overwrite=T)
  #rasterize(dhwpt,r008333,field="suit",fun=max,na.rm=T,filename=histsuitout[i],overwrite=T)
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
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$dhw_coral <- a$suit * a$iscoral
  c <- count(a, vars=dhw_coral)
  y <- add_column(c,year=paste0(i*10+1845),.before=T)
  count <- rbind(y,count)
} 