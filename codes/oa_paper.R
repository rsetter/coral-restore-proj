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
#find omega arag for each point
#classify as suitable or not - save as .csv
###

setwd("E:/oa/hist")

reef <- st_read("D:/distribution/all_reef_pt.shp")
oa <- brick("D:/omega_arag/historic/co2sys_input/OA_hist.nc")
histrout <- paste("hist_",1975+(1:6)*5,"_oa.tif",sep="")
histsuitout <- paste("hist_",1975+(1:6)*5,"_suit.tif",sep="")
histsuitcsvout <- paste("hist_",1975+(1:6)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reef,ras=oa,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:6){
  x <- paste0("oa_X",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
  names(dfi) <- c("pointid","grid_code", "geometry","oa")
  dfi$suit <- with(dfi, ifelse(oa >= 3.3, 1, 0)) #if omega arag is > 3.3, it is suitable(1). otherwise not (0)
  oapt <- inner_join(reef,dfi,by="pointid")
  oadf <- as.data.frame(oapt)
  write.csv(oadf,histsuitcsvout[i], row.names=FALSE)
  rasterize(oapt,r008333,field="oa",fun=max,na.rm=T,filename=histrout[i],overwrite=T)
  rasterize(oapt,r008333,field="suit",fun=max,na.rm=T,filename=histsuitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/oa/hist",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_codex","grid_codey","oa","suit")
  b<-merge(a,filter,by="pointid")
  b$oa_coral <- b$suit * b$iscoral
  c <- count(b, vars="oa_coral")
  y <- add_column(c,year=paste0(i*5+1975),.before=T)
  count <- rbind(y,count)
} 




###rcp26

####
#extract reef suit vals
#find omega arag for each point
#classify as suitable or not - save as .csv
###

setwd("E:/oa/rcp26")

reef <- st_read("D:/distribution/all_reef_pt.shp")
oa <- brick("D:/omega_arag/co2sys_input/rcp26/OA_rcp26.nc")
rcp26rout <- paste("rcp26_",2005+(1:19)*5,"_oa.tif",sep="")
rcp26suitout <- paste("rcp26_",2005+(1:19)*5,"_suit.tif",sep="")
rcp26suitcsvout <- paste("rcp26_",2005+(1:19)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reef,ras=oa,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:19){
  x <- paste0("oa_X",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
  names(dfi) <- c("pointid","grid_code", "geometry","oa")
  dfi$suit <- with(dfi, ifelse(oa >= 3.3, 1, 0)) #if omega arag is > 3.3, it is suitable(1). otherwise not (0)
  oapt <- inner_join(reef,dfi,by="pointid")
  oadf <- as.data.frame(oapt)
  write.csv(oadf,rcp26suitcsvout[i], row.names=FALSE)
  rasterize(oapt,r008333,field="oa",fun=max,na.rm=T,filename=rcp26rout[i],overwrite=T)
  rasterize(oapt,r008333,field="suit",fun=max,na.rm=T,filename=rcp26suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/oa/rcp26",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_codex","grid_codey","oa","suit")
  b<-merge(a,filter,by="pointid")
  b$oa_coral <- b$suit * b$iscoral
  c <- count(b, vars="oa_coral")
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  count <- rbind(y,count)
} 




###rcp45

####
#extract reef suit vals
#find omega arag for each point
#classify as suitable or not - save as .csv
###

setwd("E:/oa/rcp45")

reef <- st_read("D:/distribution/all_reef_pt.shp")
oa <- brick("D:/omega_arag/co2sys_input/rcp45/OA_rcp45.nc")
rcp45rout <- paste("rcp45_",2005+(1:19)*5,"_oa.tif",sep="")
rcp45suitout <- paste("rcp45_",2005+(1:19)*5,"_suit.tif",sep="")
rcp45suitcsvout <- paste("rcp45_",2005+(1:19)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reef,ras=oa,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:19){
  x <- paste0("oa_X",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
  names(dfi) <- c("pointid","grid_code", "geometry","oa")
  dfi$suit <- with(dfi, ifelse(oa >= 3.3, 1, 0)) #if omega arag is > 3.3, it is suitable(1). otherwise not (0)
  oapt <- inner_join(reef,dfi,by="pointid")
  oadf <- as.data.frame(oapt)
  write.csv(oadf,rcp45suitcsvout[i], row.names=FALSE)
  rasterize(oapt,r008333,field="oa",fun=max,na.rm=T,filename=rcp45rout[i],overwrite=T)
  rasterize(oapt,r008333,field="suit",fun=max,na.rm=T,filename=rcp45suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/oa/rcp45",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_codex","grid_codey","oa","suit")
  b<-merge(a,filter,by="pointid")
  b$oa_coral <- b$suit * b$iscoral
  c <- count(b, vars="oa_coral")
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  count <- rbind(y,count)
} 





###rcp85

####
#extract reef suit vals
#find omega arag for each point
#classify as suitable or not - save as .csv
###

setwd("E:/oa/rcp85")

reef <- st_read("D:/distribution/all_reef_pt.shp")
oa <- brick("D:/omega_arag/co2sys_input/rcp85/OA_rcp85.nc")
rcp85rout <- paste("rcp85_",2005+(1:19)*5,"_oa.tif",sep="")
rcp85suitout <- paste("rcp85_",2005+(1:19)*5,"_suit.tif",sep="")
rcp85suitcsvout <- paste("rcp85_",2005+(1:19)*5,"_suit.csv",sep="")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reef,ras=oa,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:19){
  x <- paste0("oa_X",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "geometry",x))
  names(dfi) <- c("pointid","grid_code", "geometry","oa")
  dfi$suit <- with(dfi, ifelse(oa >= 3.3, 1, 0)) #if omega arag is > 3.3, it is suitable(1). otherwise not (0)
  oapt <- inner_join(reef,dfi,by="pointid")
  oadf <- as.data.frame(oapt)
  write.csv(oadf,rcp85suitcsvout[i], row.names=FALSE)
  rasterize(oapt,r008333,field="oa",fun=max,na.rm=T,filename=rcp85rout[i],overwrite=T)
  rasterize(oapt,r008333,field="suit",fun=max,na.rm=T,filename=rcp85suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
csvfiles <- list.files("E:/oa/rcp85",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",4)))
  names(a) <- c("pointid","grid_codex","grid_codey","oa","suit")
  b<-merge(a,filter,by="pointid")
  b$oa_coral <- b$suit * b$iscoral
  c <- count(b, vars="oa_coral")
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  count <- rbind(y,count)
} 
