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


#####
#####
###RCP85-ssp3
#####
#####

setwd("E:/suit/rcp85-ssp3")

#calculate overall suit

#load filter & points shapefile & sample raster
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
reef <- st_read("D:/distribution/all_reef_pt.shp")
r008333 <- raster("E:/bathy/gebco_2020_all.tif")



#load all rasters
#dhw, oa, storms, pop, land, sst, bathy
#check they are all same resolution and extent
#make sure all NAs -> 0's
#recalculate individual variable suits given filter



#load dhw
dhwfiles <- list.files("E:/dhw/rcp85",pattern="_suit.csv$",full.names=T)
#load land
landfiles <- list.files("E:/land_hr/ssp3-rcp85",pattern="_suit.csv$",full.names=T)
landfiles <-landfiles[-1] #we don't need 2005
#load oa
oafiles <- list.files("E:/oa/rcp85",pattern="_suit.csv$",full.names=T)
#load pop
popfiles <- list.files("E:/pop_hr/ssp3",pattern="_suit.csv$",full.names=T)
#load storms
stormfiles <- list.files("E:/storms/rcp85",pattern="_suit.csv$",full.names=T)

#name out files
csvout <- paste("rcp85_",2000+(1:19+1)*5,"_suit.csv",sep="")
pointsout <- paste("rcp85_",2000+(1:19+1)*5,"_suit",sep="")
rasterccout <- paste("rcp85_",2000+(1:19+1)*5,"_coralsuit.tif",sep="")
rasternewout <- paste("rcp85_",2000+(1:19+1)*5,"_newsuit.tif",sep="")


#load each file per year
#count overall suit sites for percentage
#csv method - multiply columns together. save csv
#join with shapefile. save
#rasterize. rotate so Pacific Ocean is centered
#print raster 
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:19){
  #dhw
  dhw <- read.csv(file=dhwfiles[i],row.names=NULL,header=F,skip=1,
                  colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",4))) #only keep pointid and suit cols
  names(dhw) <- c("pointid","dhw_suit")
  dhw[is.na(dhw)] <- 0
  #land
  land <- read.csv(file=landfiles[i],row.names=NULL,header=F,skip=1,
                  colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",2))) #only keep pointid and suit cols
  names(land) <- c("pointid","land_suit")
  #oa
  oa <- read.csv(file=oafiles[i],row.names=NULL,header=F,skip=1,
                  colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",4))) #only keep pointid and suit cols
  names(oa) <- c("pointid","oa_suit")
  oa[is.na(oa)] <- 0
  #pop
  pop <- read.csv(file=popfiles[i],row.names=NULL,header=F,skip=1,
                  colClasses=c("numeric",rep("NULL",4),"numeric",rep("NULL",2))) #only keep pointid and suit cols
  names(pop) <- c("pointid","pop_suit")
  #storms
  storms <- read.csv(file=stormfiles[i],row.names=NULL,header=F,skip=1,
                  colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",4))) #only keep pointid and suit cols
  names(storms) <- c("pointid","storm_suit")
  #join together overall suit. save csv, shapefile, and raster
  suitdf <- join_all(list(dhw,land,oa,pop,storms,filter),by="pointid")
  suitdf$suitcoral <- suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$iscoral
  suitdf$suitall <- suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$sstbathycoral
  write.csv(suitdf,csvout[i], row.names=FALSE)
  suitpt <- merge(reef,suitdf,by="pointid")
  st_write(suitpt, dsn = pointsout[i], driver = "ESRI Shapefile")
  rasterize(suitpt,r008333,field="suitcoral",fun=max,na.rm=T,filename=rasterccout[i],overwrite=T)
  rasterize(suitpt,r008333,field="suitall",fun=max,na.rm=T,filename=rasternewout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/suit/rcp85-ssp3",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
countcoral<-data.frame()
countnew<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  c <- count(a, vars="suitcoral")
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  countcoral <- rbind(y,countcoral)
  d <- count(a,vars="suitall")
  z <- add_column(d,year=paste0(i*5+2005),.before=T)
  countnew <- rbind(z,countnew)
} 






######
######
#RCP85-ssp5
######
######


setwd("E:/suit/rcp85-ssp5")

#calculate overall suit

#load filter & points shapefile & sample raster
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
reef <- st_read("D:/distribution/all_reef_pt.shp")
r008333 <- raster("E:/bathy/gebco_2020_all.tif")



#load all rasters
#dhw, oa, storms, pop, land, sst, bathy
#check they are all same resolution and extent
#make sure all NAs -> 0's
#recalculate individual variable suits given filter



#load dhw
dhwfiles <- list.files("E:/dhw/rcp85",pattern="_suit.csv$",full.names=T)
#load land
landfiles <- list.files("E:/land_hr/ssp5-rcp85",pattern="_suit.csv$",full.names=T)
landfiles <-landfiles[-1] #we don't need 2005
#load oa
oafiles <- list.files("E:/oa/rcp85",pattern="_suit.csv$",full.names=T)
#load pop
popfiles <- list.files("E:/pop_hr/ssp5",pattern="_suit.csv$",full.names=T)
#load storms
stormfiles <- list.files("E:/storms/rcp85",pattern="_suit.csv$",full.names=T)

#name out files
csvout <- paste("rcp85.5_",2000+(1:19+1)*5,"_suit.csv",sep="")
pointsout <- paste("rcp85.5_",2000+(1:19+1)*5,"_suit",sep="")
rasterccout <- paste("rcp85.5_",2000+(1:19+1)*5,"_coralsuit.tif",sep="")
rasternewout <- paste("rcp85.5_",2000+(1:19+1)*5,"_newsuit.tif",sep="")


#load each file per year
#count overall suit sites for percentage
#csv method - multiply columns together. save csv
#join with shapefile. save
#rasterize. rotate so Pacific Ocean is centered
#print raster 
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:19){
  #dhw
  dhw <- read.csv(file=dhwfiles[i],row.names=NULL,header=F,skip=1,
                  colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",4))) #only keep pointid and suit cols
  names(dhw) <- c("pointid","dhw_suit")
  dhw[is.na(dhw)] <- 0
  #land
  land <- read.csv(file=landfiles[i],row.names=NULL,header=F,skip=1,
                   colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",2))) #only keep pointid and suit cols
  names(land) <- c("pointid","land_suit")
  #oa
  oa <- read.csv(file=oafiles[i],row.names=NULL,header=F,skip=1,
                 colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",4))) #only keep pointid and suit cols
  names(oa) <- c("pointid","oa_suit")
  oa[is.na(oa)] <- 0
  #pop
  pop <- read.csv(file=popfiles[i],row.names=NULL,header=F,skip=1,
                  colClasses=c("numeric",rep("NULL",4),"numeric",rep("NULL",2))) #only keep pointid and suit cols
  names(pop) <- c("pointid","pop_suit")
  #storms
  storms <- read.csv(file=stormfiles[i],row.names=NULL,header=F,skip=1,
                     colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",4))) #only keep pointid and suit cols
  names(storms) <- c("pointid","storm_suit")
  #join together overall suit. save csv, shapefile, and raster
  suitdf <- join_all(list(dhw,land,oa,pop,storms,filter),by="pointid")
  suitdf$suitcoral <- suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$iscoral
  suitdf$suitall <- suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$sstbathycoral
  write.csv(suitdf,csvout[i], row.names=FALSE)
  suitpt <- merge(reef,suitdf,by="pointid")
  st_write(suitpt, dsn = pointsout[i], driver = "ESRI Shapefile")
  rasterize(suitpt,r008333,field="suitcoral",fun=max,na.rm=T,filename=rasterccout[i],overwrite=T)
  rasterize(suitpt,r008333,field="suitall",fun=max,na.rm=T,filename=rasternewout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/suit/rcp85-ssp5",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
countcoral<-data.frame()
countnew<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  c <- count(a, vars="suitcoral")
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  countcoral <- rbind(y,countcoral)
  d <- count(a,vars="suitall")
  z <- add_column(d,year=paste0(i*5+2005),.before=T)
  countnew <- rbind(z,countnew)
} 





#####
#####
###RCP45-ssp2
#####
#####

setwd("E:/suit/rcp45-ssp2")

#calculate overall suit

#load filter & points shapefile & sample raster
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
reef <- st_read("D:/distribution/all_reef_pt.shp")
r008333 <- raster("E:/bathy/gebco_2020_all.tif")



#load all rasters
#dhw, oa, storms, pop, land, sst, bathy
#check they are all same resolution and extent
#make sure all NAs -> 0's
#recalculate individual variable suits given filter



#load dhw
dhwfiles <- list.files("E:/dhw/rcp45",pattern="_suit.csv$",full.names=T)
#load land
landfiles <- list.files("E:/land_hr/ssp2-rcp45",pattern="_suit.csv$",full.names=T)
landfiles <-landfiles[-1] #we don't need 2005
#load oa
oafiles <- list.files("E:/oa/rcp45",pattern="_suit.csv$",full.names=T)
#load pop
popfiles <- list.files("E:/pop_hr/ssp2",pattern="_suit.csv$",full.names=T)
#load storms
stormfiles <- list.files("E:/storms/rcp45",pattern="_suit.csv$",full.names=T)

#name out files
csvout <- paste("rcp45_",2000+(1:19+1)*5,"_suit.csv",sep="")
pointsout <- paste("rcp45_",2000+(1:19+1)*5,"_suit",sep="")
rasterccout <- paste("rcp45_",2000+(1:19+1)*5,"_coralsuit.tif",sep="")
rasternewout <- paste("rcp45_",2000+(1:19+1)*5,"_newsuit.tif",sep="")


#load each file per year
#count overall suit sites for percentage
#csv method - multiply columns together. save csv
#join with shapefile. save
#rasterize. rotate so Pacific Ocean is centered
#print raster 
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:19){
  #dhw
  dhw <- read.csv(file=dhwfiles[i],row.names=NULL,header=F,skip=1,
                  colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",4))) #only keep pointid and suit cols
  names(dhw) <- c("pointid","dhw_suit")
  dhw[is.na(dhw)] <- 0
  #land
  land <- read.csv(file=landfiles[i],row.names=NULL,header=F,skip=1,
                   colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",2))) #only keep pointid and suit cols
  names(land) <- c("pointid","land_suit")
  #oa
  oa <- read.csv(file=oafiles[i],row.names=NULL,header=F,skip=1,
                 colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",4))) #only keep pointid and suit cols
  names(oa) <- c("pointid","oa_suit")
  oa[is.na(oa)] <- 0
  #pop
  pop <- read.csv(file=popfiles[i],row.names=NULL,header=F,skip=1,
                  colClasses=c("numeric",rep("NULL",4),"numeric",rep("NULL",2))) #only keep pointid and suit cols
  names(pop) <- c("pointid","pop_suit")
  #storms
  storms <- read.csv(file=stormfiles[i],row.names=NULL,header=F,skip=1,
                     colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",4))) #only keep pointid and suit cols
  names(storms) <- c("pointid","storm_suit")
  #join together overall suit. save csv, shapefile, and raster
  suitdf <- join_all(list(dhw,land,oa,pop,storms,filter),by="pointid")
  suitdf$suitcoral <- suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$iscoral
  suitdf$suitall <- suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$sstbathycoral
  write.csv(suitdf,csvout[i], row.names=FALSE)
  suitpt <- merge(reef,suitdf,by="pointid")
  st_write(suitpt, dsn = pointsout[i], driver = "ESRI Shapefile")
  rasterize(suitpt,r008333,field="suitcoral",fun=max,na.rm=T,filename=rasterccout[i],overwrite=T)
  rasterize(suitpt,r008333,field="suitall",fun=max,na.rm=T,filename=rasternewout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/suit/rcp45-ssp2",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
countcoral<-data.frame()
countnew<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  c <- count(a, vars="suitcoral")
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  countcoral <- rbind(y,countcoral)
  d <- count(a,vars="suitall")
  z <- add_column(d,year=paste0(i*5+2005),.before=T)
  countnew <- rbind(z,countnew)
} 






#####
#####
###RCP26-ssp1
#####
#####

setwd("E:/suit/rcp26-ssp1")

#calculate overall suit

#load filter & points shapefile & sample raster
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
reef <- st_read("D:/distribution/all_reef_pt.shp")
r008333 <- raster("E:/bathy/gebco_2020_all.tif")



#load all rasters
#dhw, oa, storms, pop, land, sst, bathy
#check they are all same resolution and extent
#make sure all NAs -> 0's
#recalculate individual variable suits given filter



#load dhw
dhwfiles <- list.files("E:/dhw/rcp26",pattern="_suit.csv$",full.names=T)
#load land
landfiles <- list.files("E:/land_hr/ssp1-rcp26",pattern="_suit.csv$",full.names=T)
landfiles <-landfiles[-1] #we don't need 2005
#load oa
oafiles <- list.files("E:/oa/rcp26",pattern="_suit.csv$",full.names=T)
#load pop
popfiles <- list.files("E:/pop_hr/ssp1",pattern="_suit.csv$",full.names=T)
#load storms
stormfiles <- list.files("E:/storms/rcp26",pattern="_suit.csv$",full.names=T)

#name out files
csvout <- paste("rcp26_",2000+(1:19+1)*5,"_suit.csv",sep="")
pointsout <- paste("rcp26_",2000+(1:19+1)*5,"_suit",sep="")
rasterccout <- paste("rcp26_",2000+(1:19+1)*5,"_coralsuit.tif",sep="")
rasternewout <- paste("rcp26_",2000+(1:19+1)*5,"_newsuit.tif",sep="")


#load each file per year
#count overall suit sites for percentage
#csv method - multiply columns together. save csv
#join with shapefile. save
#rasterize. rotate so Pacific Ocean is centered
#print raster 
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:19){
  #dhw
  dhw <- read.csv(file=dhwfiles[i],row.names=NULL,header=F,skip=1,
                  colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",4))) #only keep pointid and suit cols
  names(dhw) <- c("pointid","dhw_suit")
  dhw[is.na(dhw)] <- 0
  #land
  land <- read.csv(file=landfiles[i],row.names=NULL,header=F,skip=1,
                   colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",2))) #only keep pointid and suit cols
  names(land) <- c("pointid","land_suit")
  #oa
  oa <- read.csv(file=oafiles[i],row.names=NULL,header=F,skip=1,
                 colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",4))) #only keep pointid and suit cols
  names(oa) <- c("pointid","oa_suit")
  oa[is.na(oa)] <- 0
  #pop
  pop <- read.csv(file=popfiles[i],row.names=NULL,header=F,skip=1,
                  colClasses=c("numeric",rep("NULL",4),"numeric",rep("NULL",2))) #only keep pointid and suit cols
  names(pop) <- c("pointid","pop_suit")
  #storms
  storms <- read.csv(file=stormfiles[i],row.names=NULL,header=F,skip=1,
                     colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",4))) #only keep pointid and suit cols
  names(storms) <- c("pointid","storm_suit")
  #join together overall suit. save csv, shapefile, and raster
  suitdf <- join_all(list(dhw,land,oa,pop,storms,filter),by="pointid")
  suitdf$suitcoral <- suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$iscoral
  suitdf$suitall <- suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$sstbathycoral
  write.csv(suitdf,csvout[i], row.names=FALSE)
  suitpt <- merge(reef,suitdf,by="pointid")
  st_write(suitpt, dsn = pointsout[i], driver = "ESRI Shapefile")
  rasterize(suitpt,r008333,field="suitcoral",fun=max,na.rm=T,filename=rasterccout[i],overwrite=T)
  rasterize(suitpt,r008333,field="suitall",fun=max,na.rm=T,filename=rasternewout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/suit/rcp26-ssp1",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
countcoral<-data.frame()
countnew<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  c <- count(a, vars="suitcoral")
  y <- add_column(c,year=paste0(i*5+2005),.before=T)
  countcoral <- rbind(y,countcoral)
  d <- count(a,vars="suitall")
  z <- add_column(d,year=paste0(i*5+2005),.before=T)
  countnew <- rbind(z,countnew)
} 






#####
#####
###historic
#####
#####

setwd("E:/suit/hist")

#calculate overall suit

#load filter & points shapefile & sample raster
filter <- read.csv(file="E:/currentcoral/filter.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",5),"numeric","NULL","NULL","numeric"))
reef <- st_read("D:/distribution/all_reef_pt.shp")
r008333 <- raster("E:/bathy/gebco_2020_all.tif")



#load all rasters
#dhw, oa, storms, pop, land, sst, bathy
#check they are all same resolution and extent
#make sure all NAs -> 0's
#recalculate individual variable suits given filter


#use only overlapping years: 1985 to 2000
#load dhw
dhwfiles <- list.files("E:/dhw/historic",pattern="_suit.csv$",full.names=T) #from 1985 to 2015
dhwfiles <- dhwfiles[-c(5:7)] #we don't need 2005,2010 or 2015
#load land
landfiles <- list.files("E:/land_hr/historic",pattern="_suit.csv$",full.names=T) #from 1980 to 2010
landfiles <-landfiles[-c(1,6:7)] #we don't need 1980, 2005, 2010
#load oa
oafiles <- list.files("E:/oa/hist",pattern="_suit.csv$",full.names=T) #from 1980 to 2005
oafiles <- oafiles[-c(1,6)] #we don't need 1980, 2005
#load pop
popfiles <- list.files("E:/pop_hr/historic",pattern="_suit.csv$",full.names=T) #from 1970 to 2000
popfiles <- popfiles[-c(1:3)] #don't need 1970, 1975, 1980
#load storms
stormfiles <- list.files("E:/storms/historic",pattern="_suit.csv$",full.names=T) #from 1970 to 2015
stormfiles <- stormfiles[-c(1:3,8:10)] #don't need 1970,1975,1980,2005,2010,2015

#name out files
csvout <- paste("hist_",1980+(1:4)*5,"_suit.csv",sep="")
pointsout <- paste("hist",1980+(1:4)*5,"_suit",sep="")
rasterccout <- paste("hist",1980+(1:4)*5,"_coralsuit.tif",sep="")
rasternewout <- paste("hist",1980+(1:4)*5,"_newsuit.tif",sep="")


#load each file per year
#count overall suit sites for percentage
#csv method - multiply columns together. save csv
#join with shapefile. save
#rasterize. rotate so Pacific Ocean is centered
#print raster 
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:4){
  #dhw
  dhw <- read.csv(file=dhwfiles[i],row.names=NULL,header=F,skip=1,
                  colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",4))) #only keep pointid and suit cols
  names(dhw) <- c("pointid","dhw_suit")
  dhw[is.na(dhw)] <- 0
  #land
  land <- read.csv(file=landfiles[i],row.names=NULL,header=F,skip=1,
                   colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",2))) #only keep pointid and suit cols
  names(land) <- c("pointid","land_suit")
  #oa
  oa <- read.csv(file=oafiles[i],row.names=NULL,header=F,skip=1,
                 colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",4))) #only keep pointid and suit cols
  names(oa) <- c("pointid","oa_suit")
  oa[is.na(oa)] <- 0
  #pop
  pop <- read.csv(file=popfiles[i],row.names=NULL,header=F,skip=1,
                  colClasses=c("numeric",rep("NULL",4),"numeric",rep("NULL",2))) #only keep pointid and suit cols
  names(pop) <- c("pointid","pop_suit")
  #storms
  storms <- read.csv(file=stormfiles[i],row.names=NULL,header=F,skip=1,
                     colClasses=c("numeric",rep("NULL",3),"numeric",rep("NULL",4))) #only keep pointid and suit cols
  names(storms) <- c("pointid","storm_suit")
  #join together overall suit. save csv, shapefile, and raster
  suitdf <- join_all(list(dhw,land,oa,pop,storms,filter),by="pointid")
  suitdf$suitcoral <- suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$iscoral
  suitdf$suitall <- suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$sstbathycoral
  write.csv(suitdf,csvout[i], row.names=FALSE)
  suitpt <- merge(reef,suitdf,by="pointid")
  st_write(suitpt, dsn = pointsout[i], driver = "ESRI Shapefile")
  rasterize(suitpt,r008333,field="suitcoral",fun=max,na.rm=T,filename=rasterccout[i],overwrite=T)
  rasterize(suitpt,r008333,field="suitall",fun=max,na.rm=T,filename=rasternewout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/suit/hist",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
countcoral<-data.frame()
countnew<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  c <- count(a, vars="suitcoral")
  y <- add_column(c,year=paste0(i*5+1980),.before=T)
  countcoral <- rbind(y,countcoral)
  d <- count(a,vars="suitall")
  z <- add_column(d,year=paste0(i*5+1980),.before=T)
  countnew <- rbind(z,countnew)
} 