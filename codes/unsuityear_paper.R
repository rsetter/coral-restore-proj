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

#calculate year when site becomes unsuitable 
#for rcp8.5-ssp3


setwd("E:/suit/rcp85-ssp3")

#get historic csv files
csvfiles <- list.files("E:/suit/hist",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
test <- read.csv(file=csvfiles[1],row.names=NULL,header=T)
suitcoralall <- test[c("pointid","iscoral")]
suitcoralall <- suitcoralall[!(suitcoralall$iscoral==0),] #remove rows if they are not coral sites
suitcoralall <- distinct(suitcoralall) #get rid of duplicates

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="")
registerDoParallel(cl)
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  keep <- a[c("pointid","suitcoral","iscoral")]
  keep <- keep[!(keep$iscoral==0),]
  keep1 <- keep[c("pointid","suitcoral")]
  colnames(keep1) <- c("pointid",paste0("X",i*5+1980))
  keep1 <- distinct(keep1)
  suitcoralall <- left_join(suitcoralall,keep1,by="pointid")
} 
stopCluster(cl)
finish <- Sys.time()
finish-start
write.csv(suitcoralall,"E:/suit/hist_allyr1.csv", row.names=FALSE)

suitcoralall <- read.csv(file="E:/suit/hist_allyr1.csv",row.names=NULL, header=T)


#get csv files. keep suitcoral column. rename column the correct year. merge together into one dataframe
csvfiles <- list.files("E:/suit/rcp85-ssp3",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="")
registerDoParallel(cl)
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  keep <- a[c("pointid","suitcoral","iscoral")]
  keep <- keep[!(keep$iscoral==0),]
  keep1 <- keep[c("pointid","suitcoral")]
  colnames(keep1) <- c("pointid",paste0("X",i*5+2005))
  keep1 <- distinct(keep1)
  suitcoralall <- left_join(suitcoralall,keep1,by="pointid")
}
stopCluster(cl)
finish <- Sys.time()
finish-start
write.csv(suitcoralall,"E:/suit/rcp85-ssp3_allyr.csv", row.names=FALSE)

suitcoralc <- suitcoralall

#find year site becomes unsuitable

suitcoralall$year85 <- ifelse(suitcoralall$X1985 == 0, 1985,1)
suitcoralall$year90 <- ifelse(suitcoralall$year == 1, ifelse(suitcoralall$X1990 ==0,1990,1),suitcoralall$year)
suitcoralall$year95 <- ifelse(suitcoralall$year90 == 1, ifelse(suitcoralall$X1995 ==0,1995,1),suitcoralall$year90)
suitcoralall$year00 <- ifelse(suitcoralall$year95 == 1, ifelse(suitcoralall$X2000 ==0,2000,1),suitcoralall$year95)
suitcoralall$year10 <- ifelse(suitcoralall$year00 == 1, ifelse(suitcoralall$X2010 ==0,2010,1),suitcoralall$year00)
suitcoralall$year15 <- ifelse(suitcoralall$year10 == 1, ifelse(suitcoralall$X2015 ==0,2015,1),suitcoralall$year10)
suitcoralall$year20 <- ifelse(suitcoralall$year15 == 1, ifelse(suitcoralall$X2020 ==0,2020,1),suitcoralall$year15)
suitcoralall$year25 <- ifelse(suitcoralall$year20 == 1, ifelse(suitcoralall$X2025==0,2025,1),suitcoralall$year20)
suitcoralall$year30 <- ifelse(suitcoralall$year25 == 1, ifelse(suitcoralall$X2030==0,2030,1),suitcoralall$year25)
suitcoralall$year35 <- ifelse(suitcoralall$year30 == 1, ifelse(suitcoralall$X2035==0,2035,1),suitcoralall$year30)
suitcoralall$year40 <- ifelse(suitcoralall$year35 == 1, ifelse(suitcoralall$X2040==0,2040,1),suitcoralall$year35)
suitcoralall$year45 <- ifelse(suitcoralall$year40 == 1, ifelse(suitcoralall$X2045==0,2045,1),suitcoralall$year40)
suitcoralall$year50 <- ifelse(suitcoralall$year45 == 1, ifelse(suitcoralall$X2050==0,2050,1),suitcoralall$year45)
suitcoralall$year55 <- ifelse(suitcoralall$year50 == 1, ifelse(suitcoralall$X2055==0,2055,1),suitcoralall$year50)
suitcoralall$year60 <- ifelse(suitcoralall$year55 == 1, ifelse(suitcoralall$X2060==0,2060,1),suitcoralall$year55)
suitcoralall$year65 <- ifelse(suitcoralall$year60 == 1, ifelse(suitcoralall$X2065==0,2065,1),suitcoralall$year60)
suitcoralall$year70 <- ifelse(suitcoralall$year65 == 1, ifelse(suitcoralall$X2070==0,2070,1),suitcoralall$year65)
suitcoralall$year75 <- ifelse(suitcoralall$year70 == 1, ifelse(suitcoralall$X2075==0,2075,1),suitcoralall$year70)
suitcoralall$year80 <- ifelse(suitcoralall$year75 == 1, ifelse(suitcoralall$X2080==0,2080,1),suitcoralall$year75)
suitcoralall$year85 <- ifelse(suitcoralall$year80 == 1, ifelse(suitcoralall$X2085==0,2085,1),suitcoralall$year80)
suitcoralall$year90 <- ifelse(suitcoralall$year85 == 1, ifelse(suitcoralall$X2090==0,2090,1),suitcoralall$year85)
suitcoralall$year95 <- ifelse(suitcoralall$year90 == 1, ifelse(suitcoralall$X2095==0,2095,1),suitcoralall$year90)
suitcoralall$year2100 <- ifelse(suitcoralall$year95 == 1, ifelse(suitcoralall$X2100==0,2100,1),suitcoralall$year95)

suitcoralb <- suitcoralall[c("pointid","year2100")]

write.csv(suitcoralall,"E:/suit/rcp85-ssp3_yrexp.csv", row.names=FALSE)


#merge with shapefile and save
reef <- st_read("D:/distribution/all_reef_pt.shp")
suitcoralyr <- merge(reef,suitcoralb,by="pointid")
st_write(suitcoralyr, dsn = "E:/suit/rcp85-ssp3_suityr", driver = "ESRI Shapefile")



