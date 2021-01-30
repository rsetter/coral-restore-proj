library(raster)
library(ncdf4)
library(doParallel)
library(foreach)



######
#RCP26
######

path.rcp26 <- "D:/SST/CMIP5/RCP26/"
setwd(path.rcp26)

#use NOAA MMM file
noaa.mmm <- brick("ct5km_climatology_v3.1.nc", varname = "sst_clim_mmm", lvar=4)
crs(noaa.mmm) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "




#get mean monthly sst of all models

path.rcp26x <- "D:/SST/CMIP5/RCP26/grid/"
setwd(path.rcp26x)
files <- list.files(path.rcp26x,full.names=T) 
nfiles <- length(files)
rcp26_list  <- paste("sstm",1:nfiles,sep="")

for (i in 1:nfiles){
  
  x <- brick(files[i], varname = "tos", lvar=4)
  assign(rcp26_list[i],x)
}

sstm6 <- resample(sstm6,sstm1)
#sstm10 <- dropLayer(sstm10, c(1141:2328)) #has too many layers. goes up to 2199. There's something wrong with this raster! begins at 2099...
sstm11 <- dropLayer(sstm11, c(1141:3540)) #has too many layers. goes up to 2300

#make a monthly average of all models SST
meanmonth <- overlay(sstm1,sstm2,sstm3,sstm4,sstm5,sstm6,sstm7,sstm8,sstm9,sstm11,sstm12,sstm13,sstm14,sstm15,sstm16, 
                 fun=function(x){ mean(x,na.rm=T)}, filename = "E:/dhw/rcp26/sst26_mean_month.tif")




#median all models
median_monthly <- overlay(dhw1,dhw2,dhw3,dhw4,dhw5,dhw6,dhw7,dhw8,dhw9,dhw10,dhw11,dhw12,dhw13,dhw14,dhw15,dhw16,
                          fun=function(x){ median(x,na.rm=T)},filename = "E:/sst/sstRCP26_monthly_temp.nc") 

meanz_years <- unstack(median_monthly) # now a list of rasters 
fiveyr <- 60 # desired layers per stack -  5 years (or 60 months)
# uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
list <- c(rep(0,24),rep(1:(ceiling((length(meanz_years)-2) / fiveyr )), 
                        each = fiveyr, length.out = length(meanz_years)))
list <- list[-c(1141:1164)]
names(meanz_years) <- list
# make list of rasters into a list of stacks. basically each decade is now one separate stack
stacks <- lapply(unique(names(meanz_years)), function(y) {
  mean_years <- meanz_years[names(meanz_years) == y]
  stack(mean_years)
})
# find the mean for each 5 year period
stacksout <- lapply(stacks, function(w) {
  g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
})
meanz_fiveyr <- stack(stacksout)
names(meanz_fiveyr) <- c('2005','2010','2015','2020','2025','2030','2035','2040','2045','2050','2055','2060',
                         '2065','2070','2075','2080','2085','2090','2095','2100')
writeRaster(meanz_fiveyr, filename="E:/sst/sstRCP26med_mean5.tif", format="GTiff", overwrite=TRUE)



meanmonth <- brick("E:/sst/sstRCP26_monthly_temp.nc")

#disaggregate to match noaa.mmm
#calculate dhw per year
i_brickC <- calc(meanmonth, fun=function(x){x - 273.15}) #change from K to C
i_brickCr <- rotate(i_brickC)
i_brickCd <- disaggregate(i_brickCr,fact=2)
i_brick1 <- crop(i_brickCd,extent(180,180.5,-90,90))
i_brick2 <- crop(i_brickCd,extent(-179.5,180,-90,90))
extent(i_brick1) <- c(-180,-179.5,-90,90)
i_brickCc <- merge(i_brick1,i_brick2)
i_brickCdd <- disaggregate(i_brickCc,fact=10, filename = "E:/dhw/rcp26/sst26_mean_month_res025.tif")
DHW <- overlay(i_brickCdd, noaa.mmm, fun=function(x,y){return((x-y) * 4.34)}, filename = "E:/dhw/rcp26/sst26_DHW.tif") #gives us the weekly difference with MMM
DHWn <- reclassify(DHW, c(-Inf,0,0)) #reclassify anything less than 1degree difference doesn't count towards the DHW 

mvlist <- unstack(DHWn) # now a list of rasters 
grpsize <- 12 # desired layers per stack - 12 months or one year
# uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
names(mvlist) <- rep(1:(ceiling(length(mvlist) / grpsize )), 
                     each = grpsize, length.out = length(mvlist))
# make list of rasters into a list of stacks. basically each year is now one separate stack
stacks <- lapply(unique(names(mvlist)), function(y) {
  b <- mvlist[names(mvlist) == y]
  stack(b)
})
# find the accumulated DHW in a year - w is each stack/year. anything above 8degrees is "widespread bleaching and mortality"
stacksout <- lapply(stacks, function(w) {
  g <- calc(w, function(x){sum(x)}, forceapply = TRUE)
})
DHW_annual <- stack(stacksout)
writeRaster(DHW_annual, filename="E:/dhw/rcp26/sst26_DHW_annual.tif", overwrite=TRUE)



#calculate average annual dhw per 5 years
#unstack to groups of 5, average, then restack 
meanz_years <- unstack(DHW_annual) # now a list of rasters 
fiveyr <- 5 # desired layers per stack -  5 years (or 60 months)
# uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
list <- c(0,0,rep(1:(ceiling((length(meanz_years)-2) / fiveyr )), 
                  each = fiveyr, length.out = length(meanz_years)))
list <- list[-c(95:96)]
names(meanz_years) <- list
# make list of rasters into a list of stacks. basically each decade is now one separate stack
stacks <- lapply(unique(names(meanz_years)), function(y) {
  mean_years <- meanz_years[names(meanz_years) == y]
  stack(mean_years)
})
# find the accumulated DHW in a year - w is each stack/year. anything above 8degrees is "widespread bleaching and mortality"
stacksout <- lapply(stacks, function(w) {
  g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
})
meanz_fiveyr <- stack(stacksout)
names(meanz_fiveyr) <- c('2005','2010','2015','2020','2025','2030','2035','2040','2045','2050','2055','2060','2065','2070','2075','2080','2085','2090','2095','2100')
writeRaster(meanz_fiveyr, filename="dhw26_mean_fiveyr.tif", overwrite=TRUE)






######
#RCP45
######

path.rcp45 <- "D:/SST/CMIP5/RCP45/"
setwd(path.rcp45)

#use NOAA MMM file
noaa.mmm <- brick("D:/SST/CMIP5/RCP26/ct5km_climatology_v3.1.nc", varname = "sst_clim_mmm", lvar=4)
crs(noaa.mmm) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "




#get mean monthly sst of all models

path.rcp45x <- "D:/SST/CMIP5/RCP45/grid/"
setwd(path.rcp45x)
files <- list.files(path.rcp45x,full.names=T) 
nfiles <- length(files)
rcp45_list  <- paste("sstm",1:nfiles,sep="")

for (i in 1:nfiles){
  
  x <- brick(files[i], varname = "tos", lvar=4)
  assign(rcp45_list[i],x)
}


#make a monthly average of all models SST
meanmonth <- overlay(sstm1,sstm2,sstm3,sstm4,sstm5,sstm6,sstm7,sstm8,sstm9,sstm10,sstm11,sstm12,sstm13,sstm14,sstm15,sstm16,sstm17,sstm18,sstm19,
                     sstm20,sstm21,sstm22,sstm23,sstm24,sstm25,sstm26,sstm27,
                     fun=function(x){ mean(x,na.rm=T)}, filename = "E:/dhw/rcp45/sst45_mean_month.tif")




#median all models
median_monthly <- overlay(dhw1,dhw2,dhw3,dhw4,dhw5,dhw6,dhw7,dhw8,dhw9,dhw10,dhw11,dhw12,dhw13,dhw14,dhw15,dhw16,
                          dhw17,dhw18,dhw19,dhw20,dhw21,dhw22,dhw23,dhw24,dhw25,dhw26,dhw27, 
                          fun=function(x){ median(x,na.rm=T)},filename = "E:/sst/sstRCP45_monthly_temp.nc") 

meanz_years <- unstack(median_monthly) # now a list of rasters 
fiveyr <- 60 # desired layers per stack -  5 years (or 60 months)
# uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
list <- c(rep(0,24),rep(1:(ceiling((length(meanz_years)-2) / fiveyr )), 
                        each = fiveyr, length.out = length(meanz_years)))
list <- list[-c(1141:1164)]
names(meanz_years) <- list
# make list of rasters into a list of stacks. basically each decade is now one separate stack
stacks <- lapply(unique(names(meanz_years)), function(y) {
  mean_years <- meanz_years[names(meanz_years) == y]
  stack(mean_years)
})
# find the mean for each 5 year period
stacksout <- lapply(stacks, function(w) {
  g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
})
meanz_fiveyr <- stack(stacksout)
names(meanz_fiveyr) <- c('2005','2010','2015','2020','2025','2030','2035','2040','2045','2050','2055','2060',
                         '2065','2070','2075','2080','2085','2090','2095','2100')
writeRaster(meanz_fiveyr, filename="E:/sst/sstRCP45med_mean5.tif", format="GTiff", overwrite=TRUE)


meanmonth <- brick("E:/sst/sstRCP45_monthly_temp.nc")


#disaggregate to match noaa.mmm
#calculate dhw per year
i_brickC <- calc(meanmonth, fun=function(x){x - 273.15}) #change from K to C
i_brickCr <- rotate(i_brickC)
i_brickCd <- disaggregate(i_brickCr,fact=2)
i_brick1 <- crop(i_brickCd,extent(180,180.5,-90,90))
i_brick2 <- crop(i_brickCd,extent(-179.5,180,-90,90))
extent(i_brick1) <- c(-180,-179.5,-90,90)
i_brickCc <- merge(i_brick1,i_brick2)
i_brickCdd <- disaggregate(i_brickCc,fact=10, filename = "E:/dhw/rcp45/sst45_mean_month_res025.tif",overwrite=T)
DHW <- overlay(i_brickCdd, noaa.mmm, fun=function(x,y){return((x-y) * 4.34)}, filename = "E:/dhw/rcp45/sst45_DHW.tif") #gives us the weekly difference with MMM
DHWn = reclassify(DHW, c(-Inf,0,0)) #reclassify anything less than 1degree difference doesn't count towards the DHW 

mvlist <- unstack(DHWn) # now a list of rasters 
grpsize <- 12 # desired layers per stack - 12 months or one year
# uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
names(mvlist) <- rep(1:(ceiling(length(mvlist) / grpsize )), 
                     each = grpsize, length.out = length(mvlist))
# make list of rasters into a list of stacks. basically each year is now one separate stack
stacks <- lapply(unique(names(mvlist)), function(y) {
  b <- mvlist[names(mvlist) == y]
  stack(b)
})
# find the accumulated DHW in a year - w is each stack/year. anything above 8degrees is "widespread bleaching and mortality"
stacksout <- lapply(stacks, function(w) {
  g <- calc(w, function(x){sum(x)}, forceapply = TRUE)
})
DHW_annual <- stack(stacksout)
writeRaster(DHW_annual, filename="E:/dhw/rcp45/sst45_DHW_annual.tif", overwrite=TRUE)



#calculate average annual dhw per 5 years
#unstack to groups of 5 (e.g. 2006-2010,2011-2015, 2016-2020, 2020-2025), average, then restack 
meanz_years <- unstack(DHW_annual) # now a list of rasters 
fiveyr <- 5 # desired layers per stack -  5 years (or 60 months)
# uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
list <- c(0,0,rep(1:(ceiling((length(meanz_years)-2) / fiveyr )), 
                  each = fiveyr, length.out = length(meanz_years)))
list <- list[-c(95:96)]
names(meanz_years) <- list
# make list of rasters into a list of stacks. basically each decade is now one separate stack
stacks <- lapply(unique(names(meanz_years)), function(y) {
  mean_years <- meanz_years[names(meanz_years) == y]
  stack(mean_years)
})
# find the accumulated DHW in a year - w is each stack/year. anything above 8degrees is "widespread bleaching and mortality"
stacksout <- lapply(stacks, function(w) {
  g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
})
meanz_fiveyr <- stack(stacksout)
names(meanz_fiveyr) <- c('2005','2010','2015','2020','2025','2030','2035','2040','2045','2050','2055','2060','2065','2070','2075','2080','2085','2090','2095','2100')
writeRaster(meanz_fiveyr, filename="E:/dhw/rcp45/dhw45_mean_fiveyr.tif",overwrite=TRUE)




######

######
#RCP85
######


setwd("E:/dhw/rcp85/")

#use NOAA MMM file
noaa.mmm <- brick("D:/SST/CMIP5/RCP26/ct5km_climatology_v3.1.nc", varname = "sst_clim_mmm", lvar=4)
crs(noaa.mmm) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "




#get mean monthly sst of all models

path.rcp85x <- "D:/SST/CMIP5/RCP85/grid/"
setwd(path.rcp85x)
files <- list.files(path.rcp85x,full.names=T) 
nfiles <- length(files)
rcp85_list  <- paste("sstm",1:nfiles,sep="")

for (i in 1:nfiles){
  
  x <- brick(files[i], varname = "tos", lvar=4)
  assign(rcp85_list[i],x)
}

sstm18 <- dropLayer(sstm18, c(1141))
sstm20 <- dropLayer(sstm20, c(1141:3540))

#make a monthly average of all models SST
meanmonth <- overlay(sstm1,sstm2,sstm3,sstm4,sstm5,sstm6,sstm7,sstm8,sstm9,sstm10,sstm11,sstm12,sstm13,sstm14,sstm15,sstm16,sstm17,sstm18,sstm19,
                     sstm20,sstm21,sstm22,sstm23,sstm24,sstm25,sstm26,sstm27,
                     fun=function(x){ mean(x,na.rm=T)}, filename = "E:/dhw/rcp85/sst85_mean_month.tif")


#median all models
median_monthly <- overlay(dhw1,dhw2,dhw3,dhw4,dhw5,dhw6,dhw7,dhw8,dhw9,dhw10,dhw11,dhw12,dhw13,dhw14,dhw15,dhw16,
                          dhw17,dhw18,dhw19,dhw20,dhw21,dhw22,dhw23,dhw24,dhw25,dhw26,dhw27, 
                          fun=function(x){ median(x,na.rm=T)},filename = "E:/sstRCP85_monthly_temp.nc") 

meanz_years <- unstack(median_monthly) # now a list of rasters 
fiveyr <- 60 # desired layers per stack -  5 years (or 60 months)
# uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
list <- c(rep(0,24),rep(1:(ceiling((length(meanz_years)-2) / fiveyr )), 
                        each = fiveyr, length.out = length(meanz_years)))
list <- list[-c(1141:1164)]
names(meanz_years) <- list
# make list of rasters into a list of stacks. basically each decade is now one separate stack
stacks <- lapply(unique(names(meanz_years)), function(y) {
  mean_years <- meanz_years[names(meanz_years) == y]
  stack(mean_years)
})
# find the mean for each 5 year period
stacksout <- lapply(stacks, function(w) {
  g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
})
meanz_fiveyr <- stack(stacksout)
names(meanz_fiveyr) <- c('2005','2010','2015','2020','2025','2030','2035','2040','2045','2050','2055','2060',
                         '2065','2070','2075','2080','2085','2090','2095','2100')
writeRaster(meanz_fiveyr, filename="sstRCP85med_mean5.tif", format="GTiff", overwrite=TRUE)




start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

meanmonth <- brick("E:/sst/sstRCP85_monthly_temp.nc")

#disaggregate to match noaa.mmm
#calculate dhw per year
i_brickC <- calc(meanmonth, fun=function(x){x - 273.15}) #change from K to C
i_brickCr <- rotate(i_brickC)
i_brickCd <- disaggregate(i_brickCr,fact=2)
i_brick1 <- crop(i_brickCd,extent(180,180.5,-90,90))
i_brick2 <- crop(i_brickCd,extent(-179.5,180,-90,90))
extent(i_brick1) <- c(-180,-179.5,-90,90)
i_brickCc <- merge(i_brick1,i_brick2)
i_brickCdd <- disaggregate(i_brickCc,fact=10, filename = "E:/dhw/rcp85/sst85_mean_month_res025.tif")
DHW <- overlay(i_brickCdd, noaa.mmm, fun=function(x,y){return((x-y) * 4.34)}, filename = "E:/dhw/rcp85/sst85_DHW.tif") #gives us the weekly difference with MMM
DHWn <- reclassify(DHW, c(-Inf,0,0)) #reclassify anything less than 1degree difference doesn't count towards the DHW 

mvlist <- unstack(DHWn) # now a list of rasters 
grpsize <- 12 # desired layers per stack - 12 months or one year
# uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
names(mvlist) <- rep(1:(ceiling(length(mvlist) / grpsize )), 
                     each = grpsize, length.out = length(mvlist))
# make list of rasters into a list of stacks. basically each year is now one separate stack
stacks <- lapply(unique(names(mvlist)), function(y) {
  b <- mvlist[names(mvlist) == y]
  stack(b)
})
# find the accumulated DHW in a year - w is each stack/year. anything above 8degrees is "widespread bleaching and mortality"
stacksout <- lapply(stacks, function(w) {
  g <- calc(w, function(x){sum(x)}, forceapply = TRUE)
})
DHW_annual <- stack(stacksout)
writeRaster(DHW_annual, filename="E:/dhw/rcp85/sst85_DHW_annual.tif", overwrite=TRUE)

stopCluster(cl)
finish <- Sys.time()
finish-start


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
#calculate average annual dhw per 5 years
#unstack to groups of 5 (e.g. 2006-2010,2011-2015, 2016-2020, 2020-2025), average, then restack 
meanz_years <- unstack(DHW_annual) # now a list of rasters 
fiveyr <- 5 # desired layers per stack -  5 years (or 60 months)
# uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
list <- c(0,0,rep(1:(ceiling((length(meanz_years)-2) / fiveyr )), 
                  each = fiveyr, length.out = length(meanz_years)))
list <- list[-c(95:96)]
names(meanz_years) <- list
stacks <- lapply(unique(names(meanz_years)), function(y) {
  mean_years <- meanz_years[names(meanz_years) == y]
  stack(mean_years)
})
# find the accumulated DHW in a year - w is each stack/year. anything above 8degrees is "widespread bleaching and mortality"
stacksout <- lapply(stacks, function(w) {
  g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
})
meanz_fiveyr <- stack(stacksout)
names(meanz_fiveyr) <- c('2005','2010','2015','2020','2025','2030','2035','2040','2045','2050','2055','2060','2065','2070','2075','2080','2085','2090','2095','2100')
writeRaster(meanz_fiveyr, filename="dhw85_mean_fiveyr.nc", format="CDF", overwrite=TRUE)



stopCluster(cl)
finish <- Sys.time()
finish-start








#####
#hist empirical
#####


setwd("E:/dhw/historic")

#use NOAA MMM file
noaa.mmm <- brick("D:/SST/CMIP5/RCP26/ct5km_climatology_v3.1.nc", varname = "sst_clim_mmm", lvar=4)
crs(noaa.mmm) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "




#use oisst monthly
oisst <- brick("D:/SST/OISST/oisstmasked.nc")
oisstr <- rotate(oisst)

#disaggregate to match noaa.mmm
#calculate dhw per year - 1985,1990,1995,2000,2005,2010,2015
oisstd <- projectRaster(oisstr,noaa.mmm)
DHM <- overlay(oisstd, noaa.mmm, fun=function(x,y){return(x-y)}, filename = "E:/dhw/historic/hist_DHM.nc") #gives us the degreesC that temperature has exceeded MMM in a given month
DHW <- calc(DHM, fun=function(x){x * 4.34}, filename = "E:/dhw/historic/hist_DHW.nc") #gives us the weekly difference with MMM
DHWn = reclassify(DHW, c(-Inf,0,0)) #reclassify anything less than 1degree difference doesn't count towards the DHW 

mvlist <- unstack(DHWn) # now a list of rasters 
grpsize <- 12 # desired layers per stack - 12 months or one year
# uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
names(mvlist) <- rep(1:(ceiling(length(mvlist) / grpsize )), 
                     each = grpsize, length.out = length(mvlist))
# make list of rasters into a list of stacks. basically each year is now one separate stack
stacks <- lapply(unique(names(mvlist)), function(y) {
  b <- mvlist[names(mvlist) == y]
  stack(b)
})
# find the accumulated DHW in a year - w is each stack/year. anything above 8degrees is "widespread bleaching and mortality"
stacksout <- lapply(stacks, function(w) {
  g <- calc(w, function(x){sum(x)}, forceapply = TRUE)
})
DHW_annual <- stack(stacksout)
writeRaster(DHW_annual, filename="E:/dhw/historic/hist_DHW_annual.tif", overwrite=TRUE)
DHW_annual <- dropLayer(DHW_annual, c(8), filename="E:/dhw/historic/empirical/hist_DHW_annual.tif", overwrite=TRUE)


#
DHW_annual <-brick("E:/dhw/historic/empirical/hist_DHW_annual.tif") #this file has annual total accumulated DHW for 1981-2019
DHW_annuals <- dropLayer(DHW_annual, c(1:2,38)) #get rid of years don't need for 5 year mean: 1981, 1982, 2018, 2019

#calculate average annual dhw per 5 years
#unstack to groups of 5 (e.g. 2006-2010,2011-2015, 2016-2020, 2020-2025), average, then restack 
meanz_years <- unstack(DHW_annuals) # now a list of rasters 
fiveyr <- 5 # desired layers per stack -  years
# uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
names(meanz_years) <- rep(1:(ceiling(length(meanz_years) / fiveyr )), 
                          each = fiveyr, length.out = length(meanz_years))
# make list of rasters into a list of stacks. basically each decade is now one separate stack
stacks <- lapply(unique(names(meanz_years)), function(y) {
  mean_years <- meanz_years[names(meanz_years) == y]
  stack(mean_years)
})
# find the accumulated DHW in a year - w is each stack/year. anything above 8degrees is "widespread bleaching and mortality"
stacksout <- lapply(stacks, function(w) {
  g <- calc(w, function(x){mean(x)}, forceapply = TRUE) #or try median to balance out the 2005 bleaching event?
})
meanz_fiveyr <- stack(stacksout)
writeRaster(meanz_fiveyr, filename="hist_mean_fiveyr.tif", overwrite=TRUE)

#fiveyr_suitable <- reclassify(meanz_fiveyr, c(-Inf,8,1,8,Inf,0))
#writeRaster(fiveyr_suitable, filename="hist_suitable_fiveyr.nc", format="CDF", overwrite=TRUE)






###
##historical model
###


setwd("E:/dhw/historic/")

#use NOAA MMM file
noaa.mmm <- brick("D:/SST/CMIP5/RCP26/ct5km_climatology_v3.1.nc", varname = "sst_clim_mmm", lvar=4)
crs(noaa.mmm) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)


medmonth <- brick("E:/sst/SSThist_modelmedian.tif")

#disaggregate to match noaa.mmm
#calculate dhw per year
i_brickC <- calc(medmonth, fun=function(x){x - 273.15}) #change from K to C
i_brickCr <- rotate(i_brickC)
i_brickCd <- disaggregate(i_brickCr,fact=2)
i_brick1 <- crop(i_brickCd,extent(180,180.5,-90,90))
i_brick2 <- crop(i_brickCd,extent(-179.5,180,-90,90))
extent(i_brick1) <- c(-180,-179.5,-90,90)
i_brickCc <- merge(i_brick1,i_brick2)
i_brickCdd <- disaggregate(i_brickCc,fact=10, filename = "E:/dhw/historic/ssthist_med_month_res025.tif")
DHW <- overlay(i_brickCdd, noaa.mmm, fun=function(x,y){return((x-y)*4.34)}, filename = "E:/dhw/historic/ssthist_DHW.nc",overwrite=T) 
    #gives us the degreesK that temperature has exceeded MaxMonthlyMean in a given month
DHWn <- reclassify(DHW, c(-Inf,0,0)) #reclassify anything less than 1degree difference doesn't count towards the DHW 

mvlist <- unstack(DHWn) # now a list of rasters 
grpsize <- 12 # desired layers per stack - 12 months or one year
# uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
names(mvlist) <- rep(1:(ceiling(length(mvlist) / grpsize )), 
                     each = grpsize, length.out = length(mvlist))
# make list of rasters into a list of stacks. basically each year is now one separate stack
stacks <- lapply(unique(names(mvlist)), function(y) {
  b <- mvlist[names(mvlist) == y]
  stack(b)
})
# find the accumulated DHW in a year - w is each stack/year. anything above 8degrees is "widespread bleaching and mortality"
stacksout <- lapply(stacks, function(w) {
  g <- calc(w, function(x){sum(x)}, forceapply = TRUE)
})
DHW_annual <- stack(stacksout)
writeRaster(DHW_annual, filename="E:/dhw/historic/ssthist_DHW_annual.tif", overwrite=TRUE)

stopCluster(cl)
finish <- Sys.time()
finish-start


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
#calculate average annual dhw per 10 years
#unstack to groups of ten (one decade), average, then restack 
meanz_years <- unstack(DHW_annual) # now a list of rasters 
tenyr <- 10 # desired layers per stack -  years
# uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
names(meanz_years) <- rep(1:(ceiling(length(meanz_years) / tenyr )), 
                          each = tenyr, length.out = length(meanz_years))
# make list of rasters into a list of stacks. basically each decade is now one separate stack
stacks <- lapply(unique(names(meanz_years)), function(y) {
  mean_years <- meanz_years[names(meanz_years) == y]
  stack(mean_years)
})
# find the accumulated DHW in a year - w is each stack/year. anything above 8degrees is "widespread bleaching and mortality"
stacksout <- lapply(stacks, function(w) {
  g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
})
meanz_tenyr <- stack(stacksout)
names(meanz_tenyr) <- c(paste(1845+10*1:16, sep=""))
writeRaster(meanz_tenyr, filename="dhwhist_mean_tenyr.tif", overwrite=TRUE) #this is the mean number of dhw per year for a ten year period


stopCluster(cl)
finish <- Sys.time()
finish-start

