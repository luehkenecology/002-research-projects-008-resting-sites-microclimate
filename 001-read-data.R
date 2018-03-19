# Load libraries----------------------------------------------------------------------
library(lubridate)
library(rdwd)
library(stringr)
library(sp)
#library(plyr)
#library(ggplot2)
#library(scales)

# pathogen models----------------------------------------------------------------------
# bluetongue virus
btv <- function(Temp){0.0003*Temp*(Temp-10.4057)/24}

# schmallenberg virus
sbv <- function(Temp){(0.019*(Temp-13.3))/24}

# dengue virus
denv <- function(Temp){(-0.1393+0.008*Temp)/24}

# dirofilaria immitis
diro <- function(Temp){(Temp-14/130)/24}

# malaria
malaria <- function(Temp){(0.000126*Temp+(Temp-14.244)*sqrt(34.4-Temp))/24}

# west nile virus
wnv <- function(Temp){(-0.132+0.0092*Temp)/24}

# set working directory
RPROJ <- list(PROJHOME = normalizePath(getwd()))
attach(RPROJ)
rm(RPROJ)
setwd(PROJHOME)

# read file including coordinates of the sampling sites
coordinates <- read.table (file = "data/coordinates.csv",row.names=1,header=TRUE,sep=";",fill=T) # read xy-values

# extract coordinate information
gps_info <- matrix(nrow = nrow(coordinates), ncol = 2)

# download temperature.humdity----------------------------------------------------------------------
data(fileIndex)
data(metaIndex)

# download function----------------------------------------------------------------------
dwd <- function(var_x = "wind"){
  #var_x = "air_temperature"
  
  if(var_x == "solar"){
    # file with download path
    download_path <- subset(fileIndex,
                            res == "hourly" & var == var_x)
    
    # file with gps information
    gps_info <- subset(metaIndex, res == "hourly" & var == var_x)
  }else{
    # file with download path
    download_path <- subset(fileIndex,
                            res == "hourly" & var == var_x & per == "recent")
    
    # file with gps information
    gps_info <- subset(metaIndex, res == "hourly" &
                         var == var_x & per == "recent")
    }
  
  # extract station id from the download path
  download_path$Stations_id2 <- str_sub(download_path$path, -13, -9)
  
  # convert the station ids to name with 5 characters
  for(k in 1:nrow(gps_info)){
    #i<-1
    gps_info$Stations_id2[k] <- ifelse(nchar(as.character(gps_info[k,1])) < 5, 
                                       paste(paste(replicate(5-nchar(as.character(gps_info[k,1])), "0"), collapse = ""), gps_info[k,1], sep = ""),
                                       as.character(gps_info[k,1]))}
  
  # merge download path infos and gps infos
  gps_download_merge <- merge(gps_info, download_path, by = "Stations_id2")
  
  gps_download_merge$lastyear <- as.numeric(str_sub(gps_download_merge[,4], 1,-5))
  
  gps_download_merge <- subset(gps_download_merge, lastyear == 2017)
  
  result_file <- data.frame()
  gps_info_station <- matrix(nrow = nrow(coordinates), ncol = 2)
  
  for(i in 1:nrow(coordinates)){
    # i<-2
    
    # coordinates of the sampling site
    new.pos <- c(coordinates$WGS84.Y[i],
                 coordinates$WGS84.X[i])
    
    # identify the nearest weather station
    nearest.idx <- which.min(colSums((t(gps_download_merge[,6:7]) - new.pos)^2))
    
    # temporary directory
    td = tempdir()
    
    # create temporary file
    tf = tempfile(tmpdir=td, fileext=".zip")
    
    # download the file
    download.file(paste("ftp://ftp-cdc.dwd.de/pub/CDC/observations_germany/climate/",gps_download_merge$path[nearest.idx], sep=""), tf)
    
    # unzip the file to extract the file names in the zip
    unzf <- unzip(tf,  exdir = td)
    
    # unzip the data file in the zip
    e <- unzip(tf,
               files= str_sub(unzf[length(unzf)], 49,93),
               exdir = td, overwrite=T)
    
    # identify the file.path of the unzipped data file
    fpath = file.path(td, str_sub(unzf[length(unzf)], 49,93))
    
    # read the data file  
    rdata <- read.table(fpath, sep = ";", header = T)
    
    # date to as.POSIXct
    rdata$date <- as.POSIXct(as.character(rdata$MESS_DATUM), 
                             format = "%Y%m%d%H")                
    
    # extract year info
    rdata$year <- year(rdata$date)
    
    # subset the years > 2016
    rdata_sub <- subset(rdata, rdata$year>2016)
    
    # extract gps info of station
    gps_info_station[i, 1:2] <- c(as.numeric(gps_download_merge[nearest.idx,6:7]))
    
    # save gps info of station
    write.table(gps_info_station, paste("output/", var_x, "_GPS.csv"),sep=";")
    
    # sampling site information
    rdata_sub$Standort <- coordinates$Standort[i]
    rdata_sub$ID <- coordinates$ID[i]
    rdata_sub$distance <- round(spDistsN1(pts = as.matrix(gps_download_merge[,6:7]), new.pos, longlat=T)[nearest.idx], 2)
    
    result_file <- rbind(result_file, rdata_sub)
    
    # save results
    write.table(result_file, paste("output/", var_x, ".csv"), sep = ";", row.names = F)

    # progress
    print(i)
  }
}

dwd(var_x = "air_temperature")
dwd(var_x = "solar")
dwd(var_x = "wind")
dwd(var_x = "precipitation")