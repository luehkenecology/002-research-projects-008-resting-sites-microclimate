# set working directory
setwd("C:/Rdata/Test")
#setwd("C:/Users/RenkeLuehken/Google Drive/R/Restingsites_Linda/Logger Zwischenauslegungen")
getwd()
dir()

# functions
fahr_to_kelvin <- function(temp) {
  #Converts Fahrenheit to Kelvin
  kelvin <- ((temp - 32) * (5/9)) + 273.15
  kelvin
}

kelvin_to_celsius <- function(temp) {
  #Converts Kelvin to Celsius
  Celsius <- temp - 273.15
  Celsius
}

fahr_to_celsius <- function(temp) {
  #Converts Fahrenheit to Celsius using fahr_to_kelvin() and kelvin_to_celsius()
  temp_k <- fahr_to_kelvin(temp)
  result <- kelvin_to_celsius(temp_k)
  result
}

# read libraries
library(readxl)
library(plyr)
library(gtools)
library(lubridate)
library(stringr)
library(ggplot2)

# convert *.xls to *.csv
files.to.read = list.files(pattern = "xls")
for(i in 1:length(files.to.read)){
  df = read_excel(files.to.read[i])
  write.csv(df, gsub("xls","csv",files.to.read[i]))
}

# identify name of logger
A <- strsplit(files.to.read, "[.]")
B <- do.call(rbind, A)
rest_code_uni <- B[,1]
A
# identify name of logger
g <- matrix(nrow = 0, ncol = 6)

# data manipulation
for(ff in 1:length(rest_code_uni)){
  Fahrenheit = F
  
  tryCatch({
    
    all_content <- readLines(paste(rest_code_uni[ff], ".csv", sep = ""))
    
    # exluce the first two rows (uncessary information)
    skip_second <- all_content[-c(1,2)]
    
    # convert to csv
    l <- read.csv(textConnection(skip_second), sep = ",",header = F, stringsAsFactors = T)
    
    # add logger name info
    l[,1] <- rest_code_uni[ff]
    
    # if the n of cols is too small, add a number of empty cols
    if(ncol(l) < 6){
      l[,(ncol(l)+1):6]<-NA
    }
    
    # convert potential factors to numbers
    l[,5] <- as.numeric(gsub(",",".",  l[,5]))
    l[,6] <- as.numeric(gsub(",",".",  l[,6]))
    
    # if Fahrenheit, convert to C
    if(max(l[,5], na.rm=T) > 70) {
      l[,5] <- fahr_to_celsius(l[,5])
      Fahrenheit = T
    }
    
    # progress
    print(paste(round((ff/length(rest_code_uni)*100), 1), "%, ",rest_code_uni[ff], ", F to C conversion: ", 
                Fahrenheit, sep = ""))
    
    # merge everything together
    g <- smartbind(g, l[,1:6])
    
    # print error of try function
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# convert date information to POSIXct
date1 <- paste(as.Date(g[,3], format = "%m.%d.%y"), 
               g[,4], sep = " ")
head(date1)
date2 <- as.POSIXct(strptime(date1, "%Y-%m-%d %H:%M:%S"))
head(date2)
# merge final dataset
final_data1 <- data.frame(date = date2,
                          temperature = as.numeric(g[,5]),
                          humidity = as.numeric(g[,6]),
                          id = as.factor(g[,1]),
                          height = as.factor(str_sub(g[,1], -1, -1)),
                          rs = as.factor(str_sub(g[,1], -2, -2)),
                          biotype = as.factor(str_sub(g[,1],-3,-3)),
                          site = as.factor(str_sub(g[,1], -4, -4)))

head(final_data1)
#####
ff<- final_data1
head(ff)
# convert day to julian date
ff$day <- yday(ff$date)
head(ff)

#Subset Data bezüglich Zeitraum(30.04.-31.10.2017)
#und Datenlücke (Zwischenauslesung)
zeitraum<-c("30.04.2017","31.10.2017")
tagimjahr<-strptime(zeitraum,format="%d.%m.%Y")
yday(tagimjahr)                                           
dat<-subset(ff, day > 119 & day < 305)
head(dat)
summary(dat)### 404 NAs in Temperatur, wegen Zwischenauslesung
dim(dat)
final_dat<-subset(dat, !is.na(dat$temperature))
summary(final_dat)
dim(final_dat)
head(final_dat)

#erste Plots
bla <- ddply(final_dat, .(site,id,height,rs,day), summarize, 
             temperature = mean(temperature))
plot(bla$temperature~bla$height,las=1)
plot(bla$temperature~bla$site,las=1)
abline(h=mean(bla$temperature))

###### example Plot: Deviation of daily mean temperature to reference
ref<-bla$temperature[bla$id=="JRRR"]
par(mfrow=c(1,1))
plot(bla$temperature[bla$id=="JAP2"]-ref,type="l",col="red3",
        ylim=c(-2,2),las=1,axes=F,ylab="mean deviation to reference logger in °C",
        xlab="day in the year 2017",cex.lab=1.5,main="JA Wörlitzer Elbaue")
lines(bla$temperature[bla$id=="JAP1"]-ref,type="l",col="red1")
lines(bla$temperature[bla$id=="JAP3"]-ref,type="l",col="red4")
lines(bla$temperature[bla$id=="JAV1"]-ref,type="l",col="green1")
lines(bla$temperature[bla$id=="JAV2"]-ref,type="l",col="green3")
lines(bla$temperature[bla$id=="JAV3"]-ref,type="l",col="green4")    
abline(h=0)
axis(1,at=c(0,50,100,150),labels=c("120","170","220","270"))
axis(2,las=2)
box()


# dataset without reference looger
noRef <- final_dat[!(final_dat$rs %in% c("R")),]
head(noRef)
# dataset with reference looger
Ref <- subset(final_dat, final_dat$rs == "R")
head(Ref)

##########Ab Hier habe ich nicht mehr
##########am Skript weitergearbeitet


# link dataset without reference looger and dataset with reference looger
MergData <- merge(noRef, Ref, 
                  by = c("rs", "day"),
                  all.x = T)
head(MergData)
#
LOL <- ddply(MergData, .(group, size.x,tree.x,height.x,day), summarize, 
      temp_mean = mean(temperature.x)-mean(temperature.y))

LOL2 <- ddply(LOL, .(size.x,tree.x,height.x,day), summarize, 
             temp_mean = max(temp_mean))

# data between day 120 and 250
plotData <- subset(LOL2, day > 120 & day < 250)

# plot
png(file = "all_max.png",width = 8, height=6, units = 'in', res = 500)
ggplot(na.omit(plotData), aes(x = day, y = temp_mean, 
                         group = size.x, 
                         colour = size.x))+
         geom_line()+
  ylab("mean deviation from the reference logger [°C]")+
  xlab("day of the year")+
  facet_wrap(~tree.x+height.x, scales = "free")+
  theme_classic()
dev.off()


# example for one site
for(i in 1:length(unique(MergData$group))){

  MergData_sub <- subset(MergData, group == unique(MergData$group)[i])
  
  #LOL_sub <- ddply(MergData_sub, .(size.x,tree.x,height.x,day), summarize, 
  #                 dev = mean(temperature.x, na.rm = T)-mean(temperature.y, na.rm = T))
  
  LOL_sub <- ddply(MergData_sub, .(size.x,tree.x,height.x,day), summarize, 
                   temp = mean(temperature.x)-mean(temperature.y))

  LOL_sub2 <- subset(LOL_sub, day > 120 & day < 250)

  #colnames(LOL_sub2)
  p<-ggplot(LOL_sub2, aes(y=temp,x= day, 
                      group = size.x, 
                      colour = size.x))+
    geom_line()+
    ylab("mean deviation from the reference logger [°C]")+
    xlab("day of the year")+
    facet_wrap(~tree.x+height.x)+
    theme_classic()+
    ggtitle(unique(MergData$group)[i])
  
  png(file = paste(unique(MergData$group)[i],".png", sep = ""),width = 8, height=6, units = 'in', res = 500)
  print(p)
  dev.off()
}
