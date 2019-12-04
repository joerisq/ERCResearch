##################################################################################################################
## Project: Investigate Fire Occurrence Database
##
## Script purpose: Process Karen Shorts 2017 FOD Data set
##
## Date: 8th November 2019
## Author: JR
##################################################################################################################

library(RODBC)

##################################################################################################################
## Section: Set global directory locations
##################################################################################################################

# Get the data directory from readtext
ROOT_DATA_DIR <- ('C:/ERCTimeSeriesAnalysis') 
DATA_DIR <- ('Data') 
FOD_DIR <- ('FOD_Latest/RDS-2013-0009.4_ACCB/Data')
OUTPUT_DIR <- ('Output')

# specify working directory; you will need to change this line if run elsewhere
setwd(ROOT_DATA_DIR)

# Check directory
getwd()
wd <- getwd()
print(paste0('Current working dir: ', wd))

##################################################################################################################
## Section: Read FOD data from Access database
##################################################################################################################

# Check directory
initial_path <- file.path(paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/',FOD_DIR))
setwd(initial_path)
wd <- getwd()
print(paste0('Current working dir: ', wd))

# add DSN: Control Panel > Administrtative Tools > Data Sources (ODBC)
dbname <- 'FPA_FOD_20170508.accdb'
fpafod.db <- odbcConnect(dbname)
odbcGetInfo(fpafod.db)
#odbcDataSources()
sqlTables(fpafod.db, tableType="TABLE")
sqlColumns(fpafod.db, "Fires")$COLUMN_NAME

query1 <- paste("SELECT FOD_ID,NWCG_REPORTING_AGENCY,FIRE_YEAR,DISCOVERY_DATE,DISCOVERY_DOY,",
                "STAT_CAUSE_CODE,CONT_DATE,CONT_DOY,FIRE_SIZE,LATITUDE,LONGITUDE,STATE FROM Fires", sep="")

fpafod <- sqlQuery(fpafod.db, query1)
odbcClose(fpafod.db)

head(fpafod)
summary(fpafod)
length(fpafod[,1])
str(fpafod, strict.width="cut")

# Get the startday number (startdaynum day number in the year), and month (startmon) and day (startday) of each record.
library(lubridate)
fpafod$DISCOVERY_DATE <- as.Date(fpafod$DISCOVERY_DATE)
fpafod$startday <- as.numeric(format(fpafod$DISCOVERY_DATE, format="%d"))
fpafod$startmon <- as.numeric(format(fpafod$DISCOVERY_DATE, format="%m"))

# Convert acres to hectares
fpafod$AREA_HA <- fpafod$FIRE_SIZE * 0.404686

# Map data
library(maps)
plot(fpafod$LONGITUDE, fpafod$LATITUDE, ylim=c(17,80), xlim=c(-180,-55), type = "n")
map("world", add = TRUE, lwd = 2, col = "gray")
points(fpafod$LONGITUDE, fpafod$LATITUDE, pch = 16, cex=0.2, col="red")

# Time series, number of fires by diffreent causes
# 1 Lightning; 2 Equipment Use; 3 Smoking; 4 Campfire; 5 Debris Burning; 6 Railroad; 7 Arson; 8 Children; 
# 9 Miscellaneous; 10 Fireworks; 11 Power Line; 12 Structure; 13 Missing/Undefined;  
table(fpafod$STAT_CAUSE_CODE)

# Number of fires per year
hist(fpafod$FIRE_YEAR, xlim=c(1990,2015), ylim=c(0,150000),breaks=seq(1979.5,2015.5,by=1))
table(fpafod$FIRE_YEAR)

hist(fpafod$FIRE_YEAR[fpafod$STAT_CAUSE_CODE == 1], xlim=c(1990,2015), ylim=c(0,150000),
     breaks=seq(1979.5,2015.5,by=1), main="fpafod Natural")
table(fpafod$FIRE_YEAR[fpafod$STAT_CAUSE_CODE == 1])

hist(fpafod$FIRE_YEAR[fpafod$STAT_CAUSE_CODE > 1], xlim=c(1990,2015), ylim=c(0,150000),
     breaks=seq(1979.5,2015.5,by=1), main="fpafod Human")
table(fpafod$FIRE_YEAR[fpafod$STAT_CAUSE_CODE > 1])

# Get the total area burned by year
total_by_year <- tapply(fpafod$AREA_HA, fpafod$FIRE_YEAR, sum)
total_by_year

year <- as.numeric(unlist(dimnames(total_by_year)))
total_by_year <- as.numeric(total_by_year)
plot(total_by_year ~ year, pch=16, type="o", lwd=3, col="red", main="Total Area Burned (ha) (U.S. All)")

# Mean area burned by year
mean_by_year <- tapply(fpafod$AREA_HA, fpafod$FIRE_YEAR, mean)
mean_by_year

year <- as.numeric(unlist(dimnames(mean_by_year)))
mean_by_year <- as.numeric(mean_by_year)
plot(mean_by_year ~ year, pch=16, type="o", lwd=3, col="red", main="Mean Area Burned (ha) (U.S. All)")

# Start-day summaries
hist(fpafod$DISCOVERY_DOY, breaks=seq(-0.5,366.5,by=1), freq=-TRUE, ylim=c(0,12000), xlim=c(0,400))
hist(fpafod$DISCOVERY_DOY[fpafod$STAT_CAUSE_CODE == 1], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,12000), xlim=c(0,400), main="DISCOVERY_DOY (Lightning)")
hist(fpafod$DISCOVERY_DOY[fpafod$STAT_CAUSE_CODE > 1], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,12000), xlim=c(0,400), main="DISCOVERY_DOY (Human)")

# Histograms of the start day over all months
hist(fpafod$startday, breaks=seq(-0.5,31.5,by=1), freq=-TRUE, ylim=c(0,70000))

hist(fpafod$startday[fpafod$STAT_CAUSE_CODE == 1], breaks=seq(-0.5,31.5,by=1), freq=-TRUE, 
     ylim=c(0,70000), main="startday (Lighting)")
hist(fpafod$startday[fpafod$STAT_CAUSE_CODE > 1], breaks=seq(-0.5,31.5,by=1), freq=-TRUE, 
     ylim=c(0,70000),  main="startday (Human)")

# Mosaic plots  showing the distribution of fires by month and year
library(RColorBrewer)
mycols <- colors()[c(8, 5, 30, 53, 118, 72)] #
# or you could enter the color names directly
# mycols <- c("aquamarine", "antiquewhite2", "blue4", "chocolate1", "deeppink2", "cyan4")

# You could also get and store all distinct colors in the cl object and use the sample function to select colors at random
cl <- colors(distinct = TRUE)
set.seed(15887) # to set random generator seed
monthcolors <- sample(cl, 12)
mosaiccolor <- sample(cl, 11)

fpafod.tablemon <- table(fpafod$FIRE_YEAR, fpafod$startmon)
mosaicplot(fpafod.tablemon, color=monthcolors, cex.axis=0.6, las=3, main="All Fires")

fpafod.tablemon.n <- table(fpafod$FIRE_YEAR[fpafod$STAT_CAUSE_CODE == 1], 
                           fpafod$startmon[fpafod$STAT_CAUSE_CODE == 1])
mosaicplot(fpafod.tablemon.n, color=monthcolors, cex.axis=0.6, las=3, main="Natural Fires")

fpafod.tablemon.h <- table(fpafod$FIRE_YEAR[fpafod$STAT_CAUSE_CODE > 1], 
                           fpafod$startmon[fpafod$STAT_CAUSE_CODE > 1])
mosaicplot(fpafod.tablemon.h, color=monthcolors, cex.axis=0.6, las=3, main="Human Fires")

fpafod.tablecause <- table(fpafod$FIRE_YEAR, fpafod$NWCG_REPORTING_AGENCY)
mosaicplot(fpafod.tablecause, cex.axis=0.6, las=2, color=mosaiccolor, main="Fires by Organization")

# 1 Lightning; 2 Equipment Use; 3 Smoking; 4 Campfire; 5 Debris Burning; 6 Railroad; 7 Arson; 8 Children; 
# 9 Miscellaneous; 10 Fireworks; 11 Power Line; 12 Structure; 13 Missing/Undefined
mosaiccolor <- sample(cl, 13)
fpafod.tablecause <- table(fpafod$FIRE_YEAR, fpafod$STAT_CAUSE_CODE)
mosaicplot(fpafod.tablecause, cex.axis=0.6, las=2, color=mosaiccolor, main="Fires by Cause")

# 1 Lightning; 2 Equipment Use; 3 Smoking; 4 Campfire; 5 Debris Burning; 6 Railroad; 7 Arson; 8 Children; 
# 9 Miscellaneous; 10 Fireworks; 11 Power Line; 12 Structure; 13 Missing/Undefined
fpafod.tableagencycause <- table(fpafod$NWCG_REPORTING_AGENCY, fpafod$STAT_CAUSE_CODE)
mosaicplot(fpafod.tableagencycause, cex.axis=0.6, las=2, color=mosaiccolor, main="Fires by Agency and Cause")

# 1 Lightning; 2 Equipment Use; 3 Smoking; 4 Campfire; 5 Debris Burning; 6 Railroad; 7 Arson; 8 Children; 
# 9 Miscellaneous; 10 Fireworks; 11 Power Line; 12 Structure; 13 Missing/Undefined
fpafod.tablecause <- table(fpafod$FIRE_YEAR[fpafod$NWCG_REPORTING_AGENCY =="ST/C&L"], 
                           fpafod$STAT_CAUSE_CODE[fpafod$NWCG_REPORTING_AGENCY =="ST/C&L"])
mosaicplot(fpafod.tablecause, cex.axis=0.6, las=2, color=mosaiccolor, main="ST/C&L Fires by Cause")

# Cause summary histograms by cause
# 1 Lightning; 2 Equipment Use; 3 Smoking; 4 Campfire; 5 Debris Burning; 6 Railroad; 7 Arson; 8 Children; 
# 9 Miscellaneous; 10 Fireworks; 11 Power Line; 12 Structure; 13 Missing/Undefined
oldpar <- par(mfrow=c(2,2))
hist(fpafod$DISCOVERY_DOY[fpafod$STAT_CAUSE_CODE == 1], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,4000), xlim=c(0,400), xlab="", main="DISCOVERY_DOY (1 Lightning)")
hist(fpafod$DISCOVERY_DOY[fpafod$STAT_CAUSE_CODE == 2], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,4000), xlim=c(0,400), xlab="", main="DISCOVERY_DOY (2 Equipment Use)")
hist(fpafod$DISCOVERY_DOY[fpafod$STAT_CAUSE_CODE == 3], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,4000), xlim=c(0,400), xlab="", main="DISCOVERY_DOY (3 Smoking)")
hist(fpafod$DISCOVERY_DOY[fpafod$STAT_CAUSE_CODE == 4], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,4000), xlim=c(0,400), xlab="", main="DISCOVERY_DOY (4 Campfire)")

hist(fpafod$DISCOVERY_DOY[fpafod$STAT_CAUSE_CODE == 5], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,4000), xlim=c(0,400), xlab="", main="DISCOVERY_DOY (5 Debris Burning)")
hist(fpafod$DISCOVERY_DOY[fpafod$STAT_CAUSE_CODE == 6], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,4000), xlim=c(0,400), xlab="", main="DISCOVERY_DOY (6 Railroad)")
hist(fpafod$DISCOVERY_DOY[fpafod$STAT_CAUSE_CODE == 7], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,4000), xlim=c(0,400), xlab="", main="DISCOVERY_DOY (7 Arson)")
hist(fpafod$DISCOVERY_DOY[fpafod$STAT_CAUSE_CODE == 8], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,4000), xlim=c(0,400), xlab="", main="DISCOVERY_DOY (8 Children)")  

hist(fpafod$DISCOVERY_DOY[fpafod$STAT_CAUSE_CODE == 9], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,4000), xlim=c(0,400), xlab="", main="DISCOVERY_DOY (9 Miscellaneous)")
hist(fpafod$DISCOVERY_DOY[fpafod$STAT_CAUSE_CODE == 10], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,4000), xlim=c(0,400), xlab="", main="DISCOVERY_DOY (10 Fireworks)")
hist(fpafod$DISCOVERY_DOY[fpafod$STAT_CAUSE_CODE == 11], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,4000), xlim=c(0,400), xlab="", main="DISCOVERY_DOY (11 Power Line)")
hist(fpafod$DISCOVERY_DOY[fpafod$STAT_CAUSE_CODE == 12], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,4000), xlim=c(0,400), xlab="", main="DISCOVERY_DOY (12 Structure)")

hist(fpafod$DISCOVERY_DOY[fpafod$STAT_CAUSE_CODE == 13], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,4000), xlim=c(0,400), xlab="", main="DISCOVERY_DOY (13 Missing/Undefined)")
par(oldpar)

# Number of fires/year (FS and FWS)
hist(fpafod$FIRE_YEAR[fpafod$NWCG_REPORTING_AGENCY =="FS"], xlim=c(1990,2015), ylim=c(0,40000),
     breaks=seq(1979.5,2015.5,by=1))
hist(fpafod$FIRE_YEAR[fpafod$NWCG_REPORTING_AGENCY =="FWS"], xlim=c(1990,2015), ylim=c(0,40000),
     breaks=seq(1979.5,2015.5,by=1))

# Mosaic plots by agency (FS and FWS)
table(fpafod$NWCG_REPORTING_AGENCY)
fpafod.tablemon <- table(fpafod$FIRE_YEAR[fpafod$NWCG_REPORTING_AGENCY == "FS"], 
                         fpafod$startmon[fpafod$NWCG_REPORTING_AGENCY == "FS"])
mosaicplot(fpafod.tablemon, color=monthcolors, cex.axis=0.6, las=3, main="All Fires -- FS")

fpafod.tablemon <- table(fpafod$FIRE_YEAR[fpafod$NWCG_REPORTING_AGENCY == "FWS"], 
                         fpafod$startmon[fpafod$NWCG_REPORTING_AGENCY == "FWS"])
mosaicplot(fpafod.tablemon, color=monthcolors, cex.axis=0.6, las=3, main="All Fires -- FWS")

# Large(r) fires only
sizecut <- 100.
fpafod.large <- fpafod[fpafod$AREA_HA > sizecut,]
fpafod.large.tablecause <- table(fpafod.large$FIRE_YEAR, fpafod.large$NWCG_REPORTING_AGENCY)
mosaicplot(fpafod.large.tablecause, cex.axis=0.6, las=2, color=mosaiccolor, main="Large Fires by Organization")

# 1 Lightning; 2 Equipment Use; 3 Smoking; 4 Campfire; 5 Debris Burning; 6 Railroad; 7 Arson; 8 Children; 
# 9 Miscellaneous; 10 Fireworks; 11 Power Line; 12 Structure; 13 Missing/Undefined
fpafod.large.tablecause <- table(fpafod.large$FIRE_YEAR, fpafod.large$STAT_CAUSE_CODE)
mosaicplot(fpafod.large.tablecause, cex.axis=0.6, las=2, color=mosaiccolor, main="Large Fires by Cause")

# 1 Lightning; 2 Equipment Use; 3 Smoking; 4 Campfire; 5 Debris Burning; 6 Railroad; 7 Arson; 8 Children; 
# 9 Miscellaneous; 10 Fireworks; 11 Power Line; 12 Structure; 13 Missing/Undefined
fpafod.large.tableagencycause <- table(fpafod.large$NWCG_REPORTING_AGENCY, fpafod.large$STAT_CAUSE_CODE)
mosaicplot(fpafod.large.tableagencycause, cex.axis=0.6, las=2, color=mosaiccolor, 
           main="Large Fires by Agency and Cause")

# 1 Lightning; 2 Equipment Use; 3 Smoking; 4 Campfire; 5 Debris Burning; 6 Railroad; 7 Arson; 8 Children; 
# 9 Miscellaneous; 10 Fireworks; 11 Power Line; 12 Structure; 13 Missing/Undefined
fpafod.large.tablecause <- table(fpafod.large$FIRE_YEAR[fpafod.large$NWCG_REPORTING_AGENCY =="ST/C&L"], 
                                 fpafod.large$STAT_CAUSE_CODE[fpafod.large$NWCG_REPORTING_AGENCY =="ST/C&L"])
mosaicplot(fpafod.large.tablecause, cex.axis=0.6, las=2, color=mosaiccolor, main="ST/C&L Large Fires by Cause")

# Cause summary
# 1 Lightning; 2 Equipment Use; 3 Smoking; 4 Campfire; 5 Debris Burning; 6 Railroad; 7 Arson; 8 Children; 
# 9 Miscellaneous; 10 Fireworks; 11 Power Line; 12 Structure; 13 Missing/Undefined
oldpar <- par(mfrow=c(2,2))
ymax <- 500
hist(fpafod.large$DISCOVERY_DOY[fpafod.large$STAT_CAUSE_CODE == 1], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,ymax), xlim=c(0,400), xlab="", main="Large Fires DISCOVERY_DOY (1 Lightning)")
hist(fpafod.large$DISCOVERY_DOY[fpafod.large$STAT_CAUSE_CODE == 2], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,ymax), xlim=c(0,400), xlab="", main="Large Fires DISCOVERY_DOY (2 Equipment Use)")
hist(fpafod.large$DISCOVERY_DOY[fpafod.large$STAT_CAUSE_CODE == 3], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,ymax), xlim=c(0,400), xlab="", main="Large Fires DISCOVERY_DOY (3 Smoking)")
hist(fpafod.large$DISCOVERY_DOY[fpafod.large$STAT_CAUSE_CODE == 4], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,ymax), xlim=c(0,400), xlab="", main="Large Fires DISCOVERY_DOY (4 Campfire)")

hist(fpafod.large$DISCOVERY_DOY[fpafod.large$STAT_CAUSE_CODE == 5], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,ymax), xlim=c(0,400), xlab="", main="Large Fires DISCOVERY_DOY (5 Debris Burning)")
hist(fpafod.large$DISCOVERY_DOY[fpafod.large$STAT_CAUSE_CODE == 6], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,ymax), xlim=c(0,400), xlab="", main="Large Fires DISCOVERY_DOY (6 Railroad)")
hist(fpafod.large$DISCOVERY_DOY[fpafod.large$STAT_CAUSE_CODE == 7], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,ymax), xlim=c(0,400), xlab="", main="Large Fires DISCOVERY_DOY (7 Arson)")
hist(fpafod.large$DISCOVERY_DOY[fpafod.large$STAT_CAUSE_CODE == 8], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,ymax), xlim=c(0,400), xlab="", main="Large Fires DISCOVERY_DOY (8 Children)")  

hist(fpafod.large$DISCOVERY_DOY[fpafod.large$STAT_CAUSE_CODE == 9], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,ymax), xlim=c(0,400), xlab="", main="Large Fires DISCOVERY_DOY (9 Miscellaneous)")
hist(fpafod.large$DISCOVERY_DOY[fpafod.large$STAT_CAUSE_CODE == 10], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,ymax), xlim=c(0,400), xlab="", main="Large Fires DISCOVERY_DOY (10 Fireworks)")
hist(fpafod.large$DISCOVERY_DOY[fpafod.large$STAT_CAUSE_CODE == 11], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,ymax), xlim=c(0,400), xlab="", main="Large Fires DISCOVERY_DOY (11 Power Line)")
hist(fpafod.large$DISCOVERY_DOY[fpafod.large$STAT_CAUSE_CODE == 12], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,ymax), xlim=c(0,400), xlab="", main="Large Fires DISCOVERY_DOY (12 Structure)")

hist(fpafod.large$DISCOVERY_DOY[fpafod.large$STAT_CAUSE_CODE == 13], breaks=seq(-0.5,366.5,by=1), freq=-TRUE, 
     ylim=c(0,ymax), xlim=c(0,400), xlab="", main="Large Fires DISCOVERY_DOY (13 Missing/Undefined)")
par(oldpar)

proc.time()

# Read combined Canadian/US data 1982 - 2013
uscan <- read.csv(file='uscan_1986-2013.csv', header=TRUE, sep=",")
fpafod <- uscan[uscan$datasource == "fpafod", ]

# get lighting-only fires
fpafod.lt <- fpafod[fpafod$cause1 == 1 ,]
fpafod.nfires <- length(fpafod.lt[,1])
fpafod.nfires

fpafod.index <-order(fpafod.lt$area_ha, na.last=FALSE)
fpafod.lt <- fpafod.lt[fpafod.index ,]
fpafod.lt$rank <- seq(fpafod.nfires, 1)
head(cbind(fpafod.lt$area_ha,fpafod.lt$rank))

tail(cbind(fpafod.lt$area_ha,fpafod.lt$rank))

# Write a .csv files for browsing in Excel
csvpath=paste0(initial_path,'/')
write.table(fpafod.lt, paste(csvpath,"fpafod.lt.csv", sep=""), row.names=FALSE, sep=",")

# Save the R dataset
outworkspace="fpafod.lt.RData"
save.image(file=outworkspace)

# Get CNFDB data for Canada
cnfdb <- uscan[uscan$datasource == "cnfdb", ]

# get lighting-only fires
cnfdb.lt <- cnfdb[cnfdb$cause1 == 1 ,]
cnfdb.nfires <- length(cnfdb.lt[,1])
cnfdb.nfires

cnfdb.index <-order(cnfdb.lt$area_ha, na.last=FALSE)
cnfdb.lt <- cnfdb.lt[cnfdb.index ,]
cnfdb.lt$rank <- seq(cnfdb.nfires, 1)
head(cbind(cnfdb.lt$area_ha,cnfdb.lt$rank))

tail(cbind(cnfdb.lt$area_ha,cnfdb.lt$rank))

# Write a .csv files for browsing in Excel
csvpath=paste0(initial_path,'/')
write.table(cnfdb.lt, paste(csvpath,"cnfdb.lt.csv", sep=""), row.names=FALSE, sep=",")
# Save the R dataset
outworkspace="cnfdb.lt.RData"
save.image(file=outworkspace)

# Map patterns
library(maps)
options(width=110)

# Map the location of all fires
plot(0,0, ylim=c(17,80), xlim=c(-180,-55), type="n",
     xlab="longitude", ylab="latitude", main="All Lightning Fires")
map("world", add=TRUE, lwd=2, col="gray")
points(cnfdb.lt$latitude ~ cnfdb.lt$longitude, pch=16, cex=0.3, col="red")
points(fpafod.lt$latitude ~ fpafod.lt$longitude, pch=16, cex=0.3, col="red")

# NY and adjacent regions
# New York stands out, especially Long Is.:
plot(0,0, ylim=c(40,55), xlim=c(-100,-60), type="n",
       xlab="longitude", ylab="latitude", main="All Lightning Fires")
map(c("world"), add=TRUE, lwd=2, col="gray")
map(c("state"), add=TRUE, lwd=2, col="gray")
points(cnfdb.lt$latitude ~ cnfdb.lt$longitude, pch=16, cex=0.3, col="red")
points(fpafod.lt$latitude ~ fpafod.lt$longitude, pch=16, cex=0.3, col="red")

# Great Plains
# There seem to be very frew lightning fires in Kansas. There are some good things, however, in particular the prairie-forest border in Manitoba and Minnesota.
plot(0,0, ylim=c(30,55), xlim=c(-115,-90), type="n",
     xlab="longitude", ylab="latitude", main="All Lightning Fires")
map(c("world"), add=TRUE, lwd=2, col="gray")
map(c("state"), add=TRUE, lwd=2, col="gray")
points(cnfdb.lt$latitude ~ cnfdb.lt$longitude, pch=16, cex=0.3, col="red")
points(fpafod.lt$latitude ~ fpafod.lt$longitude, pch=16, cex=0.3, col="red")

# PNW Looks ok–no sign of country-border discontinuity.
plot(0,0, ylim=c(40,55), xlim=c(-130,-105), type="n",
     xlab="longitude", ylab="latitude", main="All Lightning Fires")
map("world", add=TRUE, lwd=2, col="gray")
map(c("state"), add=TRUE, lwd=2, col="gray")
points(cnfdb.lt$latitude ~ cnfdb.lt$longitude, pch=16, cex=0.3, col="red")
points(fpafod.lt$latitude ~ fpafod.lt$longitude, pch=16, cex=0.3, col="red")

# There is a hint of a little discontinuity between Alaska and Canada.
plot(0,0, ylim=c(55,70), xlim=c(-165,-135), type="n",
     xlab="longitude", ylab="latitude", main="All Lightning Fires")
map("world", add=TRUE, lwd=2, col="gray")
points(cnfdb.lt$latitude ~ cnfdb.lt$longitude, pch=16, cex=0.3, col="red")
points(fpafod.lt$latitude ~ fpafod.lt$longitude, pch=16, cex=0.3, col="red")

# Fires larger than 0.5 ha
cutsize <- 0.5
plot(0,0, ylim=c(17,80), xlim=c(-180,-55), type="n",
     xlab="longitude", ylab="latitude", main=paste("All Lightning Fires >",as.character(cutsize),"ha"))
map("world", add=TRUE, lwd=2, col="gray")
points(cnfdb.lt$latitude[cnfdb.lt$area_ha >= cutsize] ~ cnfdb.lt$longitude[cnfdb.lt$area_ha >= cutsize],
       pch=16, cex=0.3, col="red")
points(fpafod.lt$latitude[fpafod.lt$area_ha >= cutsize] ~ fpafod.lt$longitude[fpafod.lt$area_ha >= cutsize], 
       pch=16, cex=0.3, col="red")

# Fires larger than 1.0 ha
cutsize <- 1.0
plot(0,0, ylim=c(17,80), xlim=c(-180,-55), type="n",
     xlab="longitude", ylab="latitude", main=paste("All Lightning Fires >",as.character(cutsize),"ha"))
map("world", add=TRUE, lwd=2, col="gray")
points(cnfdb.lt$latitude[cnfdb.lt$area_ha >= cutsize] ~ cnfdb.lt$longitude[cnfdb.lt$area_ha >= cutsize],
       pch=16, cex=0.3, col="red")
points(fpafod.lt$latitude[fpafod.lt$area_ha >= cutsize] ~ fpafod.lt$longitude[fpafod.lt$area_ha >= cutsize], 
       pch=16, cex=0.3, col="red")

# Fires larger than 10.0 ha
cutsize <- 10.0
plot(0,0, ylim=c(17,80), xlim=c(-180,-55), type="n",
     xlab="longitude", ylab="latitude", main=paste("All Lightning Fires >",as.character(cutsize),"ha"))
map("world", add=TRUE, lwd=2, col="gray")
points(cnfdb.lt$latitude[cnfdb.lt$area_ha >= cutsize] ~ cnfdb.lt$longitude[cnfdb.lt$area_ha >= cutsize],
       pch=16, cex=0.3, col="red")
points(fpafod.lt$latitude[fpafod.lt$area_ha >= cutsize] ~ fpafod.lt$longitude[fpafod.lt$area_ha >= cutsize], 
       pch=16, cex=0.3, col="red")

# Size distributions
# The idea here is to look at the distributions (e.g. histograms) of the two data sets. 
# It’s going to turn out that the differ quite a bit. There are two possible reasons for that outcome: 
# 1) they differ becauses of unavoidable/unfixable difference in the way the data were collected, or 
# 2) they differ owing to the different environments of Canada and the U.S. (including Alaska). 
# It will also turn out that the distributions at the low end of the range of fire sizes are pretty noisy, 
# while at the upper end they are more generally similar in overall “look”.

hist(log10(fpafod.lt$area_ha), xlim=c(-4,6), ylim=c(0,120000), breaks=100,
     col="pink", border="red", main="Area (Lightning Fires)", xlab="log10 Area(ha)")
hist(log10(cnfdb.lt$area_ha), xlim=c(-4,6), ylim=c(0,120000), breaks=100,
     col=NULL, border="blue", add=TRUE)
fpafod.nfires0 <- length(fpafod.lt$area_ha); cnfdb.nfires0 <- length(cnfdb.lt$area_ha)
legend("topright", legend=c(paste("FPA-FOD",as.integer(fpafod.nfires0),"Fires"),
                            paste("CNFDB",as.integer(cnfdb.nfires0),"Fires")), lwd=2, col=c("red","blue"))

hist(log10(fpafod.lt$area_ha), xlim=c(-4,6), ylim=c(0,120000), breaks=50,
     col="pink", border="red", main="Area (Lightning Fires)", xlab="log10 Area")
hist(log10(cnfdb.lt$area_ha), xlim=c(-4,6), ylim=c(0,120000), breaks=50,
     col=NULL, border="blue", add=TRUE)
legend("topright", legend=c(paste("FPA-FOD",as.integer(fpafod.nfires0),"Fires"),
                            paste("CNFDB",as.integer(cnfdb.nfires0),"Fires")), lwd=2, col=c("red","blue"))

hist(log10(fpafod.lt$area_ha), xlim=c(-4,6), ylim=c(0,120000), breaks=200,
     col="pink", border="red", main="Area (Lightning Fires)", xlab="log10 Area")
hist(log10(cnfdb.lt$area_ha), xlim=c(-4,6), ylim=c(0,120000), breaks=200,
     col=NULL, border="blue", add=TRUE)
legend("topright", legend=c(paste("FPA-FOD",as.integer(fpafod.nfires0),"Fires"),
                            paste("CNFDB",as.integer(cnfdb.nfires0),"Fires")), lwd=2, col=c("red","blue"))

# Area-Frequency Plots
fpafod.nifires0 <- length(fpafod.lt$area_ha); cnfdb.nfires0 <- length(cnfdb.lt$area_ha)
plot(log10(fpafod.lt$area_ha), (seq(1,fpafod.nfires)), type="l", lwd=2, col="red", 
     xlim=c(-4,6), xlab="log10(Area (ha))", ylab="Cumulative Total Number of Lightning Fires", 
     main="Cumulative Total Number of Lighting Fires")
points(log10(cnfdb.lt$area_ha), (seq(1,cnfdb.nfires0)), 
       cex=0.5, pch=16, type="l", lwd=2, col="blue")
legend("topleft", legend=c(paste("FPA-FOD",as.integer(fpafod.nfires)),
                           paste("CNFDB",as.integer(cnfdb.nfires0))), lwd=2, col=c("red","blue"))

plot(log10(fpafod.lt$area_ha), (seq(1,fpafod.nfires)/fpafod.nfires), type="l", lwd=2, col="red", 
     xlim=c(-4,6), xlab="log10(Area (ha))", ylab="Cumulative Proportion", 
     main="Cumulative Proportion of Total Number of Lighting Fires")
cnfdb.nfires0 <- length(cnfdb.lt$area_ha)
points(log10(cnfdb.lt$area_ha), (seq(1,cnfdb.nfires0)/cnfdb.nfires0), 
       cex=0.5, pch=16, type="l", lwd=2, col="blue")
legend("topleft", legend=c(paste("FPA-FOD",as.integer(fpafod.nfires)),
                           paste("CNFDB",as.integer(cnfdb.nfires0))), lwd=2, col=c("red","blue"))

fpafod.cumsum <- cumsum(fpafod.lt$area_ha); cnfdb.cumsum <- cumsum(cnfdb.lt$area_ha)
plot(log10(fpafod.lt$area_ha), fpafod.cumsum, main="Cumulative Area Burned (Lightning)", ylim=c(0,8e+7),
     cex=0.5, pch=16, type="l", col="red", lwd=2, xlim=c(-4,6), xlab="log10(AREA(ha))", ylab="Cumulative Area(ha)")
points(log10(cnfdb.lt$area_ha), cnfdb.cumsum, 
       cex=0.5, pch=16, type="l", lwd=2, col="blue")
legend("topleft", legend=c(paste("FPA-FOD",as.integer(max(fpafod.cumsum)),"ha"),
                           paste("CNFDB",as.integer(max(cnfdb.cumsum)),"ha")), lwd=2, col=c("red","blue"))

fpafod.cumsum <- cumsum(fpafod.lt$area_ha); cnfdb.cumsum <- cumsum(cnfdb.lt$area_ha)
plot(log10(fpafod.lt$area_ha), (fpafod.cumsum/max(fpafod.cumsum)), 
     main="Cumulative Proportion Total Area Burned (LIghtning)",
     cex=0.5, pch=16, type="l", col="red", lwd=2, xlim=c(-4,6), 
     xlab="log10(AREA)", ylab="Cumulative Proportion")
points(log10(cnfdb.lt$area_ha), (cnfdb.cumsum/max(cnfdb.cumsum)), 
       cex=0.5, pch=16, type="l", lwd=2, col="blue")
legend("topleft", legend=c(paste("FPA-FOD",as.integer(max(fpafod.cumsum)),"ha"),
                           paste("CNFDB",as.integer(max(cnfdb.cumsum)),"ha")), lwd=2, col=c("red","blue"))

fpafod.cumsum <- cumsum(fpafod.lt$area_ha); cnfdb.cumsum <- cumsum(cnfdb.lt$area_ha)
plot(log10(fpafod.lt$area_ha), log10(fpafod.cumsum), 
     main="Log10 Cumulative Area Burned (Lightning)", ylim=c(0,8e+7),
     cex=0.5, pch=16, type="l", col="red", lwd=2, xlim=c(-4,6), xlab="log10(AREA(ha))", ylab="Cumulative Area(ha)")
points(log10(cnfdb.lt$area_ha), cnfdb.cumsum, 
       cex=0.5, pch=16, type="l", lwd=2, col="blue")
legend("topleft", legend=c(paste("FPA-FOD",as.integer(max(fpafod.cumsum)),"ha"),
                           paste("CNFDB",as.integer(max(cnfdb.cumsum)),"ha")), lwd=2, col=c("red","blue"))

#cSave the current workspace
outworkspace="all.lt.RData"
save.image(file=outworkspace)

# More distribution plots and data binning
area.thresh <- -5.0
binbeg <- -4; binend <- 6; bw <- 0.2
breaks <- seq(binbeg, binend, by=bw); mid <- seq(binbeg+bw/2, binend-bw/2, by=bw)
head(breaks); tail(breaks)
head(mid); tail(mid)

# US data, subsample and binning
# get data with areas > area.thres
farea.all <- fpafod.lt$area_ha 
fname <- "FPA-FOD Lightning" 
farea <- farea.all[farea.all > area.thresh]
farea.max <- max(farea); farea.min <- min(farea); fnfires <- length(farea)
farea.min; farea.max; fnfires; log10(farea.min); log10(farea.max)

# rank and survival curve values for clipped data
farea.index <-order(farea, na.last=FALSE)
head(farea[farea.index]);tail(farea[farea.index])

farea <- farea[farea.index]
head(farea); tail(farea)

farea.rank <- seq(fnfires, 1)
head(cbind(farea, farea.rank)); tail(cbind(farea, farea.rank))

farea.surv <- (farea.rank)/fnfires
head(cbind(farea, farea.rank, farea.surv)); tail(cbind(farea, farea.rank, farea.surv))

# bin the clipped data
#farea.hist <- hist(log10(farea), xlim=c(-4,6), breaks=100, plot=TRUE)
bin <- cut(log10(farea), breaks)
head(bin); tail(bin)

# binned variable values
fcount <- as.numeric(table(bin)); fcount

farea.ave <- as.numeric(tapply(farea, bin, mean)); farea.ave

fbin <- na.omit(data.frame(cbind(mid, 10^mid, farea.ave, fcount, farea.ave*fcount)))
names(fbin) <- c("logArea","Area.ha","aveSize","fireNum","CharFireSize")
fbin$surv <- (fnfires-cumsum(fbin$fireNum)+1)/fnfires
min.fsurv <- min(fbin$surv); max.fsurv <- max(fbin$surv)
max.fsurv; min.fsurv

inv.rel.rank <- seq(0,1, by=1/(length(fbin$surv)-1))
fbin$surv <- (fbin$surv-min.fsurv)/(max.fsurv-min.fsurv) + inv.rel.rank*min.fsurv; fbin$fsurv

head(fbin); tail(fbin)

# Canadian data, subsample and binning
# get data with areas > area.threshold
carea.all <- cnfdb.lt$area_ha[cnfdb.lt$area_ha > 0] 
cname <- "CNFDB Lightning" # 
carea <- carea.all[carea.all > area.thresh]
carea.max <- max(carea); carea.min <- min(carea); cnfires <- length(carea)
carea.min; carea.max; cnfires; log10(carea.min); log10(carea.max)

# rank and survival curve values for clipped data
carea.index <-order(carea, na.last=FALSE)
head(carea[carea.index]);tail(carea[carea.index])

carea <- carea[carea.index]
head(carea); tail(carea)

carea.rank <- seq(cnfires, 1)
head(cbind(carea, carea.rank)); tail(cbind(carea, carea.rank))

carea.surv <- (carea.rank)/cnfires
head(cbind(carea, carea.rank, carea.surv)); tail(cbind(carea, carea.rank, carea.surv))

# bin the clipped data
#carea.hist <- hist(log10(carea), xlim=c(-4,6), breaks=100, plot=TRUE)
bin <- cut(log10(carea), breaks)
head(bin); tail(bin)
# binned variable values
ccount <- as.numeric(table(bin)); ccount

carea.ave <- as.numeric(tapply(carea, bin, mean)); carea.ave

cbin <- na.omit(data.frame(cbind(mid, 10^mid, carea.ave, ccount, carea.ave*ccount)))
names(cbin) <- c("logArea","Area.ha","aveSize","fireNum","CharFireSize")
cbin$surv <- (cnfires-cumsum(cbin$fireNum)+1)/cnfires
min.csurv <- min(cbin$surv); max.csurv <- max(cbin$surv)
max.csurv; min.csurv

inv.rel.rank <- seq(0,1, by=1/(length(cbin$surv)-1))
cbin$surv <- (cbin$surv-min.csurv)/(max.csurv-min.csurv) + inv.rel.rank*min.csurv; cbin$csurv

head(cbin); tail(cbin)

# Plots of binned data
plot(log10(farea), log10(farea.surv), type="l", pch=16, xlim=c(-4,6), ylim=c(-5,0), 
     cex=0.5, lwd=3, ylab="log10(S(x))", xlab="log10(Area(ha))", col="pink", main="Impact of Binning")
points(fbin$logArea, log10(fbin$surv), type="o", pch=16, cex=0.8, lwd=2, col="red")
points(log10(carea), log10(carea.surv), type="l", pch=16, cex=0.5, lwd=3, col="lightblue")
points(cbin$logArea, log10(cbin$surv), type="o", pch=16, cex=0.8, lwd=2, col="blue")
legend("bottomleft", legend=c("FPA-FOD","FPA-FOD Binned","CNFDB","CNFDB Binned"), lwd=3, 
       col=c("pink","red","lightblue","blue"))

# Below is the “characteristic fires-size” plot of Lehsten et al.
plot(cbin$logArea, cbin$CharFireSize, xlim=c(-4,6), type="h", col="blue", lwd=3, 
     main="Characteristic Fire Sizes")
points(fbin$logArea+.05, fbin$CharFireSize, xlim=c(-4,6), type="h", col="red", lwd=3)
legend("topleft", legend=c("FPA-FOD","CNFDB"), lwd=3, col=c("red","blue"))

# The next plot is the Sachs et al. frequency-area plot, with least-squres “power-law” fits.
plot(fbin$logArea, log10(fbin$fireNum/fbin$Area.ha), xlim=c(-4,6), xlab="log10(Area(ha))",
     ylab="log10(fireNum/Area.ha)", ylim=c(-6,8), pch=16, col="red", main=paste(fname,"Frequency-Area Relationship"))
f.lm.power <- lm(log10(fbin$fireNum/fbin$Area.ha) ~ fbin$logArea)
f.lm.power2 <- lm(log10(fbin$fireNum/fbin$Area.ha)[log10(fbin$Area.ha) >=0] ~ fbin$logArea[log10(fbin$Area.ha) >=0])
abline(f.lm.power, col="pink", lwd=2); abline(f.lm.power2, col="red", lwd=2)
points(cbin$logArea, log10(cbin$fireNum/cbin$Area.ha), pch=16, col="blue")
c.lm.power <- lm(log10(cbin$fireNum/cbin$Area.ha) ~ cbin$logArea)
c.lm.power2 <- lm(log10(cbin$fireNum/cbin$Area.ha)[log10(cbin$Area.ha) >=0] ~ cbin$logArea[log10(cbin$Area.ha) >=0])
abline(c.lm.power, col="lightblue", lwd=2); abline(c.lm.power2, col="blue", lwd=2)
legend("bottomleft", legend=c(paste("FPA-FOD slope=",round(f.lm.power$coefficients[2], digits=4)), 
                              paste("FPA-FOD > 1.0 ha, slope=",round(f.lm.power2$coefficients[2], digits=4)),
                              paste("CNFDB slope=",round(c.lm.power$coefficients[2], digits=4)),
                              paste("CNFDB > 1.0 ha, slope=",round(c.lm.power2$coefficients[2], digits=4))),
       lwd=2, col=c("pink","red","lightblue","blue"))

# US Data Rank and survival-curve data
# Note that this overwrites previous data…
area.thresh = 1.0 # ha
fpafod.lt.area <- fpafod.lt$area_ha[fpafod.lt$area_ha > area.thresh]
fpafod.lt.name <- "Merged U.S. Lightning Data"
fpafod.lt.area.max <- max(fpafod.lt.area); fpafod.lt.area.min <- min(fpafod.lt.area)
fpafod.lt.nfires <- length(fpafod.lt.area)
cbind(fpafod.lt.area.max, fpafod.lt.area.min, fpafod.lt.nfires)

fpafod.lt.area.index <- order(fpafod.lt.area, na.last=FALSE)
head(fpafod.lt.area.index); tail(fpafod.lt.area.index)

# Sort the area data, and get survival-curve ordinates.
fpafod.lt.area.sort <- fpafod.lt.area[fpafod.lt.area.index]
fpafod.lt.area.rank <- seq(fpafod.lt.nfires, 1)
fpafod.lt.area.surv <- fpafod.lt.area.rank/fpafod.lt.nfires
head(cbind(fpafod.lt.area.sort, fpafod.lt.area.rank, fpafod.lt.area.surv))

tail(cbind(fpafod.lt.area.sort, fpafod.lt.area.rank, fpafod.lt.area.surv))

# Initial fits
# Initital values of parameters
a0 <- area.thresh
beta0 <- 0.40  
theta0 <-  30000.0 

plot(log10(fpafod.lt.area.sort), log10(fpafod.lt.area.surv), 
     type="l", pch=16, cex=0.5, lwd=2, col="red", xlim=c(-1,6), ylim=c(-6,0),
     xlab="log10(Area)", ylab="log10(S(x)", main=paste(fpafod.lt.name,"Fires >",as.character(a0),"ha"))
abline(0, -0.5, col="gray", lwd=2)

fpafod.lt.fx.pareto <- 1.0 - (a0/fpafod.lt.area.sort)^beta0
head(fpafod.lt.fx.pareto); tail(fpafod.lt.fx.pareto)

points(log10(fpafod.lt.area.sort), log10(1.0 - fpafod.lt.fx.pareto), col="blue", type="l", lwd=2)

fpafod.lt.fx.tappareto <- 1.0 - ((a0/fpafod.lt.area.sort)^beta0)*exp((a0-fpafod.lt.area.sort)/theta0)
head(fpafod.lt.fx.tappareto); tail(fpafod.lt.fx.tappareto)

points(log10(fpafod.lt.area.sort), log10(1.0 - fpafod.lt.fx.tappareto), type="l", pch=16, cex=0.5, col="purple", lwd=2)  
legend("bottomleft", title="by hand", legend=c("Observed","Pareto","Tapered Parto",
                                               "slope = -0.5"), lwd=2, col=c("red","blue","purple","gray"))

# Optimization
# optim() estimation Pareto
SSE.Pareto <- function(p, data, surv){
  y <- log10(surv)
  yhat <- log10((p[2]/data)^p[1])
  x <- sum((y-yhat)^2.0)
  return(x)
}
p0 <- c(beta0, a0); p0

fit.Pareto <- optim(p0, SSE.Pareto, data=fpafod.lt.area.sort, surv=fpafod.lt.area.surv )
beta1 <- fit.Pareto$par[1]; a1 <- fit.Pareto$par[2]; sse1 <- fit.Pareto$value
beta1; a1; sse1

# optim() estimation tapered Pareto
SSE.tapPareto <- function(p, data, surv){
  y <- log10(surv)
  yhat <- log10(((p[3]/data)^p[1])*exp((p[3]-data)/p[2]))
  x <- sum((y-yhat)^2.0)
  return(x)
}
p2 <- c(beta0, theta0, a0); p2

fit.tapPareto1 <- optim(p2, SSE.tapPareto, data=fpafod.lt.area.sort, surv=fpafod.lt.area.surv)
beta2 <- fit.tapPareto1$par[1]; theta2 <- fit.tapPareto1$par[2]
a2 <- fit.tapPareto1$par[3]; sse2 <- fit.tapPareto1$value
beta2; theta2; a2; sse2

# optim() estimation tapered Pareto -- a fixed
SSE.tapPareto2 <- function(p, data, surv, a){
  y <- log10(surv)
  yhat <- log10( ((a/data)^p[1])*exp((a-data)/p[2]) )
  x <- sum((y-yhat)^2.0)
  return(x)
}
p3 <- c(beta0, theta0); p3; a <- a0

fit.tapPareto2 <- optim(p3, SSE.tapPareto2, data=fpafod.lt.area.sort, a=a, surv=fpafod.lt.area.surv)
beta3 <- fit.tapPareto2$par[1]; theta3 <- fit.tapPareto2$par[2]
a3 <- a0; sse3 <- fit.tapPareto2$value
beta3; theta3; a3; sse3

# compare observed and fitted
plot(log10(fpafod.lt.area.sort), log10(fpafod.lt.area.surv), type="n", ylim=c(-5,0), xlim=c(-1,6),
     xlab="log10(Area)", ylab="log10(S(x)", main=paste(fpafod.lt.name,"Fires >",as.character(a0),"ha"))
abline(log10(fpafod.lt.area.min)/2, -0.5, col="gray", lwd=2)

fpafod.lt.fx.pareto <- 1.0 - (a1/fpafod.lt.area.sort)^beta1
head(fpafod.lt.fx.pareto); tail(fpafod.lt.fx.pareto)

points(log10(fpafod.lt.area.sort), log10(1.0 - fpafod.lt.fx.pareto), col="blue", type="l", lwd=2)

fpafod.lt.fx.tappareto <- 1.0 - ((a2/fpafod.lt.area.sort)^beta2)*exp((a2-fpafod.lt.area.sort)/theta2)
head(fpafod.lt.fx.tappareto); tail(fpafod.lt.fx.tappareto)

points(log10(fpafod.lt.area.sort), log10(1.0 - fpafod.lt.fx.tappareto), col="purple", type="l", lwd=2) 

fpafod.lt.fx.tappareto <- 1.0 - ((a3/fpafod.lt.area.sort)^beta3)*exp((a3-fpafod.lt.area.sort)/theta3)
head(fpafod.lt.fx.tappareto); tail(fpafod.lt.fx.tappareto)

points(log10(fpafod.lt.area.sort), log10(1.0 - fpafod.lt.fx.tappareto), col="pink", type="l", lwd=2) 

points(log10(fpafod.lt.area.sort), log10(fpafod.lt.area.surv), pch=16, cex=0.5, col="red")

legend("bottomleft", title="optim() fit", legend=c("FPA-FOD","Pareto","Tapered Pareto",
                                                   "Tapered Pareto (a fixed)", "slope = -0.5"), lwd=2, col=c("red","blue","purple","pink","gray"))

print(cbind(beta0, theta0, a0, beta1, a1, sse1, beta2, theta2, a2, sse2, beta3, theta3, a3, sse3))

us.beta0 <- beta0; us.theta0 <- theta0; us.a0 <- a0; us.beta1 <- beta1; us.a1 <- a1; 
us.beta2 <- beta2; us.theta2 <- theta3; us.a2 <- a2;
us.beta3 <- beta3; us.theta3 <- theta3; us.a3 <- a3

# QQ Plots
library(PtProcess)
fpafod.lt.area.prob <- ppoints(fpafod.lt.area.sort)
qp <- qpareto(fpafod.lt.area.prob, beta1, a1, lower.tail=TRUE, log.p=FALSE)
head(qp); tail(qp)

qtap <- qtappareto(fpafod.lt.area.prob, beta2, theta2, a2, lower.tail=TRUE, log.p=FALSE, tol=1e-8)
head(qtap); tail(qtap)

head(cbind(log10(fpafod.lt.area.sort),qp,qtap,log10(qp),log10(qtap),fpafod.lt.area.prob))

tail(cbind(log10(fpafod.lt.area.sort),qp,qtap,log10(qp),log10(qtap),fpafod.lt.area.prob))

plot(log10(fpafod.lt.area.sort), log10(qtap), xlim=c(-1,6), ylim=c(-1,6), xlab="Observed log10(Area)", type="n", 
     ylab="log10(Theoretical Quantiles)", main=paste("QQPlot",fpafod.lt.name,"Fires >",as.character(a0),"ha"))
abline(0.0, 1.0, col="gray", lwd=2)
points(log10(fpafod.lt.area.sort),log10(qp), type="l", col="blue", lwd=2)
points(log10(fpafod.lt.area.sort),log10(qtap), type="l", col="purple", lwd=2)
legend("topleft", legend=c("Pareto","Tapered Pareto"), lwd=2, col=c("blue","purple"))

# Canadian data Rank and survival-curve data
# Note that this overwrites previous data…
area.thresh = 1.0 # ha
cnfdb.lt.area <- cnfdb.lt$area_ha[cnfdb.lt$area_ha > area.thresh]
cnfdb.lt.name <- "Canada Lightning Data"
cnfdb.lt.area.max <- max(cnfdb.lt.area); cnfdb.lt.area.min <- min(cnfdb.lt.area)
cnfdb.lt.nfires <- length(cnfdb.lt.area)
cbind(cnfdb.lt.area.max, cnfdb.lt.area.min, cnfdb.lt.nfires)

cnfdb.lt.area.index <- order(cnfdb.lt.area, na.last=FALSE)
head(cnfdb.lt.area.index); tail(cnfdb.lt.area.index)

# Sort the area data, and get survival-curve ordinates.
cnfdb.lt.area.sort <- cnfdb.lt.area[cnfdb.lt.area.index]
cnfdb.lt.area.rank <- seq(cnfdb.lt.nfires, 1)
cnfdb.lt.area.surv <- cnfdb.lt.area.rank/cnfdb.lt.nfires
head(cbind(cnfdb.lt.area.sort, cnfdb.lt.area.rank, cnfdb.lt.area.surv))

tail(cbind(cnfdb.lt.area.sort, cnfdb.lt.area.rank, cnfdb.lt.area.surv))

# Initial fits
# intital values of parameters
a0 <- area.thresh
beta0 <- 0.4 
theta0 <-  30000.0 

plot(log10(cnfdb.lt.area.sort), log10(cnfdb.lt.area.surv), 
     type="l", pch=16, cex=0.5, lwd=2, col="red", xlim=c(-1,6), ylim=c(-6,0),
     xlab="log10(Area)", ylab="log10(S(x)", main=paste(cnfdb.lt.name,"Fires >",as.character(a0),"ha"))
abline(0, -0.5, col="gray", lwd=2)

cnfdb.lt.fx.pareto <- 1.0 - (a0/cnfdb.lt.area.sort)^beta0
head(cnfdb.lt.fx.pareto); tail(cnfdb.lt.fx.pareto)

points(log10(cnfdb.lt.area.sort), log10(1.0 - cnfdb.lt.fx.pareto), col="blue", type="l", lwd=2)

cnfdb.lt.fx.tappareto <- 1.0 - ((a0/cnfdb.lt.area.sort)^beta0)*exp((a0-cnfdb.lt.area.sort)/theta0)
head(cnfdb.lt.fx.tappareto); tail(cnfdb.lt.fx.tappareto)

points(log10(cnfdb.lt.area.sort), log10(1.0 - cnfdb.lt.fx.tappareto), type="l", pch=16, cex=0.5, col="purple", lwd=2)  
legend("bottomleft", title="by hand", legend=c("Observed","Pareto","Tapered Parto",
                                               "slope = -0.5"), lwd=2, col=c("red","blue","purple","gray"))
# Optimization
# optim() estimation Pareto
SSE.Pareto <- function(p, data, surv){
  y <- log10(surv)
  yhat <- log10((p[2]/data)^p[1])
  x <- sum((y-yhat)^2.0)
  return(x)
}
p0 <- c(beta0, a0); p0

fit.Pareto <- optim(p0, SSE.Pareto, data=cnfdb.lt.area.sort, surv=cnfdb.lt.area.surv )
beta1 <- fit.Pareto$par[1]; a1 <- fit.Pareto$par[2]; sse1 <- fit.Pareto$value
beta1; a1; sse1

# optim() estimation tapered Pareto
SSE.tapPareto <- function(p, data, surv){
  y <- log10(surv)
  yhat <- log10(((p[3]/data)^p[1])*exp((p[3]-data)/p[2]))
  x <- sum((y-yhat)^2.0)
  return(x)
}
p2 <- c(beta0, theta0, a0); p2

fit.tapPareto1 <- optim(p2, SSE.tapPareto, data=cnfdb.lt.area.sort, surv=cnfdb.lt.area.surv)
beta2 <- fit.tapPareto1$par[1]; theta2 <- fit.tapPareto1$par[2]
a2 <- fit.tapPareto1$par[3]; sse2 <- fit.tapPareto1$value
beta2; theta2; a2; sse2

# optim() estimation tapered Pareto -- a fixed
SSE.tapPareto2 <- function(p, data, surv, a){
  y <- log10(surv)
  yhat <- log10( ((a/data)^p[1])*exp((a-data)/p[2]) )
  x <- sum((y-yhat)^2.0)
  return(x)
}
p3 <- c(beta0, theta0); p3; a <- a0

fit.tapPareto2 <- optim(p3, SSE.tapPareto2, data=cnfdb.lt.area.sort, a=a, surv=cnfdb.lt.area.surv)
beta3 <- fit.tapPareto2$par[1]; theta3 <- fit.tapPareto2$par[2]
a3 <- a0; sse3 <- fit.tapPareto2$value
beta3; theta3; a3; sse3

# compare observed and fitted
plot(log10(cnfdb.lt.area.sort), log10(cnfdb.lt.area.surv), type="n", ylim=c(-5,0), xlim=c(-1,6),
     xlab="log10(Area)", ylab="log10(S(x)", main=paste(cnfdb.lt.name,"Fires >",as.character(a0),"ha"))
abline(log10(cnfdb.lt.area.min)/2, -0.5, col="gray", lwd=2)

cnfdb.lt.fx.pareto <- 1.0 - (a1/cnfdb.lt.area.sort)^beta1
head(cnfdb.lt.fx.pareto); tail(cnfdb.lt.fx.pareto)

points(log10(cnfdb.lt.area.sort), log10(1.0 - cnfdb.lt.fx.pareto), col="blue", type="l", lwd=2)

cnfdb.lt.fx.tappareto <- 1.0 - ((a2/cnfdb.lt.area.sort)^beta2)*exp((a2-cnfdb.lt.area.sort)/theta2)
head(cnfdb.lt.fx.tappareto); tail(cnfdb.lt.fx.tappareto)

points(log10(cnfdb.lt.area.sort), log10(1.0 - cnfdb.lt.fx.tappareto), col="purple", type="l", lwd=2) 

cnfdb.lt.fx.tappareto <- 1.0 - ((a3/cnfdb.lt.area.sort)^beta3)*exp((a3-cnfdb.lt.area.sort)/theta3)
head(cnfdb.lt.fx.tappareto); tail(cnfdb.lt.fx.tappareto)

points(log10(cnfdb.lt.area.sort), log10(1.0 - cnfdb.lt.fx.tappareto), col="pink", type="l", lwd=2) 

points(log10(cnfdb.lt.area.sort), log10(cnfdb.lt.area.surv), pch=16, cex=0.5, col="red")

legend("bottomleft", title="optim() fit", legend=c("FPA-FOD","Pareto","Tapered Pareto",
                                                   "Tapered Pareto (a fixed)", "slope = -0.5"), lwd=2, col=c("red","blue","purple","pink","gray"))
print(cbind(beta0, theta0, a0, beta1, a1, sse1, beta2, theta2, a2, sse2, beta3, theta3, a3, sse3))

can.beta0 <- beta0; can.theta0 <- theta0; can.a0 <- a0; can.beta1 <- beta1; can.a1 <- a1; 
can.beta2 <- beta2; can.theta2 <- theta3; can.a2 <- a2;
can.beta3 <- beta3; can.theta3 <- theta3; can.a3 <- a3

# QQ Plots
library(PtProcess)
cnfdb.lt.area.prob <- ppoints(cnfdb.lt.area.sort)
qp <- qpareto(cnfdb.lt.area.prob, beta1, a1, lower.tail=TRUE, log.p=FALSE)
head(qp); tail(qp)

qtap <- qtappareto(cnfdb.lt.area.prob, beta2, theta2, a2, lower.tail=TRUE, log.p=FALSE, tol=1e-8)
head(qtap); tail(qtap)

head(cbind(log10(cnfdb.lt.area.sort),qp,qtap,log10(qp),log10(qtap),cnfdb.lt.area.prob))

tail(cbind(log10(cnfdb.lt.area.sort),qp,qtap,log10(qp),log10(qtap),cnfdb.lt.area.prob))

plot(log10(cnfdb.lt.area.sort), log10(qtap), xlim=c(-1,6), ylim=c(-1,6), xlab="Observed log10(Area)", type="n", 
     ylab="log10(Theoretical Quantiles)", main=paste("QQPlot",cnfdb.lt.name,"Fires >",as.character(a0),"ha"))
abline(0.0, 1.0, col="gray", lwd=2)
points(log10(cnfdb.lt.area.sort),log10(qp), type="l", col="blue", lwd=2)
points(log10(cnfdb.lt.area.sort),log10(qtap), type="l", col="purple", lwd=2)
legend("topleft", legend=c("Pareto","Tapered Pareto"), lwd=2, col=c("blue","purple"))

# US & Canada Data
uscan.lt <- uscan[uscan$cause1 == 1,]

# Rank and survival-curve data
# Note that this overwrites previous data…
area.thresh = 1.0 # ha
uscan.lt.area <- uscan.lt$area_ha[uscan.lt$area_ha > area.thresh]
uscan.lt.name <- "U.S. & Canada Lightning Data"
uscan.lt.area.max <- max(uscan.lt.area); uscan.lt.area.min <- min(uscan.lt.area)
uscan.lt.nfires <- length(uscan.lt.area)
cbind(uscan.lt.area.max, uscan.lt.area.min, uscan.lt.nfires)

uscan.lt.area.index <- order(uscan.lt.area, na.last=FALSE)
head(uscan.lt.area.index); tail(uscan.lt.area.index)

# Sort the area data, and get survival-curve ordinates.
uscan.lt.area.sort <- uscan.lt.area[uscan.lt.area.index]
uscan.lt.area.rank <- seq(uscan.lt.nfires, 1)
uscan.lt.area.surv <- uscan.lt.area.rank/uscan.lt.nfires
head(cbind(uscan.lt.area.sort, uscan.lt.area.rank, uscan.lt.area.surv))

tail(cbind(uscan.lt.area.sort, uscan.lt.area.rank, uscan.lt.area.surv))

# Initial fits
# Initital values of parameters
a0 <- area.thresh
beta0 <- 0.40  
theta0 <-  30000.0 

plot(log10(uscan.lt.area.sort), log10(uscan.lt.area.surv), 
     type="l", pch=16, cex=0.5, lwd=2, col="red", xlim=c(-1,6), ylim=c(-6,0),
     xlab="log10(Area)", ylab="log10(S(x)", main=paste(uscan.lt.name,"Fires >",as.character(a0),"ha"))
abline(0, -0.5, col="gray", lwd=2)

uscan.lt.fx.pareto <- 1.0 - (a0/uscan.lt.area.sort)^beta0
head(uscan.lt.fx.pareto); tail(uscan.lt.fx.pareto)

points(log10(uscan.lt.area.sort), log10(1.0 - uscan.lt.fx.pareto), col="blue", type="l", lwd=2)

uscan.lt.fx.tappareto <- 1.0 - ((a0/uscan.lt.area.sort)^beta0)*exp((a0-uscan.lt.area.sort)/theta0)
head(uscan.lt.fx.tappareto); tail(uscan.lt.fx.tappareto)

points(log10(uscan.lt.area.sort), log10(1.0 - uscan.lt.fx.tappareto), type="l", pch=16, cex=0.5, col="purple", lwd=2)  
legend("bottomleft", title="by hand", legend=c("Observed","Pareto","Tapered Parto",
                                               "slope = -0.5"), lwd=2, col=c("red","blue","purple","gray"))
# Optimization
# optim() estimation Pareto
SSE.Pareto <- function(p, data, surv){
  y <- log10(surv)
  yhat <- log10((p[2]/data)^p[1])
  x <- sum((y-yhat)^2.0)
  return(x)
}
p0 <- c(beta0, a0); p0

fit.Pareto <- optim(p0, SSE.Pareto, data=uscan.lt.area.sort, surv=uscan.lt.area.surv )
beta1 <- fit.Pareto$par[1]; a1 <- fit.Pareto$par[2]; sse1 <- fit.Pareto$value
beta1; a1; sse1

# optim() estimation tapered Pareto
SSE.tapPareto <- function(p, data, surv){
  y <- log10(surv)
  yhat <- log10(((p[3]/data)^p[1])*exp((p[3]-data)/p[2]))
  x <- sum((y-yhat)^2.0)
  return(x)
}
p2 <- c(beta0, theta0, a0); p2

fit.tapPareto1 <- optim(p2, SSE.tapPareto, data=uscan.lt.area.sort, surv=uscan.lt.area.surv)
beta2 <- fit.tapPareto1$par[1]; theta2 <- fit.tapPareto1$par[2]
a2 <- fit.tapPareto1$par[3]; sse2 <- fit.tapPareto1$value
beta2; theta2; a2; sse2

# optim() estimation tapered Pareto -- a fixed
SSE.tapPareto2 <- function(p, data, surv, a){
  y <- log10(surv)
  yhat <- log10( ((a/data)^p[1])*exp((a-data)/p[2]) )
  x <- sum((y-yhat)^2.0)
  return(x)
}
p3 <- c(beta0, theta0); p3; a <- a0

fit.tapPareto2 <- optim(p3, SSE.tapPareto2, data=uscan.lt.area.sort, a=a, surv=uscan.lt.area.surv)
beta3 <- fit.tapPareto2$par[1]; theta3 <- fit.tapPareto2$par[2]
a3 <- a0; sse3 <- fit.tapPareto2$value
beta3; theta3; a3; sse3

# compare observed and fitted
plot(log10(uscan.lt.area.sort), log10(uscan.lt.area.surv), type="n", ylim=c(-5,0), xlim=c(-1,6),
     xlab="log10(Area)", ylab="log10(S(x)", main=paste(uscan.lt.name,"Fires >",as.character(a0),"ha"))
abline(log10(uscan.lt.area.min)/2, -0.5, col="gray", lwd=2)

uscan.lt.fx.pareto <- 1.0 - (a1/uscan.lt.area.sort)^beta1
head(uscan.lt.fx.pareto); tail(uscan.lt.fx.pareto)

points(log10(uscan.lt.area.sort), log10(1.0 - uscan.lt.fx.pareto), col="blue", type="l", lwd=2)

uscan.lt.fx.tappareto <- 1.0 - ((a2/uscan.lt.area.sort)^beta2)*exp((a2-uscan.lt.area.sort)/theta2)
head(uscan.lt.fx.tappareto); tail(uscan.lt.fx.tappareto)

points(log10(uscan.lt.area.sort), log10(1.0 - uscan.lt.fx.tappareto), col="purple", type="l", lwd=2) 

uscan.lt.fx.tappareto <- 1.0 - ((a3/uscan.lt.area.sort)^beta3)*exp((a3-uscan.lt.area.sort)/theta3)
head(uscan.lt.fx.tappareto); tail(uscan.lt.fx.tappareto)

points(log10(uscan.lt.area.sort), log10(1.0 - uscan.lt.fx.tappareto), col="pink", type="l", lwd=2) 

points(log10(uscan.lt.area.sort), log10(uscan.lt.area.surv), pch=16, cex=0.5, col="red")

legend("bottomleft", title="optim() fit", legend=c("FPA-FOD","Pareto","Tapered Pareto",
                                                   "Tapered Pareto (a fixed)", "slope = -0.5"), lwd=2, col=c("red","blue","purple","pink","gray"))

print(cbind(beta0, theta0, a0, beta1, a1, sse1, beta2, theta2, a2, sse2, beta3, theta3, a3, sse3))

uscan.beta0 <- beta0; uscan.theta0 <- theta0; uscan.a0 <- a0; uscan.beta1 <- beta1; uscan.a1 <- a1; 
uscan.beta2 <- beta2; uscan.theta2 <- theta3; uscan.a2 <- a2;
uscan.beta3 <- beta3; uscan.theta3 <- theta3; uscan.a3 <- a3

# QQ Plots
# library(PtProcess)
uscan.lt.area.prob <- ppoints(uscan.lt.area.sort)
qp <- qpareto(uscan.lt.area.prob, beta1, a1, lower.tail=TRUE, log.p=FALSE)
head(qp); tail(qp)

qtap <- qtappareto(uscan.lt.area.prob, beta2, theta2, a2, lower.tail=TRUE, log.p=FALSE, tol=1e-8)
head(qtap); tail(qtap)

head(cbind(log10(uscan.lt.area.sort),qp,qtap,log10(qp),log10(qtap),uscan.lt.area.prob))

tail(cbind(log10(uscan.lt.area.sort),qp,qtap,log10(qp),log10(qtap),uscan.lt.area.prob))

plot(log10(uscan.lt.area.sort), log10(qtap), xlim=c(-1,6), ylim=c(-1,6), xlab="Observed log10(Area)", type="n", 
     ylab="log10(Theoretical Quantiles)", main=paste("QQPlot",uscan.lt.name,"Fires >",as.character(a0),"ha"))
abline(0.0, 1.0, col="gray", lwd=2)
points(log10(uscan.lt.area.sort),log10(qp), type="l", col="blue", lwd=2)
points(log10(uscan.lt.area.sort),log10(qtap), type="l", col="purple", lwd=2)
legend("topleft", legend=c("Pareto","Tapered Pareto"), lwd=2, col=c("blue","purple"))

# Compare Tapered-Pareto Fits
# compare observed and fitted
plot(log10(fpafod.lt.area.sort), log10(fpafod.lt.area.surv), xlab="log10(Area)", ylab="log10(S(x)", type="n", 
     xlim=c(-1,6), main=paste("Lightning Fires >",as.character(a0),"ha"))
abline(log10(fpafod.lt.area.min)/2, -0.5, col="gray", lwd=2)

points(log10(fpafod.lt.area.sort), log10(fpafod.lt.area.surv), col="red", pch=16, cex=0.7 )
fpafod.lt.fx.tappareto <- 1.0 - ((us.a2/fpafod.lt.area.sort)^us.beta2)*exp((us.a2-fpafod.lt.area.sort)/
                                                                             us.theta2)
points(log10(fpafod.lt.area.sort), log10(1.0 - fpafod.lt.fx.tappareto), col="pink", type="l", lwd=2) 

points(log10(cnfdb.lt.area.sort), log10(cnfdb.lt.area.surv), col="blue", pch=16, cex=0.7)
cnfdb.lt.fx.tappareto <- 1.0 - ((can.a2/cnfdb.lt.area.sort)^can.beta2)*exp((can.a2-cnfdb.lt.area.sort)/
                                                                             can.theta2)
points(log10(cnfdb.lt.area.sort), log10(1.0 - cnfdb.lt.fx.tappareto), col="lightblue", type="l", lwd=2) 

points(log10(uscan.lt.area.sort), log10(uscan.lt.area.surv), col="purple", pch=16, cex=0.7)
uscan.lt.fx.tappareto <- 1.0 - ((uscan.a2/uscan.lt.area.sort)^uscan.beta2)*exp((uscan.a2-uscan.lt.area.sort)/
                                                                                 uscan.theta2)
points(log10(uscan.lt.area.sort), log10(1.0 - uscan.lt.fx.tappareto), col="plum", type="l", lwd=2) 
legend("bottomleft", legend=c("US","Canada","US & Canada"), lwd=2, col=c("red","blue", "purple"))

##################################################################################################################
## Section: Book keeping - Clean memory close file connections
##################################################################################################################

# CLEAN MEMORY
rm(list = ls(all.names = TRUE))
raster::removeTmpFiles(h = 0)
flush.console()

##################################################################################################################
##                                                                                                              ##
##            Program end section to Investigate Fire Occurrence Database                                       ##
##                                                                                                              ##
##################################################################################################################
