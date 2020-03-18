##################################################################################################################
## Project: Process wind speed data from weather station and NARR data 
##
## Script purpose: Read RAWS data and extract wind speeds and directions
##
## Date: 1st January 2020
## Author: JR
##################################################################################################################

# Required testing many libraries. Be selective in initializing libraries
library(rlang) 
library(maptools) 
library(GISTools) 
library(dplyr)
library(plyr)
library(lwgeom)
library(sf)
library(sp)
library(gdalUtils)
library(stringr)
library(smoothr)
library(units)
library(geosphere)
library(fields)
library(rasterVis) 
library(OpenImageR)
library(gdalUtils)
library(SDMTools)
library(truncnorm)
library(ggplot2)
library(rgdal)
library(raster) 
library(rgeos)
library(lubridate)
library(data.table)
library(clifro)
library(openair)
library(circular)

# Projection string for meter projection US wide
AEAProj <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-110
+x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m'

# Sheffs preferred output projections
laea_proj4 <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

# Test for valid GDAL install
gdal_setInstallation()
valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
if(require(raster) && require(rgdal) && valid_install)
{
  print('TRUE')
}
#gdal_setInstallation(ignore.full_scan=FALSE)

# Set raster options
rasterOptions(maxmemory = 5.0e+09)

#if("package:adehabitatHR" %in% search()) detach("package:adehabitatHR", unload=TRUE)

memory.limit()
# set max memory usage 
memory.size(max=5.0e+09)

##################################################################################################################
## Section: Set global directory locations
##################################################################################################################

# Get the data directory from readtext
ROOT_DATA_DIR <- 'F:/ERCTimeSeriesAnalysis'
DATA_DIR <- ('Data') 
OUTPUT_DIR <- ('Output') 
DATA_YR <- 1819
DATA_YR_OUT <- 2019

# specify working directory; you will need to change this line if run elsewhere
setwd(ROOT_DATA_DIR)

# Check directory
getwd()
wd <- getwd()
print(paste0('Current working dir: ', wd))

# Temporary raster location
raster::rasterOptions(tmpdir = wd)

# Specify output folder suffix
Exprmnt <- 3

##################################################################################################################
## Section: Read Station data
##################################################################################################################

# Re-read concatenated data.table of all pyromes station data combined
cmbndstdatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedStationDataByPyrome03/AllPyromesOrigStatDat.rds')
stationdata.DT.all.pyrs <- readRDS(cmbndstdatfile) 
sapply(stationdata.DT.all.pyrs, class) 

# Reset data for all pyromes 
stationdata.DT <- copy(stationdata.DT.all.pyrs)
#stationdata.DT <- stationdata.DT[PYROME == '008']
sapply(stationdata.DT, class) 

# Highlight blank entries separate from 0 entries - assuming all 28 fields are character
indx <- which(sapply(stationdata.DT, is.character))
for (j in indx) set(stationdata.DT, i = grep("^[[:space:]]*$", stationdata.DT[[j]]), j = j, value = NA_character_) 

# Reformat date column
cols <- c('DATE')
stationdata.DT[, DATE := lapply(.SD, function(x){as.Date(x,format='%Y%m%d')}), .SDcols = cols]

# Create Correct time formatted column
stationdata.DT[, paste0('HOUR', 1:2) := tstrsplit(TIME, "(?<=.{2})", perl = TRUE)][]
timecols <- c('HOUR1', 'HOUR2')
stationdata.DT[, TIME := do.call(paste, c(.SD, sep = ":")), .SDcols = timecols]
stationdata.DT[, TIME := paste0(TIME,':00')]
stationdata.DT[, DTME := paste0(DATE,' ',TIME)]
stationdata.DT[, DATETIME := parse_date_time(DTME, 'Ymd HMS')]
stationdata.DT[, DTME := NULL]

# Extract data variables 
seltvar <- c('STATNUM','DATE','WNDDIR','10MWSPD','GSTWNDDIR','GSTWSPD','PYROME','DATETIME')
stationdata.DT <- stationdata.DT[, ..seltvar]

# Enforce to be numeric
convrttype <- c('WNDDIR','10MWSPD','GSTWNDDIR','GSTWSPD')
mtimes <- rep(c('numeric'),each=length(convrttype))
mtimes <- setNames(mtimes, convrttype)
set_colclass(stationdata.DT,mtimes)

# Remove complete duplicates
setkey(stationdata.DT, PYROME, STATNUM, DATETIME)
DUPFRST <- duplicated(stationdata.DT, by = key(stationdata.DT))
stationdata.DT[, DUPSTDEL := DUPFRST | c(tail(DUPFRST, 0), FALSE)]

# Test duplicates
TESTDUPFRST <- stationdata.DT[DUPSTDEL == TRUE]

# Remove dupliactes
stationdata.DT <- stationdata.DT[DUPSTDEL == FALSE]
stationdata.DT[,DUPSTDEL := NULL]

# Create 'daily' 1300 hour reference time stationdata.DT <- stationdata.DT[PYROME == 127]
stationdata.DT[, REFTIME := parse_date_time(paste0(DATE,' 13:00:00'), 'Ymd HMS')]
stationdata.DT[, DIFTIME := difftime(DATETIME, REFTIME, units = 'secs')]
stationdata.DT[, DIFTIME := as.numeric(str_replace_all(DIFTIME,'secs',''))]
stationdata.DT[, DAYM := day(DATETIME)]

# Set order 
setorder(stationdata.DT, PYROME, STATNUM, DATETIME)

# Create index of date occurrency
stationdata.DT[,INDX := 1:.N, by=c('PYROME','STATNUM','DATE')]

# Test multiples set dummy
stationdata.DT[, INDXf:='n'] 

# Create diff(DINDX) shifted upwards, padding last observation with 0 to avoid cycling
stationdata.DT[, DINDX := c(diff(INDX, lag=1), 0), by=c('PYROME','STATNUM','DATE')]
stationdata.DT[ DINDX > 0, INDXf:='l']

idx = stationdata.DT[, {
  ix = tail(.I[DINDX==1L], 1);
  iy = (ix+1L)*((ix+1L) <= .I[.N] | NA) 
  list(idx = c(ix, iy))
}, by = list(PYROME,STATNUM,DATE)]$idx

DUPSSCND <- stationdata.DT[idx][, INDXf := c("l")]
if (nrow(DUPSSCND) != 0) {
  print("Multiple duplicate entries by day removed")
  stationdata.DT <- stationdata.DT[!idx][, INDXf := c("l")]
} else { 
  print("No duplicate entries found")  
}

# Clean output
stationdata.DT[ ,`:=`(REFTIME = NULL, DIFTIME = NULL, DAYM = NULL, INDX = NULL, 
                      INDXf = NULL, DATETIME = NULL, DINDX = NULL)] 
stationdata.DT[1:10,]

# Add month column
stationdata.DT[,MONTH := month(DATE)]
stationdata.DT[,MONAB := as.character(month(ymd(010101) + months(MONTH-1),label=TRUE,abbr=TRUE))]
setnames(stationdata.DT,'10MWSPD','WS')
setnames(stationdata.DT,'WNDDIR','WD')

# 0 wind direction implies unknown
stationdata.DT[, WD := WD ][WD == 0, WD := NA]
stationdata.DT <- stationdata.DT[!is.na(WD)]
stationdata.DT <- stationdata.DT[!is.na(WS)]

# Write concatenated data.table of all pyromes station data combined
cmbndstdatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedStationDataByPyrome03/AllPyromesOrigWindStatDat_DDUPED.rds')
saveRDS(stationdata.DT, cmbndstdatfile)

##################################################################################################################
## Section: Summarize station level wind data
##################################################################################################################

# Example Fire family plus input: (Check units: 1 miles/hour = 0.45 meters/sec)
#
#Wind Speed vs Dir
#6 8
#September:
#speed 	  45 	  90 	  135 	180 	225 	270 	315 	360
#5.00  	  7.07 	10.25	10.20	8.98 	13.10	15.31	7.89 	5.35
#10.00  	0.36 	0.76 	1.58 	3.70 	10.69	3.32 	0.65 	0.47
#15.00  	0.00 	0.00 	0.01 	0.10 	0.16 	0.03 	0.00 	0.00
#20.00  	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00
#25.00  	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00
#30.00  	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00
#0.00

# Read concatenated data.table of all pyromes station data combined
cmbndstdatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedStationDataByPyrome03/AllPyromesOrigWindStatDat_DDUPED.rds')
stationdata.DT <- readRDS(cmbndstdatfile) 
sapply(stationdata.DT, class) 

# Test pyrome 8 station 100423
#stationdata.DT <- stationdata.DT[PYROME == '008']
#stationdata.DT <- stationdata.DT[STATNUM == '100423']

# Clean wind speeds
stationdata.DT[,WS := WS][WS >= 30, WS := 27.5]

# Enforce to be numeric
convrttype <- c('WD','WS')
mtimes <- rep(c('numeric'),each=length(convrttype))
mtimes <- setNames(mtimes, convrttype)
set_colclass(stationdata.DT,mtimes)
sapply(stationdata.DT, class) 

# Bin directions into 45 degree bins
dirres <- 45
dir.breaks <- c(-dirres/2,
                seq(dirres/2, 360-dirres/2, by = dirres),
                360+dirres/2)  
dir.labels <- c(paste(360-dirres/2,"-",dirres/2),
                paste(seq(dirres/2, 360-3*dirres/2, by = dirres),
                      "-",
                      seq(3*dirres/2, 360-dirres/2, by = dirres)),
                paste(360-dirres/2,"-",dirres/2))

stationdata.DT[, WDBIN := factor(cut(WD,labels = FALSE, breaks = dir.breaks,
                                     include.lowest = FALSE, dig.lab = 3,
                                     LEFT = TRUE))]

# Enforce WDBIN to be numeric
convrttype <- c('WDBIN')
mtimes <- rep(c('numeric'),each=length(convrttype))
mtimes <- setNames(mtimes, convrttype)
set_colclass(stationdata.DT,mtimes)
sapply(stationdata.DT, class) 
#stationdata.DT[, WDBIN := WDBIN][WDBIN == '(-22.5,22.5]', WDBIN := '(338,382]']
stationdata.DT[, WDBIN := WDBIN][WDBIN == 1, WDBIN := 9]
stationdata.DT[, WDBIN := (WDBIN - 1)*45]

# Bin windspeeds into 5 miles per hour bins
brks <- seq(0, 30, length=7)
stationdata.DT[,WSBIN := findInterval(WS, brks)]

# Enforce WSBIN to be numeric
convrttype <- c('WSBIN')
mtimes <- rep(c('numeric'),each=length(convrttype))
mtimes <- setNames(mtimes, convrttype)
set_colclass(stationdata.DT,mtimes)
sapply(stationdata.DT, class) 
stationdata.DT[, WSBIN := WSBIN][WSBIN == 0, WSBIN := 1]
stationdata.DT[, WSBIN := WSBIN * 5]

# Finally aggregate 
Wind.stats <- stationdata.DT[, list(.N), by=c('PYROME','MONTH','WDBIN','WSBIN')]
setorder(Wind.stats, PYROME, MONTH, WDBIN, WSBIN)
Wind.stats.cst <- dcast(Wind.stats, FUN = length, PYROME + MONTH + WSBIN ~ WDBIN)
Wind.stats.cst[is.na(Wind.stats.cst)] <- 0
sapply(Wind.stats.cst, class) 
setorder(Wind.stats.cst, PYROME, MONTH, WSBIN)

# Create proportional table by pyrome
#collength <- dim(Wind.stats.cst)[2]
#selected_cols <- colnames(Wind.stats.cst)[4:collength]
#Wind.stats.props <- Wind.stats.cst[, prop.table(.SD), by = c('PYROME','MONTH'), .SDcols = selected_cols]

# Test results
#test1 <- Wind.stats.cst[, prop.table(.SD), by = c('PYROME','MONTH'), .SDcols = selected_cols]
#sum(subset(test1, MONTH == 1, select = selected_cols))

# Cross Join to expand data
WSbrks <- seq(5, 30, length=6)
Mthsbrks <- seq(1, 12, length=12)
Pyrbrks <- seq(1, 128, length=128)
Pyrbrks <- sprintf('%03d',Pyrbrks)
bins <- CJ(Pyrbrks, Mthsbrks, WSbrks)
setDT(bins)
old <- names(bins)
new <- c('PYROME', 'MONTH', 'WSBIN')
setnames(bins, old, new )
sapply(bins, class) 

# Expand original table
setkey(Wind.stats.cst, PYROME, MONTH, WSBIN)
setkey(bins, PYROME, MONTH, WSBIN)
stationdata.DT.exp <- Wind.stats.cst[bins, allow.cartesian=TRUE] 
stationdata.DT.exp[is.na(stationdata.DT.exp)] <- 0

# Renormalize data by pyrome and month
colsToSum <- colnames(stationdata.DT.exp)[4:11] 
stationdata.DT.exp[, sumcheck :=rowSums(.SD, na.rm = TRUE), .SDcols = colsToSum]
stationdata.DT.exp[ , sumcheckr :=lapply(.SD, sum, na.rm=TRUE), by=c('PYROME','MONTH'), .SDcols='sumcheck'] 
stationdata.DT.exp[, (colsToSum) := .SD * 1.0/sumcheckr, .SDcols = colsToSum]
stationdata.DT.exp[, sumcheck := NULL]
stationdata.DT.exp[, sumcheckr := NULL]
stationdata.DT.exp[, c(colsToSum) := lapply(.SD, function(x) trimws(format(round(x, 6))) ), .SDcols=colsToSum]

# Output each Pyrome for FireFamilyPlus import
DTlist <- split(stationdata.DT.exp, by = 'PYROME')
for (j in 1:length(DTlist)) {
  DTDummy <- as.data.table(DTlist[[j]])
  fileindx <- sprintf("%03d", j)
  filedir <- file.path(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedStationDataByPyrome04'))
  fileout <- paste0('OrigWindStatDat',fileindx,'.rds')
  if (!dir.exists(filedir)){
    dir.create(filedir,recursive = TRUE)
  } else {
    print("Dir already exists!")
  }
  fileout <- paste0(filedir,'/',fileout)
  
  # Write concatenated data.table
  saveRDS(DTDummy,fileout)
}


p.wr2 <- plot.windrose(data = stationdata.DT,
                       spd = 'WS',
                       dir = 'WD',
                       spdseq = c(0,5,10,15,20,25,30))

p.wr3 <- p.wr2 + facet_wrap(~MONTH ,
                            ncol = 4)
# and remove labels for clarity
p.wr3 <- p.wr3 + theme(axis.text.x = element_blank(),
                       axis.title.x = element_blank())


library(openair)
windRose(stationdata.DT,
         wd = 'WD',
         ws = 'WS', angle = 22.5, width = 0.2, grid.line = 5, paddle = F)

# one windRose for each year
windRose(mydata,type = "month")

##################################################################################################################
## Section: Read NARR Gridded data
##################################################################################################################

# Re-read combined input data
cmbndgddatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedNARRGridDataByPyrome03/AllPyromesOrigNARRGridDat.rds')
griddeddata.DT.all.pyrs <- readRDS(cmbndgddatfile) 
sapply(griddeddata.DT.all.pyrs, class) 

# Reset data for all pyromes 
griddeddata.DT <- copy(griddeddata.DT.all.pyrs)
#griddeddata.DT <- griddeddata.DT[PYROME == '001']

# Update year, rain duration and amount to numeric fields
set_colclass(griddeddata.DT, c(YEAR='numeric', RnDr='numeric', Rain='numeric'))
class(griddeddata.DT$RnDr)
sapply(griddeddata.DT, class)

# Highlight blank entries separate from 0 entries - assuming all 25 fields are character
indx <- which(sapply(griddeddata.DT, is.character))
for (j in indx) set(griddeddata.DT, i = grep("^[[:space:]]*$", griddeddata.DT[[j]]), j = j, value = NA_character_) 

# Reformat and creat date column
griddeddata.DT[, D4YR := 0]
griddeddata.DT[, D4YR := D4YR][(YEAR <= 12), D4YR := D4YR + 2000][]
griddeddata.DT[, D4YR := D4YR][(D4YR < 1900), D4YR := D4YR + 1900][]
griddeddata.DT[, D4YR := D4YR + YEAR][]
set_colclass(griddeddata.DT, c(D4YR='character'))
cols <- c('D4YR', 'MONTH', 'DAY')
griddeddata.DT[, DATE := do.call(paste, c(.SD, sep = "")), .SDcols = cols]
cols <- c('DATE')
griddeddata.DT[, DATE := lapply(.SD, function(x){as.Date(x,format='%Y%m%d')}), .SDcols = cols]
griddeddata.DT[, DTME := paste0(DATE,' 13:00:00')]
griddeddata.DT[, DATETIME := parse_date_time(DTME, 'Ymd HMS')]
griddeddata.DT[, DTME := NULL]

# Extract data variables 
seltvar <- c('SIGName','DATE','WDIR','GSpd','PYROME','DATETIME')
griddeddata.DT <- griddeddata.DT[, ..seltvar]

# Enforce to be numeric
convrttype <- c('WDIR','GSpd')
mtimes <- rep(c('numeric'),each=length(convrttype))
mtimes <- setNames(mtimes, convrttype)
set_colclass(griddeddata.DT,mtimes)

# Rename to be consistent with station data
oldnames <- c('SIGName', 'WDIR', 'GSpd')
newnames <- c('STATNUM','WD','WS')
setnames(griddeddata.DT, old = oldnames, new = newnames)

# Remove complete duplicates
setkey(griddeddata.DT, PYROME, STATNUM, DATETIME)
DUPFRST <- duplicated(griddeddata.DT, by = key(griddeddata.DT))
griddeddata.DT[, DUPSTDEL := DUPFRST | c(tail(DUPFRST, 0), FALSE)]

# Test duplicates
TESTDUPFRST <- griddeddata.DT[DUPSTDEL == TRUE]

# Remove duplicates
griddeddata.DT <- griddeddata.DT[DUPSTDEL == FALSE]
griddeddata.DT[,DUPSTDEL := NULL]

# Create 'daily' 1300 hour reference time stationdata.DT <- stationdata.DT[PYROME == 127]
griddeddata.DT[, REFTIME := parse_date_time(paste0(DATE,' 13:00:00'), 'Ymd HMS')]
griddeddata.DT[, DIFTIME := difftime(DATETIME, REFTIME, units = 'secs')]
griddeddata.DT[, DIFTIME := as.numeric(str_replace_all(DIFTIME,'secs',''))]
griddeddata.DT[, DAYM := day(DATETIME)]

# Set order 
setorder(griddeddata.DT, PYROME, STATNUM, DATETIME)

# Create index of date occurrency
griddeddata.DT[,INDX := 1:.N, by=c('PYROME','STATNUM','DATE')]

# Test multiples set dummy
griddeddata.DT[, INDXf:='n'] 

# Create diff(DINDX) shifted upwards, padding last observation with 0 to avoid cycling
griddeddata.DT[, DINDX := c(diff(INDX, lag=1), 0), by=c('PYROME','STATNUM','DATE')]
griddeddata.DT[ DINDX > 0, INDXf:='l']

idx = griddeddata.DT[, {
  ix = tail(.I[DINDX==1L], 1);
  iy = (ix+1L)*((ix+1L) <= .I[.N] | NA) 
  list(idx = c(ix, iy))
}, by = list(PYROME,STATNUM,DATE)]$idx

DUPSSCND <- griddeddata.DT[idx][, INDXf := c("l")]
if (nrow(DUPSSCND) != 0) {
  print("Multiple duplicate entries by day removed")
  griddeddata.DT <- griddeddata.DT[!idx][, INDXf := c("l")]
} else { 
  print("No duplicate entries found")  
}

# Clean output
griddeddata.DT[ ,`:=`(REFTIME = NULL, DIFTIME = NULL, DAYM = NULL, INDX = NULL, 
                      INDXf = NULL, DATETIME = NULL, DINDX = NULL)] 
griddeddata.DT[1:10,]

# Add month column
griddeddata.DT[,MONTH := month(DATE)]
griddeddata.DT[,MONAB := as.character(month(ymd(010101) + months(MONTH-1),label=TRUE,abbr=TRUE))]

# 0 wind direction implies unknown
griddeddata.DT[, WD := WD ][WD == 0, WD := NA]
griddeddata.DT <- griddeddata.DT[!is.na(WD)]
griddeddata.DT <- griddeddata.DT[!is.na(WS)]

# Write concatenated data.table of all pyromes station data combined duplicates removed
cmbndgddatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedNARRGridDataByPyrome03/AllPyromesOrigWindNARRGridDat_DDUPED.rds')
saveRDS(griddeddata.DT, cmbndgddatfile) 

##################################################################################################################
## Section: Summarize station level wind data
##################################################################################################################

# Example Fire family plus input: (Check units: 1 miles/hour = 0.45 meters/sec)
#
#Wind Speed vs Dir
#6 8
#September:
#speed 	  45 	  90 	  135 	180 	225 	270 	315 	360
#5.00  	  7.07 	10.25	10.20	8.98 	13.10	15.31	7.89 	5.35
#10.00  	0.36 	0.76 	1.58 	3.70 	10.69	3.32 	0.65 	0.47
#15.00  	0.00 	0.00 	0.01 	0.10 	0.16 	0.03 	0.00 	0.00
#20.00  	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00
#25.00  	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00
#30.00  	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00
#0.00

# Read concatenated data.table of all pyromes station data combined
cmbndgddatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedNARRGridDataByPyrome03/AllPyromesOrigWindNARRGridDat_DDUPED.rds')
griddeddata.DT <- readRDS(cmbndgddatfile) 
sapply(griddeddata.DT, class) 

# Test pyrome 8 station 100423
#griddeddata.DT <- griddeddata.DT[PYROME == '008']
#griddeddata.DT <- griddeddata.DT[STATNUM == '100423']

# Clean wind speeds
griddeddata.DT[,WS := WS][WS >= 30, WS := 27.5]

# Enforce to be numeric
convrttype <- c('WD','WS')
mtimes <- rep(c('numeric'),each=length(convrttype))
mtimes <- setNames(mtimes, convrttype)
set_colclass(griddeddata.DT,mtimes)
sapply(griddeddata.DT, class) 

# Bin directions already in data 1 thru 8
griddeddata.DT[,WDBIN := WD]
griddeddata.DT[,WDBIN := WDBIN * 45]

# Bin windspeeds into 5 miles per hour bins
brks <- seq(0, 30, length=7)
griddeddata.DT[,WSBIN := findInterval(WS, brks)]

# Enforce WSBIN to be numeric
convrttype <- c('WSBIN')
mtimes <- rep(c('numeric'),each=length(convrttype))
mtimes <- setNames(mtimes, convrttype)
set_colclass(griddeddata.DT,mtimes)
sapply(griddeddata.DT, class) 
griddeddata.DT[, WSBIN := WSBIN][WSBIN == 0, WSBIN := 1]
griddeddata.DT[, WSBIN := WSBIN * 5]

# Finally aggregate 
Wind.stats <- griddeddata.DT[, list(.N), by=c('PYROME','MONTH','WDBIN','WSBIN')]
setorder(Wind.stats, PYROME, MONTH, WDBIN, WSBIN)
Wind.stats.cst <- dcast(Wind.stats, FUN = length, PYROME + MONTH + WSBIN ~ WDBIN)
Wind.stats.cst[is.na(Wind.stats.cst)] <- 0
sapply(Wind.stats.cst, class) 
setorder(Wind.stats.cst, PYROME, MONTH, WSBIN)

# Create proportional table by pyrome
#collength <- dim(Wind.stats.cst)[2]
#selected_cols <- colnames(Wind.stats.cst)[4:collength]
#Wind.stats.props <- Wind.stats.cst[, prop.table(.SD), by = c('PYROME','MONTH'), .SDcols = selected_cols]

# Test results
#test1 <- Wind.stats.cst[, prop.table(.SD), by = c('PYROME','MONTH'), .SDcols = selected_cols]
#sum(subset(test1, MONTH == 1, select = selected_cols))

# Cross Join to expand data
WSbrks <- seq(5, 30, length=6)
Mthsbrks <- seq(1, 12, length=12)
Pyrbrks <- seq(1, 128, length=128)
Pyrbrks <- sprintf('%03d',Pyrbrks)
bins <- CJ(Pyrbrks, Mthsbrks, WSbrks)
setDT(bins)
old <- names(bins)
new <- c('PYROME', 'MONTH', 'WSBIN')
setnames(bins, old, new )
sapply(bins, class) 

# Expand original table
setkey(Wind.stats.cst, PYROME, MONTH, WSBIN)
setkey(bins, PYROME, MONTH, WSBIN)
griddeddata.DT.exp <- Wind.stats.cst[bins, allow.cartesian=TRUE] 
griddeddata.DT.exp[is.na(griddeddata.DT.exp)] <- 0

# Renormalize data by pyrome and month
colsToSum <- colnames(griddeddata.DT.exp)[4:11] 
griddeddata.DT.exp[, sumcheck :=rowSums(.SD, na.rm = TRUE), .SDcols = colsToSum]
griddeddata.DT.exp[ , sumcheckr :=lapply(.SD, sum, na.rm=TRUE), by=c('PYROME','MONTH'), .SDcols='sumcheck'] 
griddeddata.DT.exp[, (colsToSum) := .SD * 1.0/sumcheckr, .SDcols = colsToSum]
griddeddata.DT.exp[, sumcheck := NULL]
griddeddata.DT.exp[, sumcheckr := NULL]
griddeddata.DT.exp[, c(colsToSum) := lapply(.SD, function(x) trimws(format(round(x, 6))) ), .SDcols=colsToSum]

# Output each Pyrome for FireFamilyPlus import
DTlist <- split(griddeddata.DT.exp, by = 'PYROME')
for (j in 1:length(DTlist)) {
  DTDummy <- as.data.table(DTlist[[j]])
  fileindx <- sprintf("%03d", j)
  filedir <- file.path(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedNARRGridDataByPyrome03'))
  fileout <- paste0('OrigWindNARRGridDat',fileindx,'.rds')
  if (!dir.exists(filedir)){
    dir.create(filedir,recursive = TRUE)
  } else {
    print("Dir already exists!")
  }
  fileout <- paste0(filedir,'/',fileout)
  
  # Write concatenated data.table
  saveRDS(DTDummy,fileout)
}

##################################################################################################################
## Section: Define miscellaneous functions
##################################################################################################################

# Remove NAs from data.table
remve_na = function(DT) {
  # either of the following for loops
  
  # by name :
  #for (j in names(DT))
    #set(DT,which(is.na(DT[[j]])),j,0)
  
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,0)
}

# Join data to DT of consecutive days covering tim eperiod of data
Merged.DTs <- Reduce(function(x,y) merge(x,y,all=TRUE) , filesappend)

# Remove entries in station data where all variables are NA  
drop_rows_all_na <- function(x, pct=1) x[!rowSums(is.na(x)) >= ncol(x)*pct,]

#Returns all items in a list that are not contained in toMatch
#toMatch can be a single item or a list of items
exclude <- function (theList, toMatch){
  return(setdiff(theList,include(theList,toMatch)))
}

#Returns all items in a list that ARE contained in toMatch
#toMatch can be a single item or a list of items
include <- function (theList, toMatch){
  matches <- unique (grep(paste(toMatch,collapse="|"), 
                          theList, value=TRUE))
  return(matches)
}

# Convert data.table classes 
#dt <- data.table(i=1:3,f=3:1)
#set_colclass(dt, c(i="numeric",f="numeric" ))
#class(dt$i)
#sapply(dt, class)
".." <- function (x) 
{
  stopifnot(inherits(x, "character"))
  stopifnot(length(x) == 1)
  get(x, parent.frame(4))
}

set_colclass <- function(x, class){
  stopifnot(all(class %in% c("integer", "numeric", "double","factor","character")))
  for(i in intersect(names(class), names(x))){
    f <- get(paste0("as.", class[i]))
    x[, (..("i")):=..("f")(get(..("i")))]
  }
  invisible(x)
}

# Test if year is leap or not
is.leapyear <- function(year){
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
}

# Test sequence of numbers for consecutive differences of 1, can be changed
seqle <- function(x,incr=1) { 
  if(!is.numeric(x)) x <- as.numeric(x) 
  n <- length(x)  
  #y <- x[-1L] != x[-n] + incr 
  y <- abs(x[-1L] - x[-n] - incr) > .Machine$double.eps ^ 0.5
  i <- c(which(y|is.na(y)),n) 
  list(lengths = diff(c(0L,i)),
       values = x[head(c(0L,i)+1L,-1L)]) 
} 

##################################################################################################################
## Section: Book keeping - Clean memory close file connections
##################################################################################################################

# CLEAN MEMORY
rm(list = ls(all.names = TRUE))
raster::removeTmpFiles(h = 0)
flush.console()

##################################################################################################################
##                                                                                                              ##
##            Program end section to Process wind speed data from weather station and NARR data                 ##
##                                                                                                              ##
##################################################################################################################






