##################################################################################################################
## Project: Extracting precipitation data from RAWS data 
##
## Script purpose: Read .fw13 formatted RAWS data and extract precipitation amounts and durations
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
ROOT_DATA_DIR <- 'D:/ERCTimeSeriesAnalysis'
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
## Section: Read both Station and NARR files defining closest of each to every Pyrome centroid
##################################################################################################################

# Read grid point list of all pyromes 
gridfile <- paste0('D:/ERCTimeSeriesAnalysis/Data/Pyromes_2019/','Top4CloseGridPointsAllPyromesCWAtrbs.rds')
NARRgrd.DT <- readRDS(gridfile)
#setorder(NARRgrd.DT, grid_x, grid_y)
#NARRgrd.DT[PYROME == '010']

# Read station list of all pyromes 
stationfile <- paste0('D:/ERCTimeSeriesAnalysis/Data/Pyromes_2019/','Top5CloseStationsAllPyromes_2019.rds')
STATgrd.DT <- readRDS(stationfile)
STATgrd.DT[PYROME == '999']
STATgrd.DT[, PYROME := NULL]

# Fix issue of station not falling within pyrome, coded as 999
newpyrome <- rep(1:128, each=5)
STATgrd.DT <- cbind(STATgrd.DT,newpyrome)
setnames(STATgrd.DT, 'newpyrome', 'PYROME')
STATgrd.DT[, PYROME := sprintf("%03d", as.integer(PYROME))] 

##################################################################################################################
## Section: Edit STATION location RAWs data .fw13 formatted
##################################################################################################################

# File format https://fam.nwcg.gov/fam-web/weatherfirecd/13.htm
# 
# Field | Field Name                                                                        | Columns 
# 01    Record type (W13)                                                                     01-03 
# 02    Station number                                                                        04-09
# 03    Observation date (YYYYMMDD)                                                           10-17
# 04    Oservation time (0000-2359)                                                           18-21
# 05    Observation type (O=Published NFDRS, R=RAWS, F=Forecast, X=Other). If type is R and 
#       observation time is the station’s standard NFDRS observation time(RS) and the State of 
#       weather code is not blank, State of weather code and Wet Flag for this observation were 
#       estimated by the WIMS RAWS Gateway routines. (Note: NFDR processing programs should be 
#                                              modified to compute NFDRS indexes for the 
#                                              traditional “O” and the “R” at RS time if 
#                                              the SOW is not null.)                          22-22
# 06    State of the weather code                                                             23-23
# 07    Dry bulb temperature ( F or C dependant on Column 63)                                 24-26
# 08    Atmospheric moisture (Wet bulb, RH or dew point dependant on column 62)               27-29 
# 09    Wind direction (Azimuth measured from true north. 0 means no wind 
#       direction, 360 is north)                                                              30-32
# 10    Average windspeed over a 10-minute period (miles or kilometers per hour 
#       based on Measurement Type code)                                                       33-35
# 11    Measured 10-hour time lag fuel moisture                                               36-37 
# 12    Maximum Temperature (degrees Fahrenheit or degrees Celsius based on 
#       Measurement Type code [col. 63])                                                      38-40
# 13    Minimum Temperature (degrees Fahrenheit or degrees Celsius based on 
#       Measurement Type code [col. 63])                                                      41-43 
# 14    Maximum relative humidity (percent)                                                   44-46 
# 15    Minimum relative humidity (percent)                                                   47-49
# 16    Precipitation duration (hours)                                                        50-51
# 17    Precipitation amount based on Measurement Type code [col. 63]. 
#       Blanks=no precipitation. U.S. measurement: inches with implied decimal 
#       nn.nnn format; trace shown as 00005. Metric measurement: 
#       measured in millimeters, no implied decimal; trace shown as 00001                     52-56
# 18    Wet flag (Y/N)                                                                        57-57
# 19    Herbaceous greenness factor (0-20)                                                    58-59 
# 20    Shrub greenness factor (0-20)                                                         60-61
# 21    Moisture Type code (1=Wet bulb, 2=Relative Humidity, 3=Dewpoint)                      62-62
# 22    Measurement Type code: 1=U.S.,2=Metric. Affects temperature (Fahrenheit or 
#       Celsius), wind (miles or kilometers per hour), and precipitation (decimal inches or 
#       millimeters - 1 = US (Precipitation amount – 24 hour), 2 = Metric (Precipitation 
#       amount – 24 hour), 3 = US (Precipitation amount – hourly), 4 = Metric 
#       (Preceiptation amount – hourly)                                                       63-63
# 23    Season code (1=Winter, 2=Spring, 3=Summer, 4=Fall)                                    64-64
# 24    Solar radiation (watts per square meter)                                              65-68
# 25    Wind direction of peak gust during the hour, degrees. Zero means 
#       no wind direction, 360 is north                                                       69-71
# 26    Speed of peak gust during the hour. (miles or kilometers per hour based on 
#       Measurement Type code)                                                                72-74
# 27    Snow Flag (Y/N). Signals fuels over the fire danger rating area are snow covered      75-75

# Pass data from text DT to final DT
Headernames <- c('RTYPE','STATNUM','DATE','TIME','OBTYPE','SOW','TEMP','ATMOIST','WNDDIR','10MWSPD',
                 '10HRFUEL','MaxT','MinT','MxRH','MnRH','RnDr','RAIN','WFLG','HGRSFACT','SHGRSFACT',
                 'MOISTTYP','MEASTYP', 'SEACODE','SRAD','GSTWNDDIR','GSTWSPD','SNWFLG')
                 
colconvrt <- c('TEMP','ATMOIST','MaxT','MinT','MxRH','MnRH')

# Column data entered by hand from description above
end_col <- 75 #cumsum(cols)
start_col <- c(1,4,10,18,22,23,24,27,30,33,36,38,41,44,47,50,52,57,58,60,62,63,64,65,69,72,75) #end_col - start_col + 1
end_col <- c(3,9,17,21,22,23,26,29,32,35,37,40,43,46,49,51,56,57,59,61,62,63,64,68,71,74,75)
start_end <- cbind(start_col, end_col) # matrix of start and end positions
column_widths <- end_col - start_col + 1
COLWIDS <- as.data.frame(rbind(column_widths))
COLWIDS <- setDT(COLWIDS)
setnames(COLWIDS, old = 1:ncol(COLWIDS), new = Headernames)
COLWIDS <- as.list(COLWIDS)

# Read station point list for all pyromes 
stationfile <- paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/','Pyromes_2019/Top5CloseStationsAllPyromes_2019.rds')
stationlistfile.dt.rnkd <- readRDS(stationfile) 
sapply(stationlistfile.dt.rnkd, class) 

# Select subset of pyromes
#pyrssub <- c('001','005', '008', '009', '010', '011', '012', '013', '014', 
#             '015', '020', '026', '075', '101', '118', '122', '127')

pyrssub <- sprintf('%03d',seq(1:128))

# Specify output location
outloc <- paste0('Expmnt_',sprintf("%03d", as.integer(Exprmnt)),'/')

# Initialize output file of combined station data for each pyrome
stationdata.DT <- data.table(NULL)

# Combined station data for all pyromes
stationdata.DT.all.pyrs <- data.table(NULL)

# Select pyrome of interest to speed processing
#stationlistfile.dt.rnkd <- stationlistfile.dt.rnkd[PYROME == '001']

# Loop through stationlistfile.dt.rnkd to process each weather station file by Pyrome
for (iLine in 1:nrow(stationlistfile.dt.rnkd)) { 
  # Evaluate pyrome code
  PYROMEIN <- unlist(stationlistfile.dt.rnkd[iLine,PYROME])
  idx_one_in_two <- match(PYROMEIN, unlist(pyrssub)) 
  # Select subset of pyromes
  if (!is.na(idx_one_in_two)) { 
    # Set working directory
    PYROMEIN <- pyrssub[idx_one_in_two]
    initial_path <- file.path(paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/','Pyromes_2019/StationWeatherData/',PYROMEIN))
    setwd(initial_path)
    wd <- getwd()
    print(paste0('Current working dir: ', wd))   
    
    # File to process for climate conditioning
    filenamein <- stationlistfile.dt.rnkd[iLine,FLNAME]
    distrankin <- stationlistfile.dt.rnkd[iLine,distrank]
    stationid <- stationlistfile.dt.rnkd[iLine, STATID] 
    fileindx <- sprintf("%03d", distrankin)
    filename <- paste0(fileindx,'_wx',stationid,'.fw13')
    filein <- paste0(initial_path,'/',filename)     
 
    # Read data to data.table line by line
    grdfile.DT <- fread(filein, colClasses = "character",
                        sep = "\n", header = FALSE, verbose = TRUE)
    
    # Function to pass data by column widths
    textin <- lapply(grdfile.DT, function(x) {
      apply(start_end, 1, function(y) substr(x, y[1], y[2])) 
    })
    grdfile.DT <- data.table(textin$V1)
    setnames(grdfile.DT, old = 1:ncol(grdfile.DT), new = Headernames)
    
    # Drop blank rows remove whitespace
    grdfile.DT <- subset(grdfile.DT[RTYPE  != "",])
    grdfile.DT[, PYROME := PYROMEIN]

    # Append pyrome data to master data DT
    if (distrankin == 1) {
      # Initialize DT for current pyrome
      stationdata.DT <- grdfile.DT

      # Append all pyromes 
      filesappend <- list(stationdata.DT.all.pyrs, grdfile.DT) 
      stationdata.DT.all.pyrs <- rbindlist(filesappend, use.names=TRUE, fill=TRUE)       
    } else {
      # Write 5 station locations
      filesappend <- list(stationdata.DT, grdfile.DT) 
      stationdata.DT <- rbindlist(filesappend, use.names=TRUE, fill=TRUE) 
      
      # Append all pyromes 
      filesappend <- list(stationdata.DT.all.pyrs, grdfile.DT) 
      stationdata.DT.all.pyrs <- rbindlist(filesappend, use.names=TRUE, fill=TRUE) 
    }
    rm(grdfile.DT)
  }
}

# Write concatenated data.table of all pyromes station data combined
cmbndstdatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedStationDataByPyrome/AllPyromesOrigStatDat.rds')
saveRDS(stationdata.DT.all.pyrs, cmbndstdatfile)

# Re-read combined input data
cmbndstdatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedStationDataByPyrome/AllPyromesOrigStatDat.rds')
stationdata.DT.all.pyrs <- readRDS(cmbndstdatfile) 
sapply(stationdata.DT.all.pyrs, class) 

# Reset data for all pyromes 
stationdata.DT <- copy(stationdata.DT.all.pyrs)
#stationdata.DT <- stationdata.DT[PYROME == '001']

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

# Update rain duration and amount to numeric fields
#set_colclass(stationdata.DT, c(RnDr='numeric', RAIN='numeric'))
#class(stationdata.DT$RnDr)
#sapply(stationdata.DT, class)
#setorder(stationdata.DT, -RAIN)

# Extract data variables 
seltvar <- c('STATNUM','DATE','RnDr','RAIN','TEMP','ATMOIST','MaxT','MinT','MxRH','MnRH','MEASTYP','PYROME','DATETIME')
stationdata.DT <- stationdata.DT[, ..seltvar]

# Enforce to be numeric
convrttype <- c('RnDr','RAIN','TEMP','ATMOIST','MaxT','MinT','MxRH','MnRH','MEASTYP')
mtimes <- rep(c('numeric'),each=length(convrttype))
mtimes <- setNames(mtimes, convrttype)
set_colclass(stationdata.DT,mtimes)

# Convert decimal inches 00000, format nn.nnn to number
for (j in names(stationdata.DT)) set(stationdata.DT,which(is.na(stationdata.DT[[j]])),j,0)
for(i in 4L){ set(stationdata.DT, i=NULL, j=i, value=0.001*stationdata.DT[['RAIN']]) }

# Convert pyrome station data to data frame for plotting
sapply(stationdata.DT, class)
stationdata.DF <- copy(stationdata.DT)
stationdata.DF <- setDF(stationdata.DF)

ggplot() +
  # plot all simulations
  geom_line(data = stationdata.DF, aes(x = DATE, y = RAIN, color = STATNUM, group=STATNUM), alpha = 1, size = 0.65) +
  labs(x="Date",y=expression('Precipitation Amount (inches)')) +
  theme(text = element_text(size=10))

ggplot(stationdata.DF, aes(x = DATE, y = RAIN)) + 
  geom_line() + 
  facet_wrap(~ STATNUM, scales = 'free_y', ncol = 1)

##################################################################################################################
## Section: Clean stationdata.DT where days have multiple enetries but different times
##          Require closest data to 1300 hours, to represent day, but averages evaluated over day
##################################################################################################################

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

# Write concatenated data.table of all pyromes station data combined
cmbndstdatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedStationDataByPyrome/AllPyromesOrigStatDat_DDUPED.rds')
saveRDS(stationdata.DT, cmbndstdatfile)

##################################################################################################################
## Section: Create combined time series of precipitation for each pyrome from station data
##################################################################################################################

# Total number of pyromes
numpyr <- 128
pyromes <- seq(1,numpyr,1)

# Loop through pyromes 
stationdata.DT[, PYROME := sprintf("%03d", as.integer(PYROME))] 
for (ipyr in pyromes) { 
  # Evaluate pyrome code
  PYROMEIN <- sprintf("%03d", as.integer(pyromes[ipyr]))
  idx_one_in_two <- match(PYROMEIN, unlist(pyrssub))
    # Select subset of pyromes
  if (!is.na(idx_one_in_two)) { 
    # Get data for this pyrome
    stationdata.DT.pyr <- stationdata.DT[PYROME == PYROMEIN]
    setcolorder(stationdata.DT.pyr, 'DATE')
    cols <- names(stationdata.DT.pyr)
    numcols <- length(cols)

    # Get min and max date values for complete timeseries
    mindate <- min(stationdata.DT.pyr[,DATE])
    maxdate <- max(stationdata.DT.pyr[,DATE])  
    cont_date_time <- as.data.table(data.frame(seq(as.Date(mindate), as.Date(maxdate), 'days')))
    setnames(cont_date_time, 1, 'DATE')

    # Split data.atble by station merge with continuous data/time data.table
    DTlist <- split(stationdata.DT.pyr, by = 'STATNUM')
    data_cont <- cont_date_time
    for (j in 1:length(DTlist)) {
      DTDummy <- as.data.table(DTlist[j])
      setnames(DTDummy, 1:numcols, cols)
      fileindx <- sprintf("%03d", j)
      if (j == 1) {
        setnames(DTDummy, cols[2:numcols], paste0(cols[2:numcols],'_',fileindx))
        setnames(DTDummy, 'PYROME_001', 'PYROME')        
      } else {
        DTDummy[,PYROME := NULL]
        numcolsin <- numcols - 1
        setnames(DTDummy, cols[2:numcolsin], paste0(cols[2:numcolsin],'_',fileindx))
      }
      data_cont <- merge(data_cont, DTDummy, by = 'DATE', all=TRUE)
    }

    # Extract RAIN *********************************************************************
    pattern <- c('RAIN')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]

    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.RAIN <- data_cont[, .(RAIN = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.RAIN)) set(data_cont.RAIN, which(is.nan(data_cont.RAIN[[j]])), j, 0)
    for(i in 2L){ set(data_cont.RAIN, i=NULL, j=i, value=1000*data_cont.RAIN[['RAIN']]) }
    setDT(data_cont.RAIN)
    # Extract RAIN *********************************************************************    

    # Extract RnDr ********************************************************************* 
    pattern <- c('RnDr')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]
    
    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.RnDr <- data_cont[, .(RnDr = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.RnDr)) set(data_cont.RnDr, which(is.nan(data_cont.RnDr[[j]])), j, 0)
    data_cont.RnDr[, c('RnDr') := lapply(.SD, function(x) (  round(x, 0)) ), .SDcols=c('RnDr')]
    # Extract RnDr ********************************************************************* 
    
    # Extract Temp ********************************************************************* 
    pattern <- c('TEMP')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]
    
    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.TEMP <- data_cont[, .(TEMP = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.TEMP)) set(data_cont.TEMP, which(is.nan(data_cont.TEMP[[j]])), j, 0)
    data_cont.TEMP[, c('TEMP') := lapply(.SD, function(x) (  round(x, 0)) ), .SDcols=c('TEMP')]    
    # Extract Temp *********************************************************************     

    # Extract ATMOIST (Basically RH) ********************************************************************* 
    pattern <- c('ATMOIST')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]
    
    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.ATMOIST <- data_cont[, .(ATMOIST = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.ATMOIST)) set(data_cont.ATMOIST, which(is.nan(data_cont.ATMOIST[[j]])), j, 0)
    data_cont.ATMOIST[, c('ATMOIST') := lapply(.SD, function(x) (  round(x, 0)) ), .SDcols=c('ATMOIST')]    
    # Extract ATMOIST (Basically RH) *********************************************************************   
    
    # Extract MaxT ********************************************************************* 
    pattern <- c('MaxT')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]
    
    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.MaxT <- data_cont[, .(MaxT = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.MaxT)) set(data_cont.MaxT, which(is.nan(data_cont.MaxT[[j]])), j, 0)
    data_cont.MaxT[, c('MaxT') := lapply(.SD, function(x) (  round(x, 0)) ), .SDcols=c('MaxT')]    
    # Extract MaxT *********************************************************************      
    
    # Extract MaxT ********************************************************************* 
    pattern <- c('MinT')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]
    
    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.MinT <- data_cont[, .(MinT = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.MinT)) set(data_cont.MinT, which(is.nan(data_cont.MinT[[j]])), j, 0)
    data_cont.MinT[, c('MinT') := lapply(.SD, function(x) (  round(x, 0)) ), .SDcols=c('MinT')]    
    # Extract MaxT *********************************************************************      
    
    # Extract MaxT ********************************************************************* 
    pattern <- c('MxRH')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]
    
    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.MxRH <- data_cont[, .(MxRH = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.MxRH)) set(data_cont.MxRH, which(is.nan(data_cont.MxRH[[j]])), j, 0)
    data_cont.MxRH[, c('MxRH') := lapply(.SD, function(x) (  round(x, 0)) ), .SDcols=c('MxRH')]    
    # Extract MaxT *********************************************************************       
    
    # Extract MaxT ********************************************************************* 
    pattern <- c('MnRH')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]
    
    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.MnRH <- data_cont[, .(MnRH = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.MnRH)) set(data_cont.MnRH, which(is.nan(data_cont.MnRH[[j]])), j, 0)
    data_cont.MnRH[, c('MnRH') := lapply(.SD, function(x) (  round(x, 0)) ), .SDcols=c('MnRH')]    
    # Extract MaxT *********************************************************************      

    # Merge final data
    DT_cmb <- data_cont.RAIN[data_cont.RnDr,][data_cont.TEMP,][data_cont.ATMOIST,]
    DT_cmb <- DT_cmb[data_cont.MaxT,][data_cont.MinT,][data_cont.MxRH,][data_cont.MnRH,]
    #setorder(data_cont.FIN, -RAIN) data_cont.RAIN

    # Output file with ammendments for testing
    filedir <- file.path(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedStationDataByPyrome'))
    fileout <- paste0('AvgdCombdStatDat_',PYROMEIN,'.rds')
    if (!dir.exists(filedir)){
      dir.create(filedir,recursive = TRUE)
    } else {
      print("Dir already exists!")
    }
    fileout <- paste0(filedir,'/',fileout)
    
    # Write concatenated data.table
    saveRDS(DT_cmb,fileout)
    
    # Clean memory
    rm(data_cont.RAIN,data_cont.RnDr,data_cont.TEMP,data_cont.ATMOIST,data_cont.MaxT,data_cont.MinT)
    rm(data_cont.MxRH,data_cont.MnRH,DT_cmb,data_cont,cont_date_time,stationdata.DT.pyr)    
    
  } # Select pyromes
} # All pyromes   

##################################################################################################################
## Section: Edit NARR grid location data .fwx formatted
##################################################################################################################

# File format ./ERCTimeSeriesAnalysis/Data/DownloadedGridData/README_GRIDDED.pdf
# Field | Field Name                                                                     | Columns 
# 1     Station Number (RRRCCC; RRR=row id from south, CCC=column id from west)          01-06 
# 2     Year                                                                             07-08
# 3     Month                                                                            09-10
# 4     Day                                                                              11-12
# 5     State of the Weather (0=clear; 3=overcast; 6=precipitation)                      13-13
# 6     Dry Bulb Temperature (F)                                                         14-16
# 7     Relative Humidity (%)                                                            17-19
# ****  BLANK NOT READ                                                                   20-27 
# 11    Wind Direction (8-pt)                                                            28-28
# 12    Wind Speed (mph)                                                                 29-31
# ****  BLANK NOT READ                                                                   32-38 
# 16    Maximum Temperature (F)                                                          39-41
# 17    Minimum Temperature (F)                                                          42-44 
# 18    Maximum RH (%)                                                                   45-47 
# 19    Minimum RH (%)                                                                   48-50
# 20    Season Code                                                                      51-51
# 21    Precipitation Duration (Hrs)                                                     52-53
# 22    Precipitation Amount (hundredths of an inch)                                     54-57
# ****  BLANK NOT READ                                                                   58-60 
# 24    Relative Humidity Variable (2=%)                                                 61-61

# Pass data from text DT to final DT
Headernames <- c('SIGName','YEAR','MONTH','DAY','SOW','TEMP','RH','BLNK1',
                 'WDIR','GSpd','BLNK2','MaxT','MinT','MxRH','MnRH','SEACD',
                 'RnDr','Rain','BLNK3','RHVAR')

colconvrt <- c('TEMP','RH','MaxT','MinT','MxRH','MnRH')

# Column data entered by hand from description above
end_col <- 61 #cumsum(cols)
start_col <- c(1,7,9,11,13,14,17,20,28,29,32,39,42,45,48,51,52,54,58,61) #end_col - start_col + 1
end_col <- c(6,8,10,12,13,16,19,27,28,31,38,41,44,47,50,51,53,57,60,61)
start_end <- cbind(start_col, end_col) # matrix of start and end positions
column_widths <- end_col - start_col + 1
COLWIDS <- as.data.frame(rbind(column_widths))
COLWIDS <- setDT(COLWIDS)
setnames(COLWIDS, old = 1:ncol(COLWIDS), new = Headernames)
COLWIDS <- as.list(COLWIDS)

# Read grid point for of all pyromes 
gridfile <- paste0('D:/ERCTimeSeriesAnalysis/Data/Pyromes_2019/','Top4CloseGridPointsAllPyromesCWAtrbs.rds')
griddeddata.DT.rnkd <- readRDS(gridfile)
sapply(griddeddata.DT.rnkd, class)
griddeddata.DT.rnkd[,filename := (as.character(filename))]
sapply(griddeddata.DT.rnkd, class)
griddeddata.DT.rnkd[, distrank := sprintf("%03d", as.integer(distrank))] 
sapply(griddeddata.DT.rnkd, class)

# Select subset of pyromes
#pyrssub <- c('001','005', '008', '009', '010', '011', '012', '013', '014', 
#             '015', '020', '026', '075', '101', '118', '122', '127')

pyrssub <- sprintf('%03d',seq(1:128))

# Specify output location
outloc <- paste0('Expmnt_',sprintf("%03d", as.integer(Exprmnt)),'/')

# Initialize output file of combined station data for each pyrome
griddeddata.DT <- data.table(NULL)

# Combined station data for all pyromes
griddeddata.DT.all.pyrs <- data.table(NULL)

# Select pyrome of interest to speed processing
#griddeddata.DT.rnkd <- griddeddata.DT.rnkd[PYROME == '001']

# Loop through griddeddata.DT.rnkd to process each NARR grid location file by Pyrome
for (iLine in 1:nrow(griddeddata.DT.rnkd)) { 
  # Evaluate pyrome code
  PYROMEIN <- unlist(griddeddata.DT.rnkd[iLine,PYROME])
  idx_one_in_two <- match(PYROMEIN, unlist(pyrssub)) 
  # Select subset of pyromes
  if (!is.na(idx_one_in_two)) {
    # Set working directory
    PYROMEIN <- pyrssub[idx_one_in_two]
    initial_path <- file.path(paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/','Pyromes_2019/GriddedWeatherData/',PYROMEIN))
    setwd(initial_path)
    wd <- getwd()
    print(paste0('Current working dir: ', wd))   
    
    # File to process for climate conditioning
    filenamein <- griddeddata.DT.rnkd[iLine,filename]
    distrankin <- griddeddata.DT.rnkd[iLine,distrank]
    filein <- paste0(initial_path,'/',distrankin,'_',filenamein)
    
    # Read data to data.table line by line
    grdfile.DT <- fread(filein, colClasses = "character",
                        sep = "\n", header = FALSE, verbose = TRUE)
    
    # Function to pass data by column widths
    textin <- lapply(grdfile.DT, function(x) {
      apply(start_end, 1, function(y) substr(x, y[1], y[2])) 
    })
    grdfile.DT <- data.table(textin$V1)
    setnames(grdfile.DT, old = 1:ncol(grdfile.DT), new = Headernames)
    
    # Drop blank rows remove whitespace
    grdfile.DT <- subset(grdfile.DT[SIGName  != "",])  
    grdfile.DT[, PYROME := PYROMEIN]

    # Append pyrome data to master data DT
    filesappend <- list(griddeddata.DT, grdfile.DT) 
    griddeddata.DT <- rbindlist(filesappend, use.names=TRUE, fill=TRUE) 
      
    # Append all pyromes 
    filesappend <- list(griddeddata.DT.all.pyrs, grdfile.DT) 
    griddeddata.DT.all.pyrs <- rbindlist(filesappend, use.names=TRUE, fill=TRUE) 
    rm(grdfile.DT)
  }
}

# Write concatenated data.table of all pyromes station data combined
cmbndgddatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedNARRGridDataByPyrome/AllPyromesOrigNARRGridDat.rds')
saveRDS(griddeddata.DT.all.pyrs, cmbndgddatfile)
#rm(griddeddata.DT.all.pyrs)

# Re-read combined input data
cmbndgddatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedNARRGridDataByPyrome/AllPyromesOrigNARRGridDat.rds')
griddeddata.DT.all.pyrs <- readRDS(cmbndgddatfile) 
sapply(griddeddata.DT.all.pyrs, class) 

# Reset data for all pyromes 
griddeddata.DT <- copy(griddeddata.DT.all.pyrs)
#griddeddata.DT <- griddeddata.DT[PYROME == '001']

# Update year, rain duration and amount to numeric fields
set_colclass(griddeddata.DT, c(YEAR='numeric', RnDr='numeric', Rain='numeric'))
class(griddeddata.DT$RnDr)
sapply(griddeddata.DT, class)

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
seltvar <- c('SIGName','DATE','RnDr','Rain','TEMP','RH','MaxT','MinT','MxRH','MnRH','PYROME','DATETIME')
griddeddata.DT <- griddeddata.DT[, ..seltvar]

# Enforce to be numeric
convrttype <- c('RnDr','Rain','TEMP','RH','MaxT','MinT','MxRH','MnRH')
mtimes <- rep(c('numeric'),each=length(convrttype))
mtimes <- setNames(mtimes, convrttype)
set_colclass(griddeddata.DT,mtimes)

# Convert decimal inches 00000, format nn.nnn to number
for (j in names(griddeddata.DT)) set(griddeddata.DT,which(is.na(griddeddata.DT[[j]])),j,0)
for(i in 4L){ set(griddeddata.DT, i=NULL, j=i, value=0.01*griddeddata.DT[['Rain']]) }

# Rename to be consistent with station data
oldnames <- c('SIGName','Rain')
newnames <- c('STATNUM','RAIN')
setnames(griddeddata.DT, old = oldnames, new = newnames)

# Convert pyrome station data to data frame for plotting
sapply(griddeddata.DT, class)
griddeddata.DF <- copy(griddeddata.DT)
griddeddata.DF <- setDF(griddeddata.DF)

ggplot() +
  # plot all simulations
  geom_line(data = griddeddata.DF, aes(x = DATE, y = RAIN, color = STATNUM, group=STATNUM), alpha = 1, size = 0.65) +
  labs(x="Date",y=expression('Precipitation Amount (inches)')) +
  theme(text = element_text(size=10))

ggplot(griddeddata.DF, aes(x = DATE, y = RAIN)) + 
  geom_line() + 
  facet_wrap(~ STATNUM, scales = 'free_y', ncol = 1)

##################################################################################################################
## Section: Clean stationdata.DT where days have multiple enetries but different times
##          Require closest data to 1300 hours, to represent day, but averages evaluated over day
##################################################################################################################

# Remove complete duplicates
setkey(griddeddata.DT, PYROME, STATNUM, DATETIME)
DUPFRST <- duplicated(griddeddata.DT, by = key(griddeddata.DT))
griddeddata.DT[, DUPSTDEL := DUPFRST | c(tail(DUPFRST, 0), FALSE)]

# Test duplicates
TESTDUPFRST <- griddeddata.DT[DUPSTDEL == TRUE]

# Remove dupliactes
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

# Write concatenated data.table of all pyromes station data combined duplicates removed
cmbndgddatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedNARRGridDataByPyrome/AllPyromesOrigNARRGridDat_DDUPED.rds')
saveRDS(griddeddata.DT, cmbndgddatfile)

##################################################################################################################
## Section: Create combined time series of precipitation for each pyrome from station data
##################################################################################################################

# Total number of pyromes
numpyr <- 128
pyromes <- seq(1,numpyr,1)

# Loop through pyromes 
griddeddata.DT[, PYROME := sprintf("%03d", as.integer(PYROME))] 
for (ipyr in pyromes) { 
  # Evaluate pyrome code
  PYROMEIN <- sprintf("%03d", as.integer(pyromes[ipyr]))
  idx_one_in_two <- match(PYROMEIN, unlist(pyrssub))
  # Select subset of pyromes
  if (!is.na(idx_one_in_two)) { 
    # Get data for this pyrome
    griddeddata.DT.pyr <- griddeddata.DT[PYROME == PYROMEIN]
    setcolorder(griddeddata.DT.pyr, 'DATE')
    cols <- names(griddeddata.DT.pyr)
    numcols <- length(cols)
    
    # Get min and max date values for complete timeseries
    mindate <- min(griddeddata.DT.pyr[,DATE])
    maxdate <- max(griddeddata.DT.pyr[,DATE])  
    cont_date_time <- as.data.table(data.frame(seq(as.Date(mindate), as.Date(maxdate), 'days')))
    setnames(cont_date_time, 1, 'DATE')
    
    # Split data.atble by station merge with continuous data/time data.table
    DTlist <- split(griddeddata.DT.pyr, by = 'STATNUM')
    data_cont <- cont_date_time
    for (j in 1:length(DTlist)) {
      DTDummy <- as.data.table(DTlist[j])
      setnames(DTDummy, 1:numcols, cols)
      fileindx <- sprintf("%03d", j)
      if (j == 1) {
        setnames(DTDummy, cols[2:numcols], paste0(cols[2:numcols],'_',fileindx))
        setnames(DTDummy, 'PYROME_001', 'PYROME')        
      } else {
        DTDummy[,PYROME := NULL]
        numcolsin <- numcols - 1
        setnames(DTDummy, cols[2:numcolsin], paste0(cols[2:numcolsin],'_',fileindx))
      }
      data_cont <- merge(data_cont, DTDummy, by = 'DATE', all=TRUE)
    }
    
    # Extract RAIN *********************************************************************
    pattern <- c('RAIN')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]
    
    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.RAIN <- data_cont[, .(RAIN = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.RAIN)) set(data_cont.RAIN, which(is.nan(data_cont.RAIN[[j]])), j, 0)
    for(i in 2L){ set(data_cont.RAIN, i=NULL, j=i, value=1000*data_cont.RAIN[['RAIN']]) }
    setDT(data_cont.RAIN)
    # Extract RAIN *********************************************************************    
    
    # Extract RnDr ********************************************************************* 
    pattern <- c('RnDr')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]
    
    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.RnDr <- data_cont[, .(RnDr = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.RnDr)) set(data_cont.RnDr, which(is.nan(data_cont.RnDr[[j]])), j, 0)
    data_cont.RnDr[, c('RnDr') := lapply(.SD, function(x) (  round(x, 0)) ), .SDcols=c('RnDr')]
    # Extract RnDr ********************************************************************* 
    
    # Extract Temp ********************************************************************* 
    pattern <- c('TEMP')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]
    
    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.TEMP <- data_cont[, .(TEMP = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.TEMP)) set(data_cont.TEMP, which(is.nan(data_cont.TEMP[[j]])), j, 0)
    data_cont.TEMP[, c('TEMP') := lapply(.SD, function(x) (  round(x, 0)) ), .SDcols=c('TEMP')]    
    # Extract Temp *********************************************************************     
    
    # Extract ATMOIST (Basically RH) ********************************************************************* 
    pattern <- c('ATMOIST')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]
    
    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.ATMOIST <- data_cont[, .(ATMOIST = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.ATMOIST)) set(data_cont.ATMOIST, which(is.nan(data_cont.ATMOIST[[j]])), j, 0)
    data_cont.ATMOIST[, c('ATMOIST') := lapply(.SD, function(x) (  round(x, 0)) ), .SDcols=c('ATMOIST')]    
    # Extract ATMOIST (Basically RH) *********************************************************************   
    
    # Extract MaxT ********************************************************************* 
    pattern <- c('MaxT')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]
    
    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.MaxT <- data_cont[, .(MaxT = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.MaxT)) set(data_cont.MaxT, which(is.nan(data_cont.MaxT[[j]])), j, 0)
    data_cont.MaxT[, c('MaxT') := lapply(.SD, function(x) (  round(x, 0)) ), .SDcols=c('MaxT')]    
    # Extract MaxT *********************************************************************      
    
    # Extract MaxT ********************************************************************* 
    pattern <- c('MinT')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]
    
    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.MinT <- data_cont[, .(MinT = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.MinT)) set(data_cont.MinT, which(is.nan(data_cont.MinT[[j]])), j, 0)
    data_cont.MinT[, c('MinT') := lapply(.SD, function(x) (  round(x, 0)) ), .SDcols=c('MinT')]    
    # Extract MaxT *********************************************************************      
    
    # Extract MaxT ********************************************************************* 
    pattern <- c('MxRH')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]
    
    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.MxRH <- data_cont[, .(MxRH = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.MxRH)) set(data_cont.MxRH, which(is.nan(data_cont.MxRH[[j]])), j, 0)
    data_cont.MxRH[, c('MxRH') := lapply(.SD, function(x) (  round(x, 0)) ), .SDcols=c('MxRH')]    
    # Extract MaxT *********************************************************************       
    
    # Extract MaxT ********************************************************************* 
    pattern <- c('MnRH')
    tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(data_cont))))
    indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
    cols <- colnames(data_cont)[indxs]
    
    # Fix zero elements as NA
    data_cont[, (cols) :=  lapply(.SD, function(x) ifelse(x <=0, NA, x)), .SDcols = cols]
    data_cont.MnRH <- data_cont[, .(MnRH = rowMeans(.SD, na.rm=TRUE)), .SDcols = cols, by = DATE]
    for (j in 1:ncol(data_cont.MnRH)) set(data_cont.MnRH, which(is.nan(data_cont.MnRH[[j]])), j, 0)
    data_cont.MnRH[, c('MnRH') := lapply(.SD, function(x) (  round(x, 0)) ), .SDcols=c('MnRH')]    
    # Extract MaxT *********************************************************************      
    
    # Merge final data
    DT_cmb <- data_cont.RAIN[data_cont.RnDr,][data_cont.TEMP,][data_cont.ATMOIST,]
    DT_cmb <- DT_cmb[data_cont.MaxT,][data_cont.MinT,][data_cont.MxRH,][data_cont.MnRH,]
    #setorder(data_cont.FIN, -RAIN) data_cont.RAIN
    
    # Output file with ammendments for testing
    filedir <- file.path(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedNARRGridDataByPyrome'))
    fileout <- paste0('AvgdCombdNARRGridDat_',PYROMEIN,'.rds')
    if (!dir.exists(filedir)){
      dir.create(filedir,recursive = TRUE)
    } else {
      print("Dir already exists!")
    }
    fileout <- paste0(filedir,'/',fileout)
    
    # Write concatenated data.table
    saveRDS(DT_cmb,fileout)
    
    # Clean memory
    rm(data_cont.RAIN,data_cont.RnDr,data_cont.TEMP,data_cont.ATMOIST,data_cont.MaxT,data_cont.MinT)
    rm(data_cont.MxRH,data_cont.MnRH,DT_cmb,data_cont,cont_date_time,stationdata.DT.pyr)    
    
  } # Select pyromes
} # All pyromes 

##################################################################################################################
## Section: Analyze differences between averaged station and NARR gridded data
##################################################################################################################




# Select files to process
filedir <- file.path(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedStationDataByPyrome'))
setwd(filedir)
wd <- getwd()
print(paste0('Current working dir: ', wd))  

# Select files to process
pattern <- c('AvgdCombdStatDat')
filesin <- list.files(pattern=paste0(pattern, collapse="|"),recursive = TRUE, full.names = TRUE)  
filenamesin <- basename(filesin)
rdsfilesin <- include(filenamesin,'.rds')
rdsfilesin <- rdsfilesin[order(as.integer(substring(rdsfilesin, 5, 8)),  decreasing = TRUE, rdsfilesin)]
rdsfilesin.srtd <- rdsfilesin[!duplicated(as.integer(substring(rdsfilesin, 5, 8)))]
if (length(rdsfilesin.srtd) >= 4) {
  rdsfilesin <- rdsfilesin.srtd[-c(4:length(rdsfilesin.srtd))]
} else {
  rdsfilesin <- rdsfilesin.srtd   
}




# Plot ggplot 
legend_title <- 'Climate Type'
p <- ggplot(stationdata.DT.plt.DF, aes(x = DATE, y = RAIN)) + 
  geom_line() + 
  facet_wrap(~ STATNUM, scales = 'free_y', ncol = 1)
png(width=1000, height=500)
print(p)
dev.off()
ggsave(paste0('Pyrome_',PYROMEIN,'_all.png'), p, width = 9, height = 9)

##################################################################################################################
## Section: Define miscellaneous functions
##################################################################################################################

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
##            Program end section to Process weather weather station data in .fw13 format                       ##
##                                                                                                              ##
##################################################################################################################






