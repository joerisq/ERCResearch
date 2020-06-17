##################################################################################################################
## Project: Prepare final climate conditioned averaged station weather fiels for FireFamilplus use
##
## Script purpose: Ingest cliamte conditioned averaged station files, reformat as old .fwx formatted files
##                 for import to FireFamily plus software
##
## Date: 19th March 2020
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
ROOT_DATA_DIR <- 'C:/ERCTimeSeriesAnalysis'
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
VersonOut <- 7

##################################################################################################################
## Section: Read climate conditioned averaged station data per pyrome in raw format and station list 
##################################################################################################################

# Read concatenated data.table of all pyromes combined - list of data.frames
#cmbndstnrdatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/AvgdCombdStationNARR.rds')
cmbndstnrdatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/IDW_AvgdCombdStationNARR.rds')
wthrdatallpyrs.DT <- readRDS(cmbndstnrdatfile) 
#sapply(wthrdatallpyrs.DT, class) 

# Read experimental station data to clean ERC output 
#cmbndstnrdatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/AvgdCombdStation.rds')
#wthrdatallpyrs.DT <- readRDS(cmbndstnrdatfile) 

# Read experimental NARR grid data to celan ERC output 
#cmbndstnrdatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/AvgdCombdNARR.rds')
#wthrdatallpyrs.DT <- readRDS(cmbndstnrdatfile) 

# Read concatenated data.table of all pyromes and closest stations
stationlistfile <- paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/','Pyromes_2019/','Top5CloseStationsAllPyromes_',DATA_YR_OUT,'.rds')
stationlistfile.dt.rnkd.master <- readRDS(stationlistfile)
stationlistfile.dt.rnkd.master <- stationlistfile.dt.rnkd.master[NESDIS != '________']
stationlistfile.dt.rnkd.master <- stationlistfile.dt.rnkd.master[STATID != 149901]
stationlistfile.dt.rnkd.master[PYROME == '088']

# Read Pyrome shapefile to calculate Pyrome centroids 
Pyrms_shp <- readOGR(paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/','Pyrome_20150605.shp')) 
Pyrms_shp <- spTransform(Pyrms_shp, CRS('+proj=longlat +ellps=WGS84 +datum=WGS84'))
centroidS <- gCentroid(Pyrms_shp,byid=TRUE)
leg_centers <- SpatialPointsDataFrame(gCentroid(Pyrms_shp, byid=TRUE), 
                                      Pyrms_shp@data, match.ID=FALSE)
# Create data frame from centroids
leg_centers.cntds <- data.frame(leg_centers@data, x = coordinates(leg_centers)[,1], y = coordinates(leg_centers)[,2])
leg_centers.cntds <- subset(leg_centers.cntds, select=c(PYROME, x, y))
setDT(leg_centers.cntds)
setnames(leg_centers.cntds, c('x','y'),c('LONGITUDE','LATITUDE'))

##################################################################################################################
## Section: Reformat data to .fwx formatted
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

# Set climate condition flag, TRUE implies to create climate conditioned output
CCDFlag <- TRUE
if (CCDFlag == FALSE) {
  filsufx <- 'BASE'
} else {
  filsufx <- 'CCD'
}

# Load in file containing the climate conditioning change factors and risQ 2-deg grid information
#   - v1: contains preliminary conditioning of temperature and RH
#   - files are on AWS
trgtyr <- 2060
percte <- 0.5
percte <- 'med'
scnio <- 'rcp85'
runspec <- paste0(filsufx,'_',trgtyr,percte,scnio)
chngfacts <- paste0('change_factors_',trgtyr,'_Q_',percte)
change_factors_tas <- readRDS(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/ClimateConditioningInput/Climate_Conditioning_Change_Factors_tas_v4.rds'))
change_factors_pr <- readRDS(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/ClimateConditioningInput/Climate_Conditioning_Change_Factors_pr_v4.rds'))
change_factors_rh <- readRDS(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/ClimateConditioningInput/Climate_Conditioning_Change_Factors_rh_v4.rds'))

# Month sequence
mnthsin <- seq(1,12,1)

# Copy Firefamily plus input templates C:\ERCTimeSeriesAnalysis\Output\Templates
softwareflesloc <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/Templates')
softwareflesin.mstr <- list.files(softwareflesloc,pattern=('*.*'),recursive = TRUE, full.names = TRUE) 

# Loop through each pyrome with list of pyrome data.frames 
for (i in 1:length(wthrdatallpyrs.DT)) { 

  # Specify pyrome as character
  pyrnms <- names(wthrdatallpyrs.DT)
  pyrmin <- substring(pyrnms[i],8,10)
  pyrmin <- as.integer(pyrmin)
  PYROMEIN <- sprintf("%03d", as.integer(pyrmin))

  # Create data.table for current pyrome data 
  wthrdat.dt <- as.data.table(wthrdatallpyrs.DT[[pyrmin]])
  setorder(wthrdat.dt, DATE)
  
  # Rename odd field names
  setnames(wthrdat.dt, 'ATMOIST', 'RH')
  setnames(wthrdat.dt, 'RAIN', 'Rain')
  
  # Split date field into Year, month and day fields
  wthrdat.dt[, YEAR := as.numeric(format(DATE, format = "%y"))] 
  wthrdat.dt[, YEAR := sprintf("%02d", as.integer(YEAR))]   
  wthrdat.dt[, MONTH := sprintf("%02d", as.integer(month(DATE)))]
  wthrdat.dt[, DAY := sprintf("%02d", as.integer(day(DATE)))]
  wthrdat.dt[, DATE := NULL]

  drop_rows_all_na <- function(x, pct=0.01) x[!rowSums(is.na(x)) >= ncol(x)*pct,]
  wthrdat.dt <- drop_rows_all_na(wthrdat.dt,0.001)
  
  ##################################### COMPARE ORIGINAL AND FINAL PYROME DATA  
  # open original file
  #origpyrdta <- readRDS(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/AvgdCombdStatDat_',PYROMEIN,'.rds'))
  #setorder(origpyrdta, DATE)
  
  # Rename odd field names
  #setnames(origpyrdta, 'ATMOIST', 'RH')
  #setnames(origpyrdta, 'RAIN', 'Rain')
  # Split date field into Year, month and day fields
  #origpyrdta[, YEAR := as.numeric(format(DATE, format = "%y"))] 
  #origpyrdta[, YEAR := sprintf("%02d", as.integer(YEAR))]   
  #origpyrdta[, MONTH := sprintf("%02d", as.integer(month(DATE)))]
  #origpyrdta[, DAY := sprintf("%02d", as.integer(day(DATE)))]
  #origpyrdta[, DATE := NULL]

  #drop_rows_all_na <- function(x, pct=0.01) x[!rowSums(is.na(x)) >= ncol(x)*pct,]
  #origpyrdta <- drop_rows_all_na(origpyrdta,0.001)  
  #rm(wthrdat.dt)
  #wthrdat.dt <- origpyrdta
  ##################################### COMPARE ORIGINAL AND FINAL PYROME DATA  
  
  if (CCDFlag == TRUE) {
  ##################################### CLIMATE CONDITION PYROME DATA
    # Get pyrome centroids
    pyrlon <- leg_centers.cntds[PYROME == pyrmin, LONGITUDE]
    pyrlat <- leg_centers.cntds[PYROME == pyrmin, LATITUDE]  
  
    cc_grdfile.DT <- condition_pyrome_3.1_V01(plat = pyrlat, plon = pyrlon, Year = eval(trgtyr), scenario = eval(scnio),
                                    scenario_obj = eval(percte), pyrome_data = wthrdat.dt, change_factors_tas = change_factors_tas,
                                    change_factors_pr = change_factors_pr, change_factors_rh = change_factors_rh)
    # Copy output
    wthrdat.dt <- copy(cc_grdfile.DT)
  }
  ##################################### CLIMATE CONDITION PYROME DATA
  
  # Initailize data.table to contain output Station id, station name, longitude, latitude and state name  
  colnms <- names(wthrdat.dt)
  coldifs <- as.vector(setdiff(Headernames,colnms))
  rowsin <- nrow(wthrdat.dt)
  wthrdat.dt.dmy <- data.table(1)[,`:=`(c(coldifs),'')][,V1:=NULL]
  wthrdat.dt.dmy <- rbind(wthrdat.dt.dmy, list(id = 1:(rowsin - 1)), fill = TRUE)
  collength <- dim(wthrdat.dt.dmy)[2] 
  Cols <- colnames(wthrdat.dt.dmy)
  convert_to_char <- Cols
  wthrdat.dt.dmy[, c(convert_to_char) := lapply(.SD, as.character), .SDcols=convert_to_char]
  sapply(wthrdat.dt.dmy, class)
  wthrdat.dt.dmy[, id := NULL]
  
  # Bind data.table containing differenced columns with input data
  wthrdat.dt <- do.call(cbind, lapply(list(wthrdat.dt.dmy, wthrdat.dt), setDT))
  setcolorder(wthrdat.dt, Headernames)
  Cols <- colnames(wthrdat.dt) 
  
  # Enforce to be numeric 
  convrttype <- c('RnDr','Rain','TEMP','RH','MaxT','MinT','MxRH','MnRH',
                  'WDIR','GSpd','SEACD','RHVAR','SOW')
  mtimes <- rep(c('numeric'),each=length(convrttype))
  mtimes <- setNames(mtimes, convrttype)
  set_colclass(wthrdat.dt,mtimes)
  
  # Fix RH and MxRH/MnRH to be consistent
  wthrdat.dt[,MnRH := MnRH][MnRH > RH, MnRH := RH]
  #wthrdat.dt[,RH := RH][MnRH > RH, RH := MnRH]  
  #wthrdat.dt[,MnRH := MnRH][MnRH > RH, MnRH := MnRH + (RH - MnRH)]
  #wthrdat.dt[,MnRH := MnRH][MnRH < 0, MnRH := RH]   

  # Fix SOW code based on existing NARR input files (https://wrcc.dri.edu/fpa/README_GRIDDED.pdf)
  wthrdat.dt[, SOW := 0] 
  wthrdat.dt[, SOW := SOW][Rain > 0.02, SOW := 6]
  
  # Set Relative Humidity Variable based on existing NARR input files
  wthrdat.dt[, RHVAR := 2]   
  
  # Season code interpreted as quarterly code
  wthrdat.dt[, SEACD := paste0('',ceiling(as.numeric(MONTH) / 3))]
  
  # Set wind speed and direction
  wthrdat.dt[, WDIR := 1]   
  wthrdat.dt[, GSpd := 30]  
  
  # Convert decimal rain from inches nn.nnn to format 0000 100s of inches, padded with zeros
  for(i in 18L){ set(wthrdat.dt, i=NULL, j=i, value=wthrdat.dt[['Rain']]/0.01) }
  wthrdat.dt[, Rain := trimws(format(ceiling(as.numeric(Rain))))]

  # Set station name
  stattype <- stationlistfile.dt.rnkd.master[PYROME == PYROMEIN][1, STATID]
  if (PYROMEIN == '017' ) {stattype = '040239'} 
  if (PYROMEIN == '056' ) {stattype = '053104'}   
  if (PYROMEIN == '018' ) {stattype = '040630'}   
  if (PYROMEIN == '024' ) {stattype = '421806'}   
  if (PYROMEIN == '010' ) {stattype = '241502'}    
  if (PYROMEIN == '041' ) {stattype = '020121'}     
  if (PYROMEIN == '053' ) {stattype = '418701'}   
  if (PYROMEIN == '069' ) {stattype = '343301'}  
  if (PYROMEIN == '087' ) {stattype = '238601'}   
  if (PYROMEIN == '084' ) {stattype = '411401'}     
  if (PYROMEIN == '096' ) {stattype = '210509'}   
  if (PYROMEIN == '100' ) {stattype = '202701'}  
  if (PYROMEIN == '094' ) {stattype = '231301'}     
  if (PYROMEIN == '121' ) {stattype = '336001'}   
  if (PYROMEIN == '126' ) {stattype = '409201'}    
  
  wthrdat.dt[, SIGName := stattype]    
  
  # Set blank fields to be blank
  pattern <- c('BLNK')
  tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, colnames(wthrdat.dt))))
  indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE))))
  cols <- colnames(wthrdat.dt)[indxs]
  for (j in names(wthrdat.dt)[indxs])set(wthrdat.dt,which(is.na(wthrdat.dt[[j]])),j,'')

  ##################################### INCLUDE RANDOM WIND DATA ###############  
  # Include random wind speeds
  filedir <- file.path(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/'))
  fileout <- paste0(filedir,'/','RndmWindsByPyrome',PYROMEIN,'.rds')
  winddat.DT <- readRDS(fileout) 
  
  # Columns to keep
  colstokeep <- c('MONTH','WD','WSDMY')
  winddat.DT <- winddat.DT[ , ..colstokeep]

  drop_rows_all_na <- function(x, pct=0.01) x[!rowSums(is.na(x)) >= ncol(x)*pct,]
  winddat.DT <- drop_rows_all_na(winddat.DT,0.001)
  
  # For .fwx formatted files require WD and WSDMY fields
  wndstat <- wthrdat.dt[ , list(COUNT = .N) , by = .(MONTH)]
  wthrdat.dt[ , `:=`( IDX = 1:.N ) , by = MONTH]
  setkey(wthrdat.dt, MONTH, IDX)

  # Loop through month, sample from each
  winddat.DT.fl.yr <- data.table(NULL)
  for (imnth in mnthsin) { 
    mnthin <- sprintf('%02d', as.integer(imnth))
    sizebymonth <- wndstat[MONTH == mnthin, COUNT]
    winddat.DT.mnth <- winddat.DT[MONTH == imnth]
    set.seed(imnth)
    winddat.DT.mnth.smpl <- winddat.DT.mnth[, .SD[sample(x = .N, size = sizebymonth)], by = MONTH]
    winddat.DT.mnth.smpl[ , `:=`( IDX = 1:.N ) , by = MONTH]
    # Append all months for current pyrome 
    if (imnth == 1) {
      winddat.DT.fl.yr <- winddat.DT.mnth.smpl
    } else {
      filesappend <- list(winddat.DT.fl.yr, winddat.DT.mnth.smpl) 
      winddat.DT.fl.yr <- rbindlist(filesappend, use.names=TRUE, fill=TRUE)  
    }
  }
  winddat.DT.fl.yr[,MONTH := sprintf("%02d", as.integer(MONTH))]
  
  # Merge winD data with all other variables
  setkey(winddat.DT.fl.yr, MONTH, IDX)
  wthrdat.dt <- wthrdat.dt[winddat.DT.fl.yr]
  wthrdat.dt[, WDIR := WD]   
  wthrdat.dt[, GSpd := WSDMY] 
  wthrdat.dt[ ,`:=`(WD = NULL, WSDMY = NULL, IDX = NULL)]
  ##################################### INCLUDE RANDOM WIND DATA ###############  

  # Replace all other NA etc with 0
  for (j in names(wthrdat.dt))set(wthrdat.dt,which(is.na(wthrdat.dt[[j]])),j,0)  

  # Round data to integers and remove whitespace
  convrttype <- convrttype[convrttype != 'Rain']
  wthrdat.dt[, c(convrttype) := lapply(.SD, function(x) trimws(format(ceiling(as.numeric(x)))) ), .SDcols=c(convrttype)]
  
  # Check column widths of output data, pad if incorrect
  colconvrt <- names(wthrdat.dt)
  wthrdat.dt[ , c(colconvrt) := Map(function(x, COLWIDS) {str_pad(x, COLWIDS, side = c('left'), pad = " ")},
                                     .SD, COLWIDS), .SDcols = c(colconvrt)]
  
  # Pad Rain field with zeros
  wthrdat.dt[, Rain := as.integer(Rain)]  
  wthrdat.dt[, Rain := str_pad(Rain, 4, side = c('left'), pad = '0')] 

  # Bind output columns as 1 column
  wthrdat.dt[, key_ := do.call(paste, c(.SD, sep = "")), .SDcols = names(wthrdat.dt)]
    
  # Specify output location
  outloc <- paste0('Version_V',sprintf("%02d", as.integer(VersonOut)),'/')
  filedir <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/',outloc,PYROMEIN,'/risQWthrFriskData') 
  fileout <- paste0('IDW_AvgdStatDat_',runspec,'_',PYROMEIN,'.fwx')  
  if (!dir.exists(filedir)){
    dir.create(filedir,recursive = TRUE)
  } else {
    print("Dir already exists!")
  }
  fileout <- paste0(filedir,'/',fileout)
  
  # Sort output on date
  wthrdat.dt[, DATE := paste0(YEAR,'-',MONTH,'-',DAY)] 
  wthrdat.dt[, DATE := parse_date_time(DATE, 'ymd')]
  setorder(wthrdat.dt, DATE)

  # Use fixed format, suppress row names and column names
  data <- wthrdat.dt[,key_ ]
  write.table(data, fileout, quote = FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  
  # Copy Firefamily input file templates to each directory - Reanme files before copy
  softwareflesin <- softwareflesin.mstr
  softwareflesinrnmd <- list.files(filedir,pattern=(PYROMEIN),recursive = TRUE, full.names = TRUE)
  if(length(file.exists(softwareflesinrnmd)) < 2){ 
    cat('Copying templates: ', PYROMEIN, '\n')
    for(file in softwareflesin) {
      file.copy(softwareflesin, filedir)
    }
    softwareflesinrnmd <- list.files(filedir,pattern='000',recursive = TRUE, full.names = TRUE)
    softwareflesrnm <- gsub('000',PYROMEIN, softwareflesinrnmd)
    runspec <- paste0(filsufx,trgtyr,percte,scnio)
    softwareflesrnm <- gsub(filsufx,runspec, softwareflesrnm)     
    file.rename(softwareflesinrnmd, softwareflesrnm)
  }

} # List of pyrome data.frames

##################################################################################################################
## Section: Define miscellaneous functions
##################################################################################################################

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
##            Program end section to Process weather inputs to carry out trend analysis                         ##
##                                                                                                              ##
##################################################################################################################