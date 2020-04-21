##################################################################################################################
## Project: Compile station description data files for input to fire simualator FSIM
##
## Script purpose: Process weather inputs to carry out FireFamilyPls5 analysis
##
## Date: 8th November 2019
## Author: JR
##################################################################################################################

# Required testing many libraries. Be selective in initializing libraries
library(rlang) 
library(maptools) 
library(rgdal) 
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
library(raster) 
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
ROOT_DATA_DIR <- ('C:/ERCTimeSeriesAnalysis')
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

##################################################################################################################
## Section: Read all station description files from FAM.NWCG.GOVFAM-files
##################################################################################################################

# Check directory
initial_path <- file.path(paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/','FAM.NWCG.GOVFAM-files-',DATA_YR))
setwd(initial_path)
wd <- getwd()
print(paste0('Current working dir: ', wd))

# Collect names of file of interest - wlstinv1!500100.txt individual station index files
filesin <- list.files(pattern = '^wlstinv1!', recursive = TRUE, full.names = TRUE)
filenamesin <- basename(list.files(pattern = "^wlstinv1!", recursive = TRUE))
numfiles <- length(filenamesin)

# Initailize data.table to contain output Station id, station name, longitude, latitude and state name
sationdata.dt <- data.table(1)[,`:=`(c('STATINDEX', 'STATID', 'STATNAME', 'NESDIS', 'LONGITUDE', 'LATITUDE', 'STATENAME', 'FLNAME'),0)][,V1:=NULL]
sationdata.dt <- sationdata.dt[rep(1:.N,numfiles)][,STATINDEX:=1:.N]
collength <- dim(sationdata.dt)[2] 
Cols <- colnames(sationdata.dt)
convert_to_char <- colnames(sationdata.dt)[c(2:4,7,collength)]
sationdata.dt[, c(convert_to_char) := lapply(.SD, as.character), .SDcols=convert_to_char]
sapply(sationdata.dt, class)

# Loop through indexfiles collecting Station names, IDs and geocodes
readsizeof <- 20000
nooflines <- 0
for (j in 1:length(filesin)) { 
  # Create connection
  file <- filesin[j]
  FLNAMEIN <- filenamesin[j]
  statein <- substring(file,3,4)
  
  # Store station attributes in sationdata.dt  
  sationdata.dt[STATINDEX  == j, STATENAME := statein] 
  sationdata.dt[STATINDEX  == j, FLNAME := FLNAMEIN]   
  
  if (statein != 'ak' && statein != 'hi' && statein != 'pr')  {
  con <- file(description=file, open="r")

  # Test read to get filesize
  #while((linesread <- length(readLines(con,readsizeof))) > 0 ) 
    #nooflines <- nooflines+linesread
  
  linesreadforfile <- readLines(con)
  print(j)
  
  # Store station attributes in sationdata.dt
  linein <- linesreadforfile[4]
  linein <- substring(linein, regexpr(":", linein) + 1)
  STATIDIN <- trimws(substring(linein,2,8))
  sationdata.dt[STATINDEX  == j, STATID := STATIDIN]
  
  linein <- substring(linein, regexpr(":", linein) + 1)
  STATNAMEIN <- trimws(substring(linein,2,22))
  sationdata.dt[STATINDEX  == j, STATNAME := STATNAMEIN]  
  
  linein <- substring(linein, regexpr(":", linein) + 1)
  NESDISIN <- trimws(substring(linein,2,34))
  sationdata.dt[STATINDEX  == j, NESDIS := NESDISIN]   
  
  linein <- linesreadforfile[6]
  linein <- substring(linein, regexpr("Obs", linein))
  TZONEIN <- trimws(linein)
  sationdata.dt[STATINDEX  == j, TZONE := TZONEIN]

  linein <- linesreadforfile[9]
  linein <- substring(linein, regexpr("Lat/Lon:", linein) + 8)
  bits <- trimws(unlist(strsplit(linein, ',')))
  bits <- gsub('[[:punct:]]', '', bits)
  bits <- gsub('  ', ' ', bits)
  bitstest <- grepl("^[[:digit:]]",gsub(' ', '', bits))
  if(is.na(bitstest[1])) {bitstest[1] = FALSE}
  if(is.na(bitstest[2])) {bitstest[2] = FALSE}  
  LonDDMMSS <- 0.0
  LatDDMMSS <- 0.0

  if (bitstest[1] == TRUE && bitstest[2] == TRUE) {
    linein <- linesreadforfile[9]
    linein <- substring(linein, regexpr("Lat/Lon:", linein) + 8)
    bits <- trimws(unlist(strsplit(linein, ',')))
    bits <- gsub('[[:punct:]]', '', bits)
    bits <- gsub('  ', ' ', bits) 
    LatDDMMSS <- angle2dec(bits[1],' ')
    LonDDMMSS <- -angle2dec(bits[2],' ')
  } else if (bitstest[1] == FALSE && bitstest[2] == FALSE) {
    linein <- linesreadforfile[11]
    linein <- substring(linein, regexpr("Lat/Lon:", linein) + 8)
    bits <- trimws(unlist(strsplit(linein, ',')))
    bits <- gsub('[[:punct:]]', '', bits)
    bits <- gsub('  ', ' ', bits) 
    bitstest <- grepl("^[[:digit:]]",gsub(' ', '', bits))
    if (bitstest[1] == TRUE && bitstest[2] == TRUE) {
      LatDDMMSS <- angle2dec(bits[1],' ')
      LonDDMMSS <- -angle2dec(bits[2],' ')      
    } else {
      sationdata.dt[STATINDEX  == j, LONGITUDE := 0] 
      sationdata.dt[STATINDEX  == j, LATITUDE := 0]        
    }
  } else {
    sationdata.dt[STATINDEX  == j, LONGITUDE := 0] 
    sationdata.dt[STATINDEX  == j, LATITUDE := 0]   
  }
  sationdata.dt[STATINDEX  == j, LONGITUDE := LonDDMMSS] 
  sationdata.dt[STATINDEX  == j, LATITUDE := LatDDMMSS]   
  close(con)
  }
}

# Prepare data for output
stationlistfile.rds <- paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/','FAM.NWCG.GOVFAM-files_StationList_',DATA_YR_OUT,'.rds')
stationlistfile.csv <- paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/','FAM.NWCG.GOVFAM-files_StationList_',DATA_YR_OUT,'.csv')
saveRDS(sationdata.dt, stationlistfile.rds)
fwrite(sationdata.dt, stationlistfile.csv)

##################################################################################################################
## Section: Define miscellaneous functions
##################################################################################################################

# Convert to decimal degrees proper
angle2dec <- function(angle,pattern) {
  angle <- as.character(angle)
  x <- do.call(rbind, strsplit(angle, pattern, perl = TRUE))
  x <- apply(x, 1L, function(y) {
    y <- as.numeric(y)
    y[1] + y[2]/60 + y[3]/3600
  })
  return(x)
}

##################################################################################################################
# Use the following sites to download station data. Fisrt site is most complete
## https://fam.nwcg.gov/fam-web/weatherfirecd/state_data.htm #
## https://cefa.dri.edu/raws/index.php

# Use the following sites to download gridded data.
## https://wrcc.dri.edu/fpa/

##################################################################################################################
## Section: Book keeping - Clean memory close file connections
##################################################################################################################

# CLEAN MEMORY
rm(list = ls(all.names = TRUE))
raster::removeTmpFiles(h = 0)
flush.console()

##################################################################################################################
##                                                                                                              ##
##            Program end section toProcess weather inputs to carry out FireFamilyPls5 analysis                 ##
##                                                                                                              ##
##################################################################################################################