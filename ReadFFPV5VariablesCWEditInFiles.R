##################################################################################################################
## Project: Climate conditioning ERC stream input to fire simualator FSIM
##
## Script purpose: Read weather variables created by FireFamilyplus v5 software to feed climate conditioning
##                 process
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
## Section: Read FFP5 weather output files for each Pyrome - PRE CLIMATE CONDITIONING
##################################################################################################################

# Check directory
initial_path <- file.path(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/','Pyromes'))
setwd(initial_path)
wd <- getwd()
print(paste0('Current working dir: ', wd))

# Collect names of file of interest - wlstinv1!500100.txt individual station index files
filesin <- list.files(pattern = 'FFP5Weather', recursive = TRUE, full.names = TRUE)
filenamesin <- basename(list.files(pattern = 'FFP5Weather', recursive = TRUE))
numfiles <- length(filenamesin)

# Initialize data.table to hold data for all pyromes 
wthrvars.DT <- data.table(NULL)

# Loop through indexfiles collecting weather variables by pyrome
for (j in 1:length(filesin)) { 
  # Create connection
  file <- filesin[j]
  FLNAMEIN <- filenamesin[j]
  pyromein <- substring(file,3,5)

  # Read data to data.table line by line
  ERC.DT <- fread(file, colClasses = "character",
                  sep = ',', header = FALSE, verbose = TRUE)
  
  # Delete meta data and header followed by dummy column 49
  ERC.DT[, 49 := NULL]
  
  # Pass data from text DT to final DT
  Headernames <- c('SIGName','DATE','Time','Temp','AvgT','MinT','MaxT','RH','AvRH','MnRH','MxRH','Rain','RnDr','Wind',
  'SC','ERC','BI','KBDI','IC','FM1','FM10','F100','F1000','FMH','FMW','SFlag','WDIR','SOW','FFMC',
  'DMC','DC','ISI','BUI','FWI','DSR','GDir','GSpd','SolR','WFlag','DPT','VPDM','VPDA','GSI','WAZI',
  'HrRa','LHerb','LWood','FFWI')
  
  oldnames <- names(ERC.DT)
  setnames(ERC.DT, old = oldnames, new = Headernames, skip_absent=TRUE)
  
  # Add Pyrome reference column
  ERC.DT[,PYROME := pyromein]
  
  # Append pyrome data to master data DT
  if (j == 1) {
    wthrvars.DT <- ERC.DT
  } else {
    # Write top 4 grid locations
    filesappend <- list(wthrvars.DT, ERC.DT) 
    wthrvars.DT <- rbindlist(filesappend, use.names=TRUE, fill=TRUE)  
  }
  rm(ERC.DT)
} # Loop through list of files containing weather variables for each pyrome

# Process wthrvars.DT containing data for all Pyromes processed
sapply(wthrvars.DT, class)
collength <- dim(wthrvars.DT)[2] 
convert_to_number <- colnames(wthrvars.DT)[4:collength]
wthrvars.DT[, c(convert_to_number) := lapply(.SD, as.numeric), .SDcols=convert_to_number]
sapply(wthrvars.DT, class)
wthrvars.DT[, ROWID := .I]

# Reformat date column
cols <- c('DATE')
wthrvars.DT[, DATE := lapply(.SD, function(x){gsub('/', '-', x)}), .SDcols = cols]

# Create day/month and year columns
wthrvars.DT.out <- wthrvars.DT[, as.list(str_split_fixed(DATE, '-', 3)), by=DATE] 
setnames(wthrvars.DT.out, old = c('V1', 'V2', 'V3'), new = c('MONTH', 'DAY', 'YEAR'))
setkey(wthrvars.DT,'DATE')
setkey(wthrvars.DT.out,'DATE')
wthrvars.DT <- wthrvars.DT[wthrvars.DT.out, allow=TRUE]
wthrvars.DT[, CORDDATE := paste0(YEAR,'-',MONTH,'-',DAY)]
setorder(wthrvars.DT, ROWID)
collength <- dim(wthrvars.DT)[2] 
convert_to_number <- colnames(wthrvars.DT)[collength - 1] # YEAR column
wthrvars.DT[, c(convert_to_number) := lapply(.SD, as.integer), .SDcols=convert_to_number]

# Remove leap year days - Appears no leap year in data!
wthrvars.DT <- subset(wthrvars.DT, !(DAY == '29' & MONTH == '02'))
wthrvars.DT[, DAYOFYR := as.integer(yday(CORDDATE))]
wthrvars.DT[, DAYOFYR := DAYOFYR][(DAYOFYR > 60) & is.leapyear(YEAR) == TRUE, DAYOFYR := DAYOFYR - 1][]
wthrvars.DT[, c('YEAR') := lapply(.SD, as.numeric), .SDcols='YEAR']

# Extract data variables 
seltvar <- c('PYROME','MONTH','DAY','YEAR','CORDDATE','DAYOFYR','Temp','AvgT','MinT','MaxT','RH','AvRH','MnRH','MxRH','ERC')
wthrvars.DT <- wthrvars.DT[, ..seltvar]
sapply(wthrvars.DT, class)
wthrvars.DT[, PYROME := sprintf("%03d", as.integer(PYROME))] 

##################################################################################################################
## Section: Random tests on data if needed
##################################################################################################################

# Create statistics by day of year and pyrome
wthrvars.daily.stats <- wthrvars.DT[,list(NUMOFYRS = .N,NUMOFDAYS = .N,
                                          AAmeanVL=mean(ERC, na.rm=TRUE),
                                          AAsdVL=sd(ERC, na.rm=TRUE),                                      
                                          AAminVL=min(ERC, na.rm=TRUE),
                                          AAmedianVL = median(ERC, na.rm=TRUE),
                                          AAQ80VL=quantile(ERC, .80, na.rm=TRUE),
                                          AAQ90VL=quantile(ERC, .90, na.rm=TRUE),
                                          AAQ97VL=quantile(ERC, .97, na.rm=TRUE),
                                          AAmaxVL=max(ERC)),by=c('PYROME','DAYOFYR')]
perc80base <- (wthrvars.DT[, .(Q80=quantile(ERC, .80, na.rm=TRUE)),by=c('PYROME')])

BASEQnts <- (wthrvars.DT[, .(Q80=quantile(ERC, .80, na.rm=TRUE),
                                Q90=quantile(ERC, .90, na.rm=TRUE),
                                Q97=quantile(ERC, .97, na.rm=TRUE)),by=c('PYROME')])

##################################################################################################################
## Section: Read FFP5 weather output files for each Pyrome - POST CLIMATE CONDITIONING
##################################################################################################################

# Specify output location
outloc <- paste0('Expmnt_',sprintf("%03d", as.integer(Exprmnt)),'/')

# Check directory
initial_path <- file.path(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/',outloc) )
setwd(initial_path)
wd <- getwd()
print(paste0('Current working dir: ', wd))

# Collect names of file of interest - wlstinv1!500100.txt individual station index files
filesin <- list.files(pattern = 'FFP5Weather', recursive = TRUE, full.names = TRUE)
filenamesin <- basename(list.files(pattern = 'FFP5Weather', recursive = TRUE))
numfiles <- length(filenamesin)

# Get correct files to match to base ERC steam weather files
#pattern <- c('scaled','ZEROPRCIP')
#tifsperrds <- as.matrix(sapply(pattern, function(x) grepl(x, filenamesin)))
#indxs <- which(apply(tifsperrds, 1, function(x) all(x == c(TRUE,TRUE))))
#filenamesin <- filenamesin[indxs]

# Initialize data.table to hold data for all pyromes 
wthrvars.CCD.DT <- data.table(NULL)

# Loop through indexfiles collecting weather variables by pyrome
for (j in 1:length(filesin)) { 
  # Create connection
  file <- filesin[j]
  FLNAMEIN <- filenamesin[j]
  pyromein <- substring(file,3,5)
  
  # Read data to data.table line by line
  ERC.DT <- fread(file, colClasses = "character",
                  sep = ',', header = FALSE, verbose = TRUE)
  
  # Delete meta data and header followed by dummy column 49
  ERC.DT[, 49 := NULL]
  
  # Pass data from text DT to final DT
  Headernames <- c('SIGName','DATE','Time','Temp','AvgT','MinT','MaxT','RH','AvRH','MnRH','MxRH','Rain','RnDr','Wind',
                   'SC','ERC','BI','KBDI','IC','FM1','FM10','F100','F1000','FMH','FMW','SFlag','WDIR','SOW','FFMC',
                   'DMC','DC','ISI','BUI','FWI','DSR','GDir','GSpd','SolR','WFlag','DPT','VPDM','VPDA','GSI','WAZI',
                   'HrRa','LHerb','LWood','FFWI')
  
  oldnames <- names(ERC.DT)
  setnames(ERC.DT, old = oldnames, new = Headernames, skip_absent=TRUE)
  
  # Add Pyrome reference column
  ERC.DT[,PYROME := pyromein]
  
  # Append pyrome data to master data DT
  if (j == 1) {
    wthrvars.CCD.DT <- ERC.DT
  } else {
    # Write top 4 grid locations
    filesappend <- list(wthrvars.CCD.DT, ERC.DT) 
    wthrvars.CCD.DT <- rbindlist(filesappend, use.names=TRUE, fill=TRUE)  
  }
  rm(ERC.DT)
} # Loop through list of files containing weather variables for each pyrome

# Process wthrvars.DT containing data for all Pyromes processed
sapply(wthrvars.CCD.DT, class)
collength <- dim(wthrvars.CCD.DT)[2] 
convert_to_number <- colnames(wthrvars.CCD.DT)[4:collength]
wthrvars.CCD.DT[, c(convert_to_number) := lapply(.SD, as.numeric), .SDcols=convert_to_number]
sapply(wthrvars.CCD.DT, class)
wthrvars.CCD.DT[, ROWID := .I]

# Reformat date column
cols <- c('DATE')
wthrvars.CCD.DT[, DATE := lapply(.SD, function(x){gsub('/', '-', x)}), .SDcols = cols]

# Create day/month and year columns
wthrvars.CCD.DT.out <- wthrvars.CCD.DT[, as.list(str_split_fixed(DATE, '-', 3)), by=DATE] 
setnames(wthrvars.CCD.DT.out, old = c('V1', 'V2', 'V3'), new = c('MONTH', 'DAY', 'YEAR'))
setkey(wthrvars.CCD.DT,'DATE')
setkey(wthrvars.CCD.DT.out,'DATE')
wthrvars.CCD.DT <- wthrvars.CCD.DT[wthrvars.CCD.DT.out, allow=TRUE]
wthrvars.CCD.DT[, CORDDATE := paste0(YEAR,'-',MONTH,'-',DAY)]
setorder(wthrvars.CCD.DT, ROWID)
collength <- dim(wthrvars.CCD.DT)[2] 
convert_to_number <- colnames(wthrvars.CCD.DT)[collength - 1] # YEAR column
wthrvars.CCD.DT[, c(convert_to_number) := lapply(.SD, as.integer), .SDcols=convert_to_number]

# Remove leap year days - Appears no leap year in data!
wthrvars.CCD.DT <- subset(wthrvars.CCD.DT, !(DAY == '29' & MONTH == '02'))
wthrvars.CCD.DT[, DAYOFYR := as.integer(yday(CORDDATE))]
wthrvars.CCD.DT[, DAYOFYR := DAYOFYR][(DAYOFYR > 60) & is.leapyear(YEAR) == TRUE, DAYOFYR := DAYOFYR - 1][]
wthrvars.CCD.DT[, c('YEAR') := lapply(.SD, as.numeric), .SDcols='YEAR']

# Extract data variables 
seltvar <- c('PYROME','MONTH','DAY','YEAR','CORDDATE','DAYOFYR','Temp','AvgT','MinT','MaxT','RH','AvRH','MnRH','MxRH','ERC')
wthrvars.CCD.DT <- wthrvars.CCD.DT[, ..seltvar]
sapply(wthrvars.CCD.DT, class)
wthrvars.CCD.DT[, PYROME := sprintf("%03d", as.integer(PYROME))] 

##################################################################################################################
## Section: Random tests on data if needed
##################################################################################################################

# Create statistics by day of year and pyrome
wthrvars.CCD.daily.stats <- wthrvars.CCD.DT[,list(NUMOFYRS = .N,NUMOFDAYS = .N,
                                          AAmeanVL=mean(ERC, na.rm=TRUE),
                                          AAsdVL=sd(ERC, na.rm=TRUE),                                      
                                          AAminVL=min(ERC, na.rm=TRUE),
                                          AAmedianVL = median(ERC, na.rm=TRUE),
                                          AAQ80VL=quantile(ERC, .80, na.rm=TRUE),
                                          AAQ90VL=quantile(ERC, .90, na.rm=TRUE),
                                          AAQ97VL=quantile(ERC, .97, na.rm=TRUE),
                                          AAmaxVL=max(ERC)),by=c('PYROME','DAYOFYR')]
perc80ccd <- (wthrvars.CCD.DT[, .(Q80=quantile(ERC, .80, na.rm=TRUE)),by=c('PYROME')])

CCDQnts <- (wthrvars.CCD.DT[, .(Q80=quantile(ERC, .80, na.rm=TRUE),
                                  Q90=quantile(ERC, .90, na.rm=TRUE),
                                  Q97=quantile(ERC, .97, na.rm=TRUE)),by=c('PYROME')])

##################################################################################################################
## Section: Evaluate climate impacts - Compare ERC quantiles between base and climate conditioned
##################################################################################################################

# Bind base and climate conditioned ERC steams into 1 data.table
appendlist <- list()
appendlist[[1]] <- wthrvars.DT
appendlist[[2]] <- wthrvars.CCD.DT
wthrvars.bs.ccd <- rbindlist(appendlist, use.names=TRUE, id='ctype')
sapply(wthrvars.bs.ccd, class)

# Bind base and climate quantiles
appendlist <- list()
appendlist[[1]] <- BASEQnts
appendlist[[2]] <- CCDQnts
Qnts.bs.ccd <- rbindlist(appendlist, use.names=TRUE, id='ctype')
sapply(Qnts.bs.ccd, class)

# Total number of pyromes
numpyr <- 128
pyromes <- seq(1,numpyr,1)

# Select subset of pyromes
pyrssub <- c('001','005', '008', '009', '010', '011', '012', '013', '014', 
             '015', '020', '026', '075', '101', '118', '122', '127')

pyrssub <- c('008')

# Climate typs
numtypes <- 2
cltypes <- seq(1,numtypes,1)

# Create results data table
pyrstats <- data.table(NULL)
pyrstats <- data.table(1)[,`:=`(c('PYROME', 'TOTDays', 'NumYR', 'Q00Dys', 
                                  'Q80Dys', 'Q90Dys', 'Q97Dys', 'MNQ80Dys', 
                                  'MNQ90Dys', 'MNQ97Dys', 'PrcQ80DIF', 'PrcQ90DIF', 'PrcQ97DIF' ),0)][,V1:=NULL][!is.na(PYROME)]
pyrstats <- pyrstats[rep(1:.N,numpyr)][,PYROME:=1:.N]
pyrstats <- pyrstats[rep(1:.N,numtypes)][,ctype:=1:.N,by=PYROME]
pyrstats[, PYROME := sprintf("%03d", as.integer(PYROME))] 

# Loop through pyromes to comapre base and climate conditioned ERC streams
cols <- c('ERC')
for (ipyr in pyromes) { 
  PYROMEIN <- sprintf("%03d", as.integer(pyromes[ipyr]))
  idx_one_in_two <- match(PYROMEIN, unlist(pyrssub)) 
  # Select subset of pyromes
  if (!is.na(idx_one_in_two)) { 
  
  # Climate types 1 = base, 2 = conditioned
  for (ctyp in cltypes) { 

    # Percentile values
    #BASEQntsp <- Qnts.bs.ccd[ctype == ctyp & PYROME == PYROMEIN]
    BASEQntsp <- Qnts.bs.ccd[ctype == ctyp & PYROME == PYROMEIN]    
    quants <- as.list(BASEQntsp[,3:5])
    testsum <- sum(BASEQntsp[,3:5]) 
    
    # ERC values
    ERCVals <- wthrvars.bs.ccd[ctype == ctyp & PYROME == PYROMEIN]
    ERCVals[,INDX := .I]
    
    # Test quants !=0
    if (testsum != 0 & nrow(ERCVals) != 0) {
      # Compare ERC to quants
      ERCVals.out <- ERCVals[, .(vecsum = sum(purrr::map2_dbl(.x = .SD, .y = quants, function(x, y) x >= y))), .SDcols = cols, by = INDX]
      statsout <- ERCVals.out[, .N, by = list(vecsum)]
      
      # Calculate number of years
      MinYR <- min(ERCVals[,YEAR])
      MaxYR <- max(ERCVals[,YEAR])
      numYRs <- MaxYR - MinYR
      
      # Accrue results
      pyrstats[PYROME == PYROMEIN & ctype == ctyp, TOTDays := sum(statsout[,N])]      
      pyrstats[PYROME == PYROMEIN & ctype == ctyp, NumYR := numYRs]       
      pyrstats[PYROME == PYROMEIN & ctype == ctyp, Q00Dys := statsout[1,N]] 
      pyrstats[PYROME == PYROMEIN & ctype == ctyp, Q80Dys := statsout[2,N]]       
      pyrstats[PYROME == PYROMEIN & ctype == ctyp, Q90Dys := statsout[3,N]]       
      pyrstats[PYROME == PYROMEIN & ctype == ctyp, Q97Dys := statsout[4,N]]   
      pyrstats[PYROME == PYROMEIN & ctype == ctyp, MNQ80Dys := statsout[2,N]/numYRs]       
      pyrstats[PYROME == PYROMEIN & ctype == ctyp, MNQ90Dys := statsout[3,N]/numYRs]       
      pyrstats[PYROME == PYROMEIN & ctype == ctyp, MNQ97Dys := statsout[4,N]/numYRs]         
    }
    
    # Join results
    #setkey(ERCVals, INDX)
    #setkey(ERCVals.out, INDX)
    
    # Resultant joined table
    #ERCVals <- ERCVals[ERCVals.out, allow=TRUE]
    rm(ERCVals,statsout,BASEQntsp,quants,testsum)

  } # Climate type
    
    # Calculate change in percentils daya counts
    test <- pyrstats[Q00Dys != 0 & PYROME == PYROMEIN]
    setorder(test, PYROME, ctype)
    q80difin <- test[, list(ctype=ctype[-.N],pct.diff=diff(MNQ80Dys) / MNQ80Dys[-.N]),PYROME]  
    q90difin <- test[, list(ctype=ctype[-.N],pct.diff=diff(MNQ90Dys) / MNQ90Dys[-.N]),PYROME]  
    q97difin <- test[, list(ctype=ctype[-.N],pct.diff=diff(MNQ97Dys) / MNQ97Dys[-.N]),PYROME]    
    pyrstats[PYROME == PYROMEIN & ctype == ctyp, PrcQ80DIF := q80difin[,pct.diff]*100.0]       
    pyrstats[PYROME == PYROMEIN & ctype == ctyp, PrcQ90DIF := q90difin[,pct.diff]*100.0]       
    pyrstats[PYROME == PYROMEIN & ctype == ctyp, PrcQ97DIF := q97difin[,pct.diff]*100.0]      
  } # Test if pyrome is in data
} # Pyrome

pyrstats <- pyrstats[Q00Dys != 0]
setorder(pyrstats, PYROME, ctype)
test <- wthrvars.bs.ccd[ctype == 1 & PYROME == '008']
test[ERC >= 68.400]
##################################################################################################################
## Section: Evaluate climate impacts - Compare ERC data between base and climate conditioned
##################################################################################################################

# Set plots per frame
par(mfrow=c(1,1))

# Total number of pyromes
numpyr <- 128
pyromes <- seq(1,numpyr,1)

# Select subset of pyromes
pyrssub <- c('001','005', '008', '009', '010', '011', '012', '013', '014', 
             '015', '020', '026', '075', '101', '118', '122', '127')

# Evaluate labels
Date <- seq(as.Date("1981/1/1"), as.Date("1981/12/31"), "day")
wthrvars.CCD.DT.dt <- wthrvars.CCD.DT[PYROME == '008']
wthrvars.CCD.DT.dt <- unique(wthrvars.CCD.DT.dt, by = 'DAYOFYR')
seltvar <- c('DAYOFYR')
wthrvars.CCD.DT.dt <- wthrvars.CCD.DT.dt[, ..seltvar]
wthrvars.CCD.DT.dt.DF <- setDF(wthrvars.CCD.DT.dt)
wthrvars.CCD.DT.dt.DF <- cbind(wthrvars.CCD.DT.dt.DF,Date)
wthrvars.CCD.DT.dt.DF$MonthDay <- format(wthrvars.CCD.DT.dt.DF$Date, "%d-%b")

# only show every third label... otherwise it's too crowded
wthrvars.CCD.DT.dt.DF$MonthDay[as.numeric(row.names(wthrvars.CCD.DT.dt.DF))%%12!=0] <- ""
labelsin <- wthrvars.CCD.DT.dt.DF$MonthDay
length(labelsin)

# Bind base and climate conditioned ERC steams into 1 data.table
appendlist <- list()
appendlist[[1]] <- wthrvars.daily.stats
appendlist[[2]] <- wthrvars.CCD.daily.stats
wthrvars.daily.stats.all = rbindlist(appendlist, use.names=TRUE, id='type')

# Remove superfluous column
rm.col = c('NUMOFYRS', 'NUMOFDAYS') 
wthrvars.daily.stats.all[,(rm.col):=NULL]

# Loop through pyromes to compare base and climate conditioned ERC streams
for (ipyr in pyromes) { 
  # Evaluate pyrome code
  PYROMEIN <- sprintf("%03d", as.integer(pyromes[ipyr]))
  idx_one_in_two <- match(PYROMEIN, unlist(pyrssub)) 
  # Select subset of pyromes
  if (!is.na(idx_one_in_two)) { 
    # Split data by class for plotting - Wind related functions
    wthrvars.DT.in.all <- wthrvars.daily.stats.all
    wthrvars.DT.pyr <- wthrvars.DT.in.all[PYROME == PYROMEIN] 
    wthrvars.DT.pyr <- wthrvars.DT.pyr[, PYROME := NULL]
    wthrvars.DT.pyr.t <- dcast(melt(wthrvars.DT.pyr, id.vars = 'DAYOFYR'), variable ~ DAYOFYR)
    split_list <- split(wthrvars.DT.pyr.t, by=c('variable'), flatten=FALSE)
    
    # Reference DT's by class name, i.e RES1 etc
    namesin <- names(split_list)
    for(i in seq_along(namesin)){
      assign(namesin[i], split_list[[i]])
    }
    
    # Extract individaul DTs or datframes examples
    wthrvars.DT.pyr_df <- as.data.frame(split_list[[1]]) # as dataframe
    wthrvars.DT.pyr_dt <- as.data.table(split_list[[1]])# as data.table

    # Set data to be plotted - Wind realted functions 
    namesin <- namesin[2:length(namesin)]
    wthrvars.DT.pyr.in <- melt(wthrvars.DT.pyr, measure.vars = patterns('AA'), value.name = 'VL', variable.name='TIME')

    # set output directory
    setwd(paste0('D:/ERCTimeSeriesAnalysis/Output/Expmnt_',sprintf("%03d", as.integer(Exprmnt)),'/plots'))
    wd <- getwd()
    
    # Loop throough list of function names
    for (i in 1:length(namesin)) { 
      # Get data frame, craete table for testing
      test <- namesin[i]
      if(test == 'AAmeanVL') { 
        cat(test,i,"  ")
        wthrvars.DT.pyr.plot <- wthrvars.DT.pyr.in[TIME == test]
        maxin <- max(wthrvars.DT.pyr.plot[, VL]) 
        maxin <- round(maxin+1,0)
        perc80basein <- perc80base[PYROME == PYROMEIN][,Q80]
        perc80ccdin <- perc80ccd[PYROME == PYROMEIN][,Q80]

        # Plot ggplot 
        legend_title <- 'Climate Type'
        p <- ggplot() + geom_line(data=wthrvars.DT.pyr.plot, aes(x=factor(DAYOFYR), y=VL, group = factor(type)
                                                                 , colour = factor(type)), size=1) + 
          theme(legend.position = c(0.9, 0.8)) + 
          scale_color_discrete(labels = c('Current', 'Conditioned'), name = legend_title) +
          scale_x_discrete(labels=labelsin) + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8)) +
          scale_y_continuous(breaks = pretty(wthrvars.DT.pyr.plot$VL, n = 10)) +
          #geom_hline(yintercept = perc80basein) +
          #geom_hline(yintercept = perc80ccdin) +          
          labs(title = paste0('Pyrome ',ipyr,' - Fire Family Plus v5 ERC data'),
             subtitle = paste0(substr(test, 3, 6),' by Day of Year'),
             x = 'Day of Year',
             y = 'ERC') + labs(caption = 'USFS/risQ Climate Conditioning project')
      
        png(width=1000, height=500)
        print(p)
        dev.off()
        ggsave(paste0('Pyrome_',PYROMEIN,'_',namesin[i],'.png'), p, width = 9, height = 9)
        rm(p,test,wthrvars.DT.pyr.in)
      }
    }   
  }
}

# Loop through pyromes to compare base and climate conditioned ERC streams
for (ipyr in pyromes) { 
  # Evaluate pyrome code
  PYROMEIN <- sprintf("%03d", as.integer(pyromes[ipyr]))
  idx_one_in_two <- match(PYROMEIN, unlist(pyrssub)) 
  # Select subset of pyromes
  if (!is.na(idx_one_in_two)) { 
    # Split data by class for plotting - Wind related functions
    wthrvars.DT.in.all <- wthrvars.daily.stats.all
    wthrvars.DT.pyr <- wthrvars.DT.in.all[PYROME == PYROMEIN] 
    wthrvars.DT.pyr <- wthrvars.DT.pyr[, PYROME := NULL]
    wthrvars.DT.pyr.t <- dcast(melt(wthrvars.DT.pyr, id.vars = 'DAYOFYR'), variable ~ DAYOFYR)
    split_list <- split(wthrvars.DT.pyr.t, by=c('variable'), flatten=FALSE)
    
    # Reference DT's by class name, i.e RES1 etc
    namesin <- names(split_list)
    for(i in seq_along(namesin)){
      assign(namesin[i], split_list[[i]])
    }
    
    # Extract individaul DTs or datframes examples
    wthrvars.DT.pyr_df <- as.data.frame(split_list[[1]]) # as dataframe
    wthrvars.DT.pyr_dt <- as.data.table(split_list[[1]])# as data.table
    
    # Set data to be plotted - Wind realted functions 
    namesin <- namesin[2:length(namesin)]
    wthrvars.DT.pyr.in <- melt(wthrvars.DT.pyr, measure.vars = patterns('AA'), value.name = 'VL', variable.name='TIME')
    
    # set output directory
    setwd(paste0('D:/ERCTimeSeriesAnalysis/Output/Expmnt_',sprintf("%03d", as.integer(Exprmnt)),'/plots'))
    wd <- getwd()
    
    # Loop throough list of function names
    #for (i in 1:length(namesin)) { 
      # Get data frame, craete table for testing
      #test <- namesin[i]
      #if(test == 'AAmeanVL') { 
        #cat(test,i,"  ")
        wthrvars.DT.pyr.plot.max <- wthrvars.DT.pyr.in[TIME == namesin[8]]
        maxin <- max(wthrvars.DT.pyr.plot.max[, VL]) 
        maxin <- round(maxin+6,0)

        wthrvars.DT.pyr.plot.mean <- wthrvars.DT.pyr.in[TIME == namesin[1]]
        
        wthrvars.DT.pyr.plot.min <- wthrvars.DT.pyr.in[TIME == namesin[3]]
        
        perc80basein <- perc80base[PYROME == PYROMEIN][,Q80]
        perc80ccdin <- perc80ccd[PYROME == PYROMEIN][,Q80]
        
        # Plot ggplot 
        legend_title <- 'Climate Type'
        p <- ggplot() + geom_line(data=wthrvars.DT.pyr.plot.max, aes(x=factor(DAYOFYR), y=VL, group = factor(type)
                                                                 , colour = factor(type)), size=1) + 
          geom_line(data=wthrvars.DT.pyr.plot.mean, aes(x=factor(DAYOFYR), y=VL, group = factor(type)
                                                   , colour = factor(type)), size=1) + 
          geom_line(data=wthrvars.DT.pyr.plot.min, aes(x=factor(DAYOFYR), y=VL, group = factor(type)
                                                        , colour = factor(type)), size=1) +           
          theme(legend.position = c(0.92, 0.9)) + 
          scale_color_discrete(labels = c('Current climate', 'Conditioned climate')) +
          theme(legend.title=element_blank()) +
          #scale_color_discrete(labels = c('Current', 'Conditioned'), name = legend_title) +          
          theme(legend.text=element_text(size=7)) +
          scale_x_discrete(labels=labelsin) + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8)) +
          #scale_y_continuous(breaks = pretty(wthrvars.DT.pyr.plot.max$VL, n = 10)) +
          scale_y_continuous(limits = c(0, maxin), breaks = seq(0, maxin, by = 10)) +
          geom_hline(yintercept = perc80basein, colour="red") +
          geom_hline(yintercept = perc80ccdin, colour="Turquoise")  +          
          labs(title = paste0('Pyrome ',ipyr,' - Fire Family Plus v5 ERC data'),
               #subtitle = paste0(namesin[i],' by Day of Year'),
               x = 'Day of Year',
               y = 'ERC') + labs(caption = 'USFS/risQ Climate Conditioning Project')
        
        png(width=1000, height=500)
        print(p)
        dev.off()
        ggsave(paste0('Pyrome_',PYROMEIN,'_all.png'), p, width = 9, height = 9)
        rm(p,test,wthrvars.DT.pyr.in)

      #}
    #}   
  }
}

# Read grid point list of all pyromes 
gridfile <- paste0('C:/ERCTimeSeriesAnalysis/Data/Pyromes/','Top4CloseGridPointsAllPyromesCWAtrbs.rds')
WFgrds_shp.DT <- readRDS(gridfile)
setorder(WFgrds_shp.DT, grid_x, grid_y)
WFgrds_shp.DT[PYROME == '010']

##################################################################################################################
## Section: Edit NARR grid location FireFamilPlu plus input weather file .fwx formatted
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

# Load in file containing the climate conditioning change factors and risQ 2-deg grid information
#   - v1: contains preliminary conditioning of temperature and RH
#   - files are on AWS
trgtyr <- 2045
percte <- 0.5
percte <- 'med'
scnio <- 'rcp85'
runspec <- paste0(trgtyr,percte,scnio)
chngfacts <- paste0('change_factors_',trgtyr,'_Q_',percte)
change_factors_file <- readRDS(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/ClimateConditioningInput/Climate_Conditioning_Change_Factors_v1.rds'))
risQ_grid <- change_factors$risQ_grid # the lat-lon and ID #'s for the risQ grid
chngfactsin <- change_factors[chngfacts][[1]]

# source the climate conditioning function - can be found in the wildfire github repo 
#source('C:/GitRepositories/ERCResearch/Scripts/condition_pyrome.R')

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

# Read grid point list of all pyromes 
gridfile <- paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/','Pyromes_2019/Top4CloseGridPointsAllPyromes.rds')
WFgrds.DT <- readRDS(gridfile) 
sapply(WFgrds.DT, class)
WFgrds.DT[,filename := (as.character(filename))]
sapply(WFgrds.DT, class)
WFgrds.DT[, distrank := sprintf("%03d", as.integer(distrank))] 
sapply(WFgrds.DT, class)

# Select subset of pyromes
pyrssub <- c('001','005', '008', '009', '010', '011', '012', '013', '014', 
             '015', '020', '026', '075', '101', '118', '122', '127')

# Specify output location
outloc <- paste0('Expmnt_',sprintf("%03d", as.integer(Exprmnt)),'/')

# Loop through WFgrds_shp.DT to process each NARR griddede weather file by Pyrome
for (iLine in 1:nrow(WFgrds.DT)) { 
  # Evaluate pyrome code
  PYROMEIN <- unlist(WFgrds.DT[iLine,PYROME])
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
    filenamein <- WFgrds.DT[iLine,filename]
    distrankin <- WFgrds.DT[iLine,distrank]
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

    # Climate contidition the variables Temp, MaxT, MinT and RH, MxRH, MnRH
    #grdfile.DT <- condition_pyrome(plat = WFgrds.DT$latitude[iLine], plon = WFgrds.DT$longitude[iLine], scenario = eval(scnio), 
                               #pyrome_data = grdfile.DT, scenario_obj = chngfactsin, grid_pts = risQ_grid)
    
    cc_grdfile.DT <- condition_pyrome(plat = WFgrds.DT$latitude[iLine], plon = WFgrds.DT$longitude[iLine], Year = eval(trgtyr), scenario = eval(scnio),
                                      scenario_obj = eval(percte), pyrome_data = grdfile.DT, change_factors = change_factors_file)    

    # Output file after climate conditioning
    filenamein <- WFgrds.DT[iLine,filename]
    distrankin <- WFgrds.DT[iLine,distrank]
    filedir <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/',outloc,PYROMEIN) 
    fileout <- paste0(distrankin,'_CCD_',runspec,'_',filenamein)
    if (!dir.exists(filedir)){
      dir.create(filedir,recursive = TRUE)
    } else {
      print("Dir already exists!")
    }
    fileout <- paste0(filedir,'/',fileout)
  
    # Check column widths of output data, pad if incorrect
    colconvrt <- names(cc_grdfile.DT)
    cc_grdfile.DT[ , c(colconvrt) := Map(function(x, COLWIDS) {str_pad(x, COLWIDS, side = c('left'), pad = " ")},
                             .SD, COLWIDS), .SDcols = c(colconvrt)]

    # Bind output columns as 1 column
    cc_grdfile.DT[, key_ := do.call(paste, c(.SD, sep = "")), .SDcols = names(cc_grdfile.DT)]
  
    # Use tabs, suppress row names and column names
    data <- cc_grdfile.DT[,key_ ]
    write.table(data, fileout, quote = FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  }
}

##################################################################################################################
## Section: Define miscellaneous functions
##################################################################################################################

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