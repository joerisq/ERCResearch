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
library(plyr)
library(dplyr)
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

# Temporary raster location
raster::rasterOptions(tmpdir = wd)

# Specify output folder suffix
Exprmnt <- 1

##################################################################################################################
## Section: Read Station data
##################################################################################################################

# Re-read concatenated data.table of all pyromes station data combined
cmbndstdatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/AllPyromesOrigStatDat.rds')
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
stationdata.DT[, TIME := TIME ][TIME == '  :45:00', TIME := '13:45:00'] 
stationdata.DT[, DTME := paste0(DATE,' ',TIME)]
stationdata.DT[, DATETIME := parse_date_time(DTME, 'Ymd HMS')]
stationdata.DT[, DTME := NULL]
setorder(stationdata.DT, PYROME)

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
stationdata.DT[,MONAB := month.abb[MONTH]]
stationdata.DT[,YRQRT := paste0('Q',ceiling(as.numeric(MONTH) / 3))]
setnames(stationdata.DT,'10MWSPD','WS')
setnames(stationdata.DT,'WNDDIR','WD')

# 0 wind direction implies unknown
stationdata.DT[, WD := WD ][WD == 0, WD := NA]
stationdata.DT <- stationdata.DT[!is.na(WD)]
stationdata.DT <- stationdata.DT[!is.na(WS)]

# Write concatenated data.table of all pyromes station data combined
cmbndstdatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/AllPyromesOrigWindStatDat_DDUPED.rds')
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
cmbndstdatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/AllPyromesOrigWindStatDat_DDUPED.rds')
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
  filedir <- file.path(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/'))
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

##################################################################################################################
## Section: Create random windspeeds and directions for each pyrome, save to file
##################################################################################################################

# Total number of pyromes
numpyr <- 128
pyromes <- seq(1,numpyr,1)

# Select subset of pyromes
pyrssub <- sprintf('%03d',seq(1:128))

# Month sequence
mnthsin <- seq(1,12,1)

# Wind speed bin sequence
wsin <- seq(1,8,1)

# Define wind speed column locations within data
str_col <- c(2,8,14,20,26,32,38,42) 
end_col <- c(7,13,19,25,31,37,43,49)
start_end <- cbind(str_col, end_col) 
start_end.dt <- setDT(melt(data = as.data.table(start_end),
           id.vars = c('str_col','end_col')))

# Loop through pyromes 
for (ipyr in pyromes) { 

  # Evaluate pyrome code
  PYROMEIN <- sprintf("%03d", as.integer(pyromes[ipyr]))
  idx_one_in_two <- match(PYROMEIN, unlist(pyrssub))
  
  # Select subset of pyromes
  if (!is.na(idx_one_in_two)) { 

  # Create winds direction and speed CDF to generate random draws of both
  fileindx <- sprintf("%03d", ipyr)
  filedir <- file.path(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/'))
  fileout <- paste0('OrigWindStatDat',fileindx,'.rds')
  #fileout <- 'AllPyromesOrigWindStatDat_DDUPED.rds'
  fileout <- paste0(filedir,'/',fileout)  
  wnds.DT <- readRDS(fileout) 
  sapply(wnds.DT, class)
  
  # Enforce to be numeric
  convrttype <- c(as.character(seq(45,360,45)))
  mtimes <- rep(c('numeric'),each=length(convrttype))
  mtimes <- setNames(mtimes, convrttype)
  set_colclass(wnds.DT,mtimes)
  
  # Fix errors in input files
  is.na(wnds.DT)<-sapply(wnds.DT, is.nan)
  wnds.DT[is.na(wnds.DT)] <- 0

  # Enforce to be numeric
  convrttype <- c(as.character(seq(45,360,45)))
  mtimes <- rep(c('numeric'),each=length(convrttype))
  mtimes <- setNames(mtimes, convrttype)
  set_colclass(wnds.DT,mtimes)
  
  # Wind directions by Months as PDF 
  WD.dist.mnth <- wnds.DT[ , lapply(.SD, sum, na.rm=TRUE), .SDcols=c(4:11), by = 'MONTH']   
  
  # Winds Speed by month and direction as PDF
  WS.dist.mnthd <- dcast(wnds.DT, MONTH  ~ WSBIN, fun=sum, value.var=c(convrttype))
  
    # Loop through months for current pyrome
    num.samples <- 15000
    rndm.dirs.dt.pyr <- data.table(NULL)
    for (imnth in mnthsin) {
    
      # Get direction data for current month
      collength <- dim(WD.dist.mnth)[2]
      imnth.WD.dist <- as.vector(unlist(WD.dist.mnth[imnth,.SD,.SDcols=c(2:collength)]))
      imnth.WD.dist <- ave(imnth.WD.dist,cumsum(imnth.WD.dist))
      
      # Smooth wind direction distribution
      if (length(which(imnth.WD.dist != 0)) < 4 & length(which(imnth.WD.dist != 0)) > 0) {
        x <- imnth.WD.dist
        seq <- seq(1, length(x), 1)
        a_mean <- sapply(seq, function(i) {mean(x[i:(i-1)])})
        y <- rev(x)
        b_mean <- sapply(seq, function(i) {mean(y[i:(i-1)])})   
        test <- (a_mean + rev(b_mean))/2.0
        test <- test/sum(test)
        testc <- (test*0.25 + x*0.75)/2.0
        testc <- testc/sum(testc)  
        imnth.WD.dist <- testc
        cat('Pyrome', ipyr, 'Month', imnth, '\n')
      }

      #names.wd.vec <- names(imnth.WD.dist)
      if (sum(imnth.WD.dist) == 0) imnth.WD.dist[] <- 1/8
      wd.vec <- imnth.WD.dist
      names(wd.vec) <- 1:8
      wd.samples  <- numeric(num.samples)
      
      # Clean renormalize
      wd.vec <- wd.vec/sum(wd.vec)

        # Comput random draws for current month
        set.seed(ipyr*imnth)
        for (i in seq_len(num.samples) ) {
          wd.samples[i] <- discrete.inv.transform.sample(wd.vec)
        }
      
        # Set plots per frame
        par(mfrow=c(1,2))
        barplot(wd.vec, main='True Probability Mass Function')
        barplot(table(wd.samples), main='Empirical Probability Mass Function')   
    
        # Save current months simulated directions for current pyrome
        rndm.dirs.dt <- data.table(NULL)
        rndm.dirs.dt <- data.table(WD = wd.samples)[, INDX := .I]   
        rndm.dirs.dt[,PYROME := PYROMEIN]
        rndm.dirs.dt[,MONTH := imnth]  
        setkey(rndm.dirs.dt, INDX)
      
        # Get wind speed data for current month
        collength <- dim(WS.dist.mnthd)[2]
        imnth.WS.dist <- as.vector(unlist(WS.dist.mnthd[imnth,.SD,.SDcols=c(1:collength)]))

        # Loop through each of the speed bins for current month
        rndm.winds.dt <- data.table(NULL)
        for (iws in wsin) { 
        
          # Get wind speed data for each wind direction in current month
          strt <- start_end.dt[iws,str_col]
          end <- start_end.dt[iws,end_col]
          imnth.ws.dist.iws <- as.vector(unlist(imnth.WS.dist[c(strt:end)]))
          
          # Smooth wind direction distribution
          if (length(which(imnth.ws.dist.iws != 0)) < 3 & length(which(imnth.ws.dist.iws != 0)) > 0) {
            x <- imnth.ws.dist.iws
            seq <- seq(1, length(x), 1)
            a_mean <- sapply(seq, function(i) {mean(x[i:(i-1)])})
            y <- rev(x)
            b_mean <- sapply(seq, function(i) {mean(y[i:(i-1)])})   
            test <- (a_mean + rev(b_mean))/2.0
            test <- test/sum(test)
            testc <- (test*0.25 + x*0.75)/2.0
            testc <- testc/sum(testc)            
            imnth.ws.dist.iws <- testc
            cat('Pyrome', ipyr, 'Month', imnth, 'Wind speed', iws, '\n')
          }

          #names.wd.vec <- names(imnth.WD.dist)
          ws.vec <- imnth.ws.dist.iws
          names(ws.vec) <- 1:6
          ws.samples  <- numeric(num.samples)
          
          # Clean renormalize
          if (sum(ws.vec) == 0) ws.vec[] <- expm1(seq(6,1,-1))
          ws.vec <- ws.vec/sum(ws.vec)

            # Compute random draws for current month
            set.seed(ipyr*imnth*iws)
            for (i in seq_len(num.samples) ) {
              ws.samples[i] <- discrete.inv.transform.sample(ws.vec)
            }
            
            # Save current simulated wind speeds for each directions for current pyrome
            colnme <- paste0('WD',sprintf("%03d", as.integer(45*iws)),'ALLWSBNS')
            rndm.winds.dt <- data.table(DUMMY = ws.samples)[, INDX := .I]  
            setnames(rndm.winds.dt,'DUMMY',colnme)
            
            setkey(rndm.winds.dt, INDX)
            rndm.dirs.dt <- rndm.dirs.dt[rndm.winds.dt,] 
            
            # Append all months for current pyrome 
            if (imnth == 1) {
              rndm.dirs.dt.pyr <- rndm.dirs.dt
            } else {
              filesappend <- list(rndm.dirs.dt.pyr, rndm.dirs.dt) 
              rndm.dirs.dt.pyr <- rbindlist(filesappend, use.names=TRUE, fill=TRUE)  
            }
            
            # Set plots per frame
            par(mfrow=c(1,2))
            barplot(ws.vec, main='True Probability Mass Function')
            barplot(table(ws.samples), main='Empirical Probability Mass Function')

        } # Wind speed bins loop
    } # Month loop
    
    # Enforce to be numeric
    convrttype <- c('WD')
    mtimes <- rep(c('integer'),each=length(convrttype))
    mtimes <- setNames(mtimes, convrttype)
    set_colclass(rndm.dirs.dt.pyr,mtimes)
    # Correct for correct column locations
    rndm.dirs.dt.pyr[,WDCRTD := WD + 4]
    # Create column containing simulated wind speed based on simulated direction in WD
    rndm.dirs.dt.pyr[, WS := .SD[[.BY[[1]]]], by=WDCRTD]
    
    # Create random wind speed within wind bin
    rndm.dirs.dt.pyr[, RNDM := rbeta(.N, 1, 1)]
    rndm.dirs.dt.pyr[, WSDMY := (((WS - 1) + RNDM)* 5.0)]   
    rndm.dirs.dt.pyr[, WSDMY := round_any(WSDMY, 0.5, ceiling)]       
    rndm.dirs.dt.pyr[, WSDMY := round_any(WSDMY, 1.0, floor)]  

    # Create random directionwithin wind bin
    rndm.dirs.dt.pyr[, WDDMY := (((WD - 1) + RNDM)* 45.0)]   
    rndm.dirs.dt.pyr[, WDDMY := round_any(WDDMY, 0.5, ceiling)]       
    rndm.dirs.dt.pyr[, WDDMY := round_any(WDDMY, 1.0, floor)]     
    
    # Drop NA records
    drop_rows_all_na <- function(x, pct=0.01) x[!rowSums(is.na(x)) >= ncol(x)*pct,]
    rndm.dirs.dt.pyr <- drop_rows_all_na(rndm.dirs.dt.pyr,0.001)

    # Write concatenated data.table
    filedir <- file.path(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/'))
    fileout <- paste0(filedir,'/','RndmWindsByPyrome',PYROMEIN,'.rds')
    saveRDS(rndm.dirs.dt.pyr,fileout)
    
  } # Select pyromes
} # All pyromes   

setDF(rndm.dirs.dt.pyr)
wrsepyr <-   windRose(rndm.dirs.dt.pyr, ws = "WSDMY", wd = "WDDMY", ws2 = NA, wd2 = NA, 
                      ws.int = 5, angle = 45, type = "MONTH", bias.corr = TRUE, cols
                      = "default", grid.line = 10, width = 1, seg = 2, auto.text 
                      = TRUE, breaks = 5, offset = 10, normalise = FALSE, max.freq = 
                        NULL, paddle = TRUE, key.header = NULL, key.footer = "(m/s)", 
                      key.position = "bottom", key = TRUE, dig.lab = 5, statistic = 
                        "prop.count", pollutant = NULL, annotate = TRUE, angle.scale = 
                        315, border = T,key.columns = 1,pch=seq(8,8),cex=2,
                      par.settings=list(fontsize=list(text=8)))

##################################################################################################################
## Section: Read NARR Gridded data
##################################################################################################################

# Re-read combined input data
cmbndgddatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/AllPyromesOrigNARRGridDat.rds')
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
cmbndgddatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/AllPyromesOrigWindNARRGridDat_DDUPED.rds')
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
cmbndgddatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/AllPyromesOrigWindNARRGridDat_DDUPED.rds')
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
  filedir <- file.path(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/'))
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

# Turn vector into data.table splitiing entry
#nVector <- c("20 37", "38 23", "39 48", "45 76", "65 44", "86 95 80")
#data.table(v1 = nVector)[, index := .I][, list(unlist(strsplit(v1, " "))), by = index]

# note: this inefficient implementation is for pedagogical purposes
# in general, consider using the rmultinom() function
discrete.inv.transform.sample <- function( p.vec ) { 
  U <- runif(1)
  if(U <= p.vec[1]){
    return(1) 
  } 
  for(state in 2:length(p.vec)) { 
    if(sum(p.vec[1:(state-1)]) < U && U <= sum(p.vec[1:state]) ) { 
      return(state) 
    } 
  }
} 

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






