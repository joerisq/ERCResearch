##################################################################################################################
## Project: Evaluate ERC Time series analysis, as M. Finneys original paper
##
## Script purpose: Process weather inputs to carry out trend analysis
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

# specify working directory; you will need to change this line if run elsewhere
setwd(ROOT_DATA_DIR)

# Check directory
getwd()
wd <- getwd()
print(paste0('Current working dir: ', wd))

##################################################################################################################
## Section: Read Californian RAWS station list and USFS Pyrome boundary file
##################################################################################################################

# Check directory
initial_path <- file.path(paste0(ROOT_DATA_DIR,'/',DATA_DIR))
setwd(initial_path)
wd <- getwd()
print(paste0('Current working dir: ', wd))

# Read Pyrome shapefile
Pyrms_shp <- readOGR(paste0(initial_path,'/','Pyrome_20150605.shp')) 
Pyrms_shp <- spTransform(Pyrms_shp, CRS('+proj=longlat +ellps=WGS84 +datum=WGS84'))
str(Pyrms_shp@data)

# Pyrome of interest
Pyrms.AOI <- Pyrms_shp[Pyrms_shp@data$PYROME == 26,]
centroid.AOI <- gCentroid(Pyrms.AOI,byid=TRUE)
plot(Pyrms.AOI)
points(centroid.AOI,pch=2, col = "Red")

# Read file containing RAWS locations for California as 2016 - Thousand Oaks lacked precision, edited by hand
RAWSFile <- paste0(initial_path,'/','CA_RAWS_Locales_2016.txt')
#RAWS.DT <- fread(RAWSFile, fill=TRUE)

# Read data to data.table line by line
RAWS.DT <- fread(RAWSFile, colClasses = "character",
                 sep = "\n", header = FALSE, verbose = TRUE)

# Column data entered by hand
end_col <- 70 #cumsum(cols)
start_col <- c(1,29,44,53,64) #end_col - cols + 1
end_col <- c(28,43,52,63,70)
start_end <- cbind(start_col, end_col) # matrix of start and end positions

# Function to pass data by column widths
text <- lapply(RAWS.DT, function(x) {
  apply(start_end, 1, function(y) substr(x, y[1], y[2])) 
})

# Pass data from text DT to final DT
RAWS.DT.final <- data.table(text$V1)
Headernames <- c('StationName','State','Elevft','LatDDMMSS','LonDDDM')
setnames(RAWS.DT.final, old = 1:ncol(RAWS.DT.final), new = Headernames)

# Delete embedded header row
RAWS.DT.final <- RAWS.DT.final[-1,] 
RAWS.DT.final[, c(1:5) := lapply(.SD, function(x) trimws(x) )]

# Convert to decinal degrees proper
angle2dec <- function(angle) {
  angle <- as.character(angle)
  x <- do.call(rbind, strsplit(angle, '(?<=..)', perl = TRUE))
  x <- apply(x, 1L, function(y) {
    y <- as.numeric(y)
    y[1] + y[2]/60 + y[3]/3600
  })
  return(x)
}

# Convert Latitudes
setDT(RAWS.DT.final)[, LATLIST := angle2dec(LatDDMMSS)]
setDT(RAWS.DT.final)[, LonLIST := substr(RAWS.DT.final[,LonDDDM], 2, 7)]
setDT(RAWS.DT.final)[, LonLIST := angle2dec(LonLIST)]
setDT(RAWS.DT.final)[, LonLIST := -1.0 * (LonLIST + 100)]

# Function to calculate distance between geo-codes in meters
dt.haversine <- function(lat_from, lon_from, lat_to, lon_to, r = 6378137){
  radians <- pi/180
  lat_to <- lat_to * radians
  lat_from <- lat_from * radians
  lon_to <- lon_to * radians
  lon_from <- lon_from * radians
  dLat <- (lat_to - lat_from)
  dLon <- (lon_to - lon_from)
  a <- (sin(dLat/2)^2) + (cos(lat_from) * cos(lat_to)) * (sin(dLon/2)^2)
  return(2 * atan2(sqrt(a), sqrt(1 - a)) * r)
}

long_orig <- coordinates(centroid.AOI)[1]
lat_orig <- coordinates(centroid.AOI)[2]
RAWS.DT.final[, dist := dt.haversine(lat_orig, long_orig, LATLIST, LonLIST)/1000.0] # convert to km
#RAWS.DT.final[, distrank:=rank(-dist,ties.method="first"),by=StationName]
RAWS.DT.final[, distrank:=rank(dist,ties.method="first")]
colVar <- 'distrank'
setorderv(RAWS.DT.final, colVar)[] 

# Select the top 5
RAWS.DT.final.rnkd <- RAWS.DT.final[distrank < 6]

# prepare the 3 components: coordinates, data, and proj4string
coords <- RAWS.DT.final.rnkd[ , c('LonLIST', 'LATLIST')]   # coordinates
data   <- RAWS.DT.final.rnkd[ , 8:9]          # data
crs    <- CRS("+init=epsg:4326") # proj4string of coords

# make the spatial points data frame object
spdf <- SpatialPointsDataFrame(coords = coords,
                               data = data, 
                               proj4string = crs)
plot(spdf, col = 'Green', add = TRUE)

# Read gridded fire weather grid point shapefile
WFgrds_shp <- readOGR(paste0(initial_path,'/','xy_gridpoints.shp')) 
#WFgrds_shp <- spTransform(WFgrds_shp, CRS('+proj=longlat +ellps=WGS84 +datum=WGS84'))
str(WFgrds_shp@data)
WFgrds_shp.DT <- setDT(WFgrds_shp@data)
long_orig <- coordinates(centroid.AOI)[1]
lat_orig <- coordinates(centroid.AOI)[2]
WFgrds_shp.DT[, dist := dt.haversine(lat_orig, long_orig, latitude, longitude)/1000.0]
#RAWS.DT.final[, distrank:=rank(-dist,ties.method="first"),by=StationName]
WFgrds_shp.DT[, distrank:=rank(dist,ties.method="first")]
colVar <- 'distrank'
setorderv(WFgrds_shp.DT, colVar)[] 

# Select the top 5
WFgrds_shp.DT.final.rnkd <- WFgrds_shp.DT[distrank < 5]

# prepare the 3 components: coordinates, data, and proj4string
coords <- WFgrds_shp.DT.final.rnkd[ , c('longitude', 'latitude')]   # coordinates
data   <- WFgrds_shp.DT.final.rnkd[ , 5:7]          # data
crs    <- CRS("+init=epsg:4326") # proj4string of coords

# make the spatial points data frame object
spdf2 <- SpatialPointsDataFrame(coords = coords,
                                data = data, 
                                proj4string = crs)
plot(spdf2, col = 'Blue', add = TRUE)

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
##            Program end section to create tested accumulated Precipitation output for hurricane event         ##
##                                                                                                              ##
##################################################################################################################


