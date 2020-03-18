##################################################################################################################
## Project: Climate conditioning ERC stream input to fire simualator FSIM
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
library(rgdal)
library(raster) 
library(rgeos)
library(ggrepel)
library(ggspatial)
library(cowplot)
library(rcartocolor)
library(ggmap)
library(maps)
library(mapdata)
library(gridExtra)
library(grid)
library(lubridate)
library(data.table)
library(openair)

# Projection string for meter projection US wide
AEAProj <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-110
+x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m'

# Sheffs preferred output projections
laea_proj4 <- '+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'

# Test for valid GDAL install
gdal_setInstallation()
valid_install <- !is.null(getOption('gdalUtils_gdalPath'))
if(require(raster) && require(rgdal) && valid_install)
{
  print('TRUE')
}
#gdal_setInstallation(ignore.full_scan=FALSE)

# Set raster options
rasterOptions(maxmemory = 5.0e+09)

#if('package:adehabitatHR' %in% search()) detach('package:adehabitatHR', unload=TRUE)

memory.limit()
# set max memory usage 
memory.size(max=5.0e+09)

##################################################################################################################
## Section: Set global directory locations
##################################################################################################################

# Get the data directory from readtext
ROOT_DATA_DIR <- ('F:/ERCTimeSeriesAnalysis')
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
#writeOGR(Pyrms_shp, initial_path, 'PYROME_LL', driver='ESRI Shapefile')
str(Pyrms_shp@data)
plot(Pyrms_shp)

# Pyrome of interest
Pyrms.AOI <- Pyrms_shp[Pyrms_shp@data$PYROME == 26,]
centroid.AOI <- gCentroid(Pyrms.AOI,byid=TRUE)
plot(Pyrms.AOI)
points(centroid.AOI,pch=2, col = 'Red')

# Read file containing RAWS locations for California as 2016 - Thousand Oaks lacked precision, edited by hand
RAWSFile <- paste0(initial_path,'/','CA_RAWS_Locales_2016.txt')
#RAWS.DT <- fread(RAWSFile, fill=TRUE)

# Read data to data.table line by line
RAWS.DT <- fread(RAWSFile, colClasses = 'character',
                 sep = '\n', header = FALSE, verbose = TRUE)

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
#RAWS.DT.final[, distrank:=rank(-dist,ties.method='first'),by=StationName]
RAWS.DT.final[, distrank:=rank(dist,ties.method='first')]
colVar <- 'distrank'
setorderv(RAWS.DT.final, colVar)[] 

# Select the top 5
RAWS.DT.final.rnkd <- RAWS.DT.final[distrank < 6]

# prepare the 3 components: coordinates, data, and proj4string
coords <- RAWS.DT.final.rnkd[ , c('LonLIST', 'LATLIST')]   # coordinates
data   <- RAWS.DT.final.rnkd[ , 8:9]          # data
crs    <- CRS('+init=epsg:4326') # proj4string of coords

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
#RAWS.DT.final[, distrank:=rank(-dist,ties.method='first'),by=StationName]
WFgrds_shp.DT[, distrank:=rank(dist,ties.method='first')]
colVar <- 'distrank'
setorderv(WFgrds_shp.DT, colVar)[] 

# Select the top 5
WFgrds_shp.DT.final.rnkd <- WFgrds_shp.DT[distrank < 5]

# prepare the 3 components: coordinates, data, and proj4string
coords <- WFgrds_shp.DT.final.rnkd[ , c('longitude', 'latitude')]   # coordinates
data   <- WFgrds_shp.DT.final.rnkd[ , 5:7]          # data
crs    <- CRS('+init=epsg:4326') # proj4string of coords

# make the spatial points data frame object
spdf2 <- SpatialPointsDataFrame(coords = coords,
                                data = data, 
                                proj4string = crs)
plot(spdf2, col = 'Blue', add = TRUE)

##################################################################################################################
## Section: Loop through Pyromes acquire station data and write Station locations to data.table
##################################################################################################################

# Read Improved station list 
stationlistfile <- paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/','FAM.NWCG.GOVFAM-files_StationList_2019.rds')
stationlistfile.dt <- readRDS(stationlistfile)
stationlistfile.dt[, FILEPRESNCE := 1]
stationlistfile.dt <- stationlistfile.dt[LONGITUDE != 0][LATITUDE != 0]

# Initialize data.table to hold data for all ptroimes and closest 5 stations
stationlistfile.dt.rnkd.master <- data.table(NULL)

# Check station file exists
for (i in 1:nrow(stationlistfile.dt)) { 
  state <-  stationlistfile.dt[i, STATENAME] 
  filename <- stationlistfile.dt[i, FLNAME] 
  filenametest <- substring(filename, 10, 15)
  stationid <- stationlistfile.dt[i, STATID] 
  if (stationid != filenametest) {
    stationid = filenametest
    cat('Update name ',i, '\n')
    stationlistfile.dt[i, STATID := stationid] 
  }
  statredindx <- stationlistfile.dt[ , STATINDEX] 
  filename <- paste0('wx',stationid,'.fw13')
  filetochck<- paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/','FAM.NWCG.GOVFAM-files-',DATA_YR,'/',state,'/',filename) 
  
  # Check file exists, if not remove record
  if(!file.exists(filetochck)) {
    stationlistfile.dt[i, FILEPRESNCE := 0]
    cat('File not present! ',i, '\n')
  } 
}
stationlistfile.dt <- stationlistfile.dt[FILEPRESNCE == 1]

# Prepare corrected data for output
stationlistfile.rds <- paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/','FAM.NWCG.GOVFAM-files_StationList_',DATA_YR_OUT,'_CORCTD.rds')
stationlistfile.csv <- paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/','FAM.NWCG.GOVFAM-files_StationList_',DATA_YR_OUT,'_CORCTD.csv')
saveRDS(stationlistfile.dt, stationlistfile.rds)
fwrite(stationlistfile.dt, stationlistfile.csv)

# Loop through PYrome region polygons
for (i in 1:nrow(Pyrms_shp)) { 
  # Create/check output directory
  PYRRegin <- Pyrms_shp[i,'PYROME']
  PYRRegincode <- unlist(Pyrms_shp@data[i, 'PYROME'])
  PyromesfOLDER <- paste0('Pyromes_2019','/','StationWeatherData','/',sprintf('%03d', PYRRegincode))
  output_dir_PYR <- file.path(ROOT_DATA_DIR, DATA_DIR, PyromesfOLDER)
  
  if (!dir.exists(output_dir_PYR)){
    dir.create(output_dir_PYR)
  } else {
    print('Dir already exists!')
  }
  
  setwd(output_dir_PYR)
  WD <- getwd()
  if (!is.null(WD)) setwd(WD) 
  print(PYRRegincode)

  # Write pyrome shape out
  writeOGR(obj=PYRRegin, dsn=output_dir_PYR, layer=paste0('Pyrome_',sprintf('%03d', PYRRegincode)), driver='ESRI Shapefile', overwrite=TRUE) 
  
  # Get centroid of Pyrome
  centroid.AOI <- gCentroid(PYRRegin,byid=TRUE)
  
  # Find closest 4 stations
  long_orig <- coordinates(centroid.AOI)[1]
  lat_orig <- coordinates(centroid.AOI)[2]
  stationlistfile.dt[, dist := dt.haversine(lat_orig, long_orig, LATITUDE, LONGITUDE)/1000.0] # convert to km
  stationlistfile.dt[, distrank:=rank(dist,ties.method='first')]
  colVar <- 'distrank'
  setorderv(stationlistfile.dt, colVar)[]
  
  plot(PYRRegin)
  points(centroid.AOI,pch=2, col = 'Red')
  # Select the top 5
  stationlistfile.dt.rnkd <- stationlistfile.dt[distrank < 6]
  stationlistfile.dt.rnkd <- stationlistfile.dt.rnkd[,PYROME := sprintf('%03d', PYRRegincode)]
  
  # prepare the 3 components: coordinates, data, and proj4string
  coords <- stationlistfile.dt.rnkd[ , c('LONGITUDE', 'LATITUDE')]   # coordinates
  data   <- stationlistfile.dt.rnkd[ , 8:11]          # data
  crs    <- CRS('+init=epsg:4326') # proj4string of coords
  
  # make the spatial points data frame object
  spdf <- SpatialPointsDataFrame(coords = coords,
                                 data = data, 
                                 proj4string = crs)
  plot(spdf, col = 'Green', add = TRUE)
  
  # Write top 5 stations
  filesappend <- list(stationlistfile.dt.rnkd.master, stationlistfile.dt.rnkd) 
  #setattr(filesappend, 'names', c('1', '2', '3'))
  stationlistfile.dt.rnkd.master <- rbindlist(filesappend, use.names=TRUE, fill=TRUE)  
  stationlistfile.rds <- paste0(output_dir_PYR,'/','Top5CloseStations_2019.rds')
  stationlistfile.csv <- paste0(output_dir_PYR,'/','Top5CloseStations_2019.csv')
  saveRDS(stationlistfile.dt.rnkd, stationlistfile.rds)
  fwrite(stationlistfile.dt.rnkd, stationlistfile.csv)
  
  # Check if data directory exists
  if (dir.exists(output_dir_PYR)) {
    # Copy station files
    for (j in 1:5) { 
    state <-  stationlistfile.dt.rnkd[distrank == j, STATENAME] 
    filename <- stationlistfile.dt.rnkd[distrank == j, FLNAME] 
    stationid <- stationlistfile.dt.rnkd[distrank == j, STATID] 
    statredindx <- stationlistfile.dt.rnkd[distrank == j, STATINDEX] 
    filetocopyDels <- paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/','FAM.NWCG.GOVFAM-files-',DATA_YR,'/',state,'/',filename) 
    fileindx <- sprintf('%03d', j)
    filename <- paste0('wx',stationid,'.fw13')
    filetocopyWthr <- paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/','FAM.NWCG.GOVFAM-files-',DATA_YR,'/',state,'/',filename) 
    filetorenme <- paste0(output_dir_PYR,'/',paste0(fileindx,'_wx',stationid,'.fw13'))     
      if (isTRUE(dir.exists(output_dir_PYR))) {
        file.copy(filetocopyDels, output_dir_PYR)
        file.copy(filetocopyWthr, output_dir_PYR)  
        file.rename(from = filename, to = paste0(fileindx,'_',filename))
        if(!file.exists(paste0(fileindx,'_',filename))) {
          stationlistfile.dt.rnkd.master <- stationlistfile.dt.rnkd.master[STATINDEX  == statredindx, PYROME := '999'] 
          print('File not present!')
        }
      } else {
        print('Dir already exists!')
      }  
    }
  } # Check if data directory exists
}

# Write concatenated data.table of all pyromes and closest stations
stationlistfile <- paste0('D:/ERCTimeSeriesAnalysis/Data/Pyromes_2019/','Top5CloseStationsAllPyromes_',DATA_YR_OUT,'.rds')
saveRDS(stationlistfile.dt.rnkd.master, stationlistfile)
stationlistfile.dt.rnkd.master[PYROME == '999']

##################################################################################################################
## Section: Loop through Pyromes and set up inset plot maps
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

centroidS <- gCentroid(Pyrms_shp,byid=TRUE)
#plot(Pyrms_shp)
#points(centroidS,pch=2, col = 'Red')
leg_centers <- SpatialPointsDataFrame(gCentroid(Pyrms_shp, byid=TRUE), 
                                      Pyrms_shp@data, match.ID=FALSE)
# Create data frame from centroids
leg_centers.cntds <- data.frame(leg_centers@data, x = coordinates(leg_centers)[,1], y = coordinates(leg_centers)[,2])
leg_centers.cntds <- subset(leg_centers.cntds, select=c(PYROME, x, y))

# Prepare US wide data for plotting
Pyrms_shp.count <- nrow(Pyrms_shp@data)
Pyrms_shp@data$id <- 1:Pyrms_shp.count
Pyrms_shp.fort <- fortify(Pyrms_shp, region='id')

# pyrome as sf object
Pyrms_shp.sf <- st_as_sf(as(Pyrms_shp, 'SpatialLinesDataFrame'))
Pyrms_shp.sf_2163 = st_transform(Pyrms_shp.sf, crs = 2163)

# Read concatenated data.table of all pyromes and closest stations
stationlistfile <- paste0('F:/ERCTimeSeriesAnalysis/Data/Pyromes_2019/','Top5CloseStationsAllPyromes_',DATA_YR_OUT,'.rds')
stationlistfile.dt.rnkd.master <- readRDS(stationlistfile)

# Read concatenated data.table of all pyromes and closest NARR grids 
gridfile <- paste0('F:/ERCTimeSeriesAnalysis/Data/Pyromes_2019/','Top4CloseGridPointsAllPyromesCWAtrbs.rds')
WFgrds_shp.DT.master <- readRDS(gridfile) 

# Read concatenated data.table of all pyromes station eind data combined
cmbndstdatfile <- paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/CombinedStationDataByPyrome03/AllPyromesOrigWindStatDat_DDUPED.rds')
stationdata.DT <- readRDS(cmbndstdatfile) 
sapply(stationdata.DT, class) 
stationdata.DT[, DATE := parse_date_time(paste0(DATE,' 13:00:00'), 'Ymd HMS')]

# Loop through Pyrome region polygons
for (i in 1:nrow(Pyrms_shp)) { 
  # Create/check output directory
  PYRRegin <- Pyrms_shp[i,'PYROME']
  PYRRegincode <- unlist(Pyrms_shp@data[i, 'PYROME'])
  
  # Get region of interest
  PYRRegin_bb = st_as_sfc(st_bbox(PYRRegin))

  # Create output directory if required
  outputDIR <- file.path(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/Plots/PyromeMaps/'))
  if (!dir.exists(outputDIR)) {dir.create(outputDIR)}
  setwd(outputDIR) 

  # Get centroid of Pyrome
  centroid.AOI <- gCentroid(PYRRegin,byid=TRUE)
  
  # Find closest 5 stations - the 3 components: coordinates, data, and proj4string
  stationlistfile.dt.rnkd <- stationlistfile.dt.rnkd.master[PYROME == sprintf('%03d', PYRRegincode)]
  coordinates(stationlistfile.dt.rnkd)<-~LONGITUDE+LATITUDE
  class(stationlistfile.dt.rnkd)
  proj4string(stationlistfile.dt.rnkd) <- CRS('+init=epsg:4326')
  stationlistfile.dt.rnkd<-spTransform(stationlistfile.dt.rnkd, CRS('+init=epsg:2163'))
  names(stationlistfile.dt.rnkd)[names(stationlistfile.dt.rnkd)=='LONGITUDE'] <- 'x'
  names(stationlistfile.dt.rnkd)[names(stationlistfile.dt.rnkd)=='LATITUDE'] <- 'y'
  spdf.st <- stationlistfile.dt.rnkd
  
  # Find closest 4 NARR Grids - the 3 components: coordinates, data, and proj4string
  WFgrds_shp.DT <- WFgrds_shp.DT.master[PYROME == sprintf('%03d', PYRRegincode)]  
  coordinates(WFgrds_shp.DT) <- ~longitude+latitude
  class(WFgrds_shp.DT)
  proj4string(WFgrds_shp.DT) <- CRS('+init=epsg:4326')
  WFgrds_shp.DT<-spTransform(WFgrds_shp.DT, CRS('+init=epsg:2163'))
  names(WFgrds_shp.DT)[names(WFgrds_shp.DT)=='longitude'] <- 'x'
  names(WFgrds_shp.DT)[names(WFgrds_shp.DT)=='latitude'] <- 'y'  
  spdf.gd <- WFgrds_shp.DT
  
  # Combine points
  pntsmerg <- list(coordinates(stationlistfile.dt.rnkd),coordinates(WFgrds_shp.DT))
  pntsmergmgd <- do.call('rbind', pntsmerg) 
  multipoints <- st_multipoint(pntsmergmgd, dim = 'XY')
  points <- st_cast(st_geometry(multipoints), 'POINT') 
  ptsbbin <- st_bbox(points)
  ptsbbinply <- st_as_sfc(ptsbbin)
  st_crs(ptsbbinply) = 2163 
  aoibbin <- st_bbox(Pyrms_shp.sf_2163[i,'PYROME']) 
  aoibbinply <- st_as_sfc(aoibbin)
  cmbnd <- st_union(ptsbbinply, aoibbinply, by_feature = FALSE)
  
  bbincmbnd <- st_bbox(cmbnd)
  xminbb = bbincmbnd[[1]]
  xmaxbb = bbincmbnd[[3]]
  yminbb = bbincmbnd[[2]]
  ymaxbb = bbincmbnd[[4]]
  
  xliminpts = c(xminbb, xmaxbb) 
  yliminpts = c(yminbb, ymaxbb)    

  # Convert to Dfs
  spdf.st.sf <- st_as_sf(spdf.st)
  spdf.gd.sf <- st_as_sf(spdf.gd)
  centroid.sf <- st_as_sf(centroid.AOI)
  
  # Prepare subset of pyromes as sf object
  PYRRegin_bb = st_as_sfc(st_bbox(Pyrms_shp.sf_2163[i,'PYROME']))
  
  Pdata <- stationlistfile.dt.rnkd.master[PYROME == sprintf('%03d', PYRRegincode)]
  statenm <- Pdata[1, STATENAME]
  statenm <- toupper(statenm)
  title <- paste0('Location of Pyrome ',sprintf('%03d', PYRRegincode),' in State of ',statenm)
  #title <- paste0(title,' (Blue RAW stations, red NARR Grids)')  
  par(mar=c(0,0,0,0))
  theme_set(theme_bw())
  
  bbin <- st_bbox(PYRRegin_bb)
  xminbb = bbin[[1]]
  xmaxbb = bbin[[3]]
  yminbb = bbin[[2]]
  ymaxbb = bbin[[4]]
  
  xlimin = c(xminbb, xmaxbb) 
  ylimin = c(yminbb, ymaxbb)  
  
  xdiff <- unlist(min(abs(xminbb),abs(xmaxbb)))
  xdiffabs <- abs(xminbb - abs(xdiff)*-sign(xminbb))
  
  ydiff <- min(abs(yminbb),abs(ymaxbb))
  ydiffabs <- abs(yminbb - abs(ydiff)*-sign(yminbb))
  
  aspectratio <- (ydiffabs / xdiffabs) 
  if(aspectratio < 0.3) {aspectratio <- 0.3}
  
  #bbox_new <- st_bbox(PYRRegin_bb) # current bounding box
  #xrange <- bbox_new$xmax - bbox_new$xmin # range of x values
  #yrange <- bbox_new$ymax - bbox_new$ymin # range of y values
  #bbox_new[1] <- bbox_new[1] - (0.5 * xrange[[1]]) # xmin - left
  #bbox_new[3] <- bbox_new[3] + (0.5 * xrange[[1]]) # xmax - right
  #bbox_new[2] <- bbox_new[2] - (0.5 * yrange[[1]]) # ymin - bottom
  #bbox_new[4] <- bbox_new[4] + (0.5 * yrange[[1]]) # ymax - top
  
  #xminbb = bbox_new[[1]]
  #xmaxbb = bbox_new[[3]]
  #yminbb = bbox_new[[2]]
  #ymaxbb = bbox_new[[4]]
  
  #xlimin = c(xminbb, xmaxbb) 
  #ylimin = c(yminbb, ymaxbb) 
  
  #bbox_new <- bbox_new %>%  # take the bounding box ...
  #  st_as_sfc() # ... and make it a sf polygon 
  
  black.bold.italic.text <- element_text(face = "bold.italic", color = "black", size = 8)
  black.bold.italic.6.text <- element_text(face = "bold.italic", color = "black", size = 6)
  locale <- ggplot() + 
    geom_sf(data = Pyrms_shp.sf_2163[i,'PYROME'], aes(fill = PYROME)) + guides(fill=FALSE) +
    geom_sf(data = spdf.gd.sf, size = 3, shape = 24, colour = 'red', fill = NA) +
    geom_sf(data = spdf.st.sf, size = 2, shape = 22, colour = 'blue', fill = 'blue') +
    geom_sf(data = centroid.sf, size = 4, shape = 13, colour = 'green', fill = 'green') + 
    coord_sf(crs = st_crs(2163), xlim = xliminpts, ylim = yliminpts) +
    theme(panel.background = element_blank()) + 
    theme(plot.title = element_text(color='black', size=8, face="bold.italic")) +
    ggtitle(title) + xlab('Longitude') + ylab('Latitude') + 
    theme(axis.text.x = black.bold.italic.6.text) +
    theme(axis.text.y = black.bold.italic.6.text) + theme(axis.title = black.bold.italic.text)

  zoomin <- ggplot() + 
    geom_sf(data = Pyrms_shp.sf_2163, color = 'darkgrey', fill = 'grey92') + 
    geom_sf(data = PYRRegin_bb, fill = NA, color = 'red', size = 0.5) +
    coord_sf(crs = st_crs(2163), xlim = c(-2150000, 2500000), ylim = c(-2150000, 
                                                                       730000),  datum = NA) +
    theme(
      panel.ontop = TRUE,   ## Note: this is to make the panel grid visible in this example
      panel.grid = element_blank(), 
      line = element_blank(), 
      rect = element_blank(), 
      text = element_blank(), 
      plot.background = element_rect(fill = 'white')) 
  
  bbin <- st_bbox(Pyrms_shp.sf_2163)
  xminbb = bbin[[1]]
  xmaxbb = bbin[[3]]
  yminbb = bbin[[2]]
  ymaxbb = bbin[[4]]
  
  xdiff <- unlist(min(abs(xminbb),abs(xmaxbb)))
  xdiffabs <- abs(xminbb - abs(xdiff)*-sign(xminbb))
  
  ydiff <- min(abs(yminbb),abs(ymaxbb))
  ydiffabs <- abs(yminbb - abs(ydiff)*-sign(yminbb))
  
  aspectratio <- (ydiffabs / xdiffabs) 
  if(aspectratio < 0.3) {aspectratio <- 0.3}

  #gg_inset_map <- ggdraw() +
   # draw_plot(locale) +
    #draw_plot(zoomin, x = 0.00, y = 0.00, width = 0.2, height = 0.2)
  
  #gg_inset_map <- locale +
    #annotation_custom(grob = zoomin, xmin = 0, xmax = 300000,
                      #ymin = 0, ymax = 300000)
  
  # Create wind rose plots by month
  stationdata.DT.PYR <- stationdata.DT[PYROME == sprintf('%03d', PYRRegincode)]
  stationdata.DT.PYR[, WS := WS * 0.45]
  setDF(stationdata.DT.PYR)

  wrsepyr <- windRose(stationdata.DT.PYR, 
           wd = 'WD',
           ws = 'WS', type = 'MONTH', angle = 45, auto.text= FALSE, annotate = FALSE,
           width = 0.2, grid.line = 10, 
           paddle = F, bias.corr=F)

  # Image export details
  filename <- paste0('OverviewMap_Pyr',sprintf('%03d', PYRRegincode),'.pdf')  
  pdf(filename, width = 6.8, height = 4.7)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  print(zoomin, vp = vplayout(1, 2))  # key is to define vplayout
  print(locale, vp = vplayout(1:2, 1))
  print(wrsepyr, vp = vplayout(2, 2))  
  dev.off()

  # Clean memory
  rm(zoomin,locale,gg_inset_map,zoomin)
}

##################################################################################################################
## Section: Loop through Pyromes acquire gridded NARR data and write grid pointa to data.table
##################################################################################################################

# Grid point attributes read earlier
gridattribs.dt <- WFgrds_shp.DT

# Initialize data.table to hold data for all pyromes and closest 5 stations
gridattribs.dt.rnkd.master <- data.table(NULL)

# Load selenium to control chrom efor file downloads
require(RSelenium)
library(RSelenium)
library(pingr)

# Start server and base website 4444:4450
sapply(9515:9515, function(p) setNames(list(pingr::ping_port('localhost', p)), p))
available.versions <- binman::list_versions('chromedriver') 
latest.version <- available.versions$win32[length(available.versions)][[1]][1]
#latest.version <- '79.0.3945.36'
rD <- rsDriver(port = 9515L, browser=c('chrome'),version = 'latest', chromever=latest.version)
remDr <- rD[['client']]

remDr$navigate('https://wrcc.dri.edu/fpa/gridded/?form_lat=&form_lon=')

# Loop through PYrome region polygons
for (i in 1:nrow(Pyrms_shp)) {
  # Create/check output directory
  PYRRegin <- Pyrms_shp[i,'PYROME']
  PYRRegincode <- unlist(Pyrms_shp@data[i, 'PYROME'])
  PyromesfOLDER <- paste0('Pyromes','/','GriddedWeatherData','/',sprintf('%03d', PYRRegincode))
  output_dir_PYR <- file.path(ROOT_DATA_DIR, DATA_DIR, PyromesfOLDER)
  
  if (!dir.exists(output_dir_PYR)){
    dir.create(output_dir_PYR)
  } else {
    print('Dir already exists!')
  }
  
  setwd(output_dir_PYR)
  WD <- getwd()
  if (!is.null(WD)) setwd(WD) 
  print(PYRRegincode)
  
  # Write pyrome shape out
  writeOGR(obj=PYRRegin, dsn=output_dir_PYR, layer=paste0('Pyrome_',sprintf('%03d', PYRRegincode)), driver='ESRI Shapefile', overwrite=TRUE) 
  
  # Get centroid of Pyrome
  centroid.AOI <- gCentroid(PYRRegin,byid=TRUE)
  
  # Find closest 4 grid points
  long_orig <- coordinates(centroid.AOI)[1]
  lat_orig <- coordinates(centroid.AOI)[2]
  gridattribs.dt[, dist := dt.haversine(lat_orig, long_orig, latitude, longitude)/1000.0] # convert to km
  gridattribs.dt[, distrank:=rank(dist,ties.method='first')]
  colVar <- 'distrank'
  setorderv(gridattribs.dt, colVar)[]
  
  plot(PYRRegin)
  points(centroid.AOI,pch=2, col = 'Red')
  # Select the top 5
  gridattribs.dt.rnkd <- gridattribs.dt[distrank < 5]
  gridattribs.dt.rnkd <- gridattribs.dt.rnkd[,PYROME := sprintf('%03d', PYRRegincode)]
  
  # prepare the 3 components: coordinates, data, and proj4string
  coords <- gridattribs.dt.rnkd[ , c('longitude', 'latitude')]   # coordinates
  data   <- gridattribs.dt.rnkd[ , 3:8]          # data
  crs    <- CRS('+init=epsg:4326') # proj4string of coords
  
  # make the spatial points data frame object
  spdf <- SpatialPointsDataFrame(coords = coords,
                                 data = data, 
                                 proj4string = crs)
  plot(spdf, col = 'Green', add = TRUE)
  
  # Write top 4 grid locations
  filesappend <- list(gridattribs.dt.rnkd.master, gridattribs.dt.rnkd) 
  gridattribs.dt.rnkd.master <- rbindlist(filesappend, use.names=TRUE, fill=TRUE)  
  gridfile.rds <- paste0(output_dir_PYR,'/','Top5CloseGridPoints_2019.rds')
  gridfile.csv <- paste0(output_dir_PYR,'/','Top5CloseGridPoints_2019.csv')
  saveRDS(gridattribs.dt.rnkd, gridfile.rds)
  fwrite(gridattribs.dt.rnkd, gridfile.csv)
  
  # Check if data directory exists
  if (dir.exists(output_dir_PYR)) {
    # Set working directory
    setwd(output_dir_PYR)
    WD <- getwd()
    if (!is.null(WD)) setwd(WD) 
    print(PYRRegincode)
    
    # Cdownload grid files
    for (j in 1:4) { 
      # Get grid x and y count to enter for download
      grid_x <- as.character(gridattribs.dt.rnkd[distrank == j, grid_x]) 
      grid_y <- as.character(gridattribs.dt.rnkd[distrank == j, grid_y]) 
      
      wxbutton <-remDr$findElement(using = 'xpath','//*[(@id = "id_x_field")]')
      wxbutton$clearElement()
      wxbutton$sendKeysToElement(list(grid_x))
      wxbutton$clickElement()
      
      wxbutton <-remDr$findElement(using = 'xpath','//*[(@id = "id_y_field")]')
      wxbutton$clearElement()
      wxbutton$sendKeysToElement(list(grid_y,key = 'enter'))
      wxbutton$clickElement()
      Sys.sleep(10)         
      
      # Files to copy and rename
      fileindx <- sprintf('%03d', j)
      grid_x_in <- sprintf('%03d', as.integer(grid_x))
      grid_y_in <- sprintf('%03d', as.integer(grid_y))
      filenamein <- paste0(grid_x_in,grid_y_in,'_gridded.fwx')
      filenameout <- paste0(fileindx,'_',grid_x_in,grid_y_in,'_gridded.fwx')
      filetocopyWthr <- paste0('C:/Users/joe66/downloads','/',filenamein) 
      if (isTRUE(dir.exists(output_dir_PYR))) {
        file.copy(filetocopyWthr, output_dir_PYR)
        file.rename(from = filenamein, to = filenameout)  

      } else {
        print('Dir already exists!')
      }  
    }
  } # Check if data directory exists 
}

# stop the selenium server
rD[['server']]$stop()
remDr$close()
rm(rD)
# If user forgets to stop server it will be garbage collected.
#gc()

# Write concatenated data.table of all pyromes and closest stations
gridfile.rds <- paste0('D:/ERCTimeSeriesAnalysis/Data/Pyromes_2019/','Top4CloseGridPointsAllPyromes.rds')
saveRDS(gridattribs.dt.rnkd.master, gridfile.rds)

##################################################################################################################
## Section: Loop through Pyromes acquire closest station to grid points
##################################################################################################################

# Read grid point list of all pyromes 
gridfile <- paste0('D:/ERCTimeSeriesAnalysis/Data/Pyromes_2019/','Top4CloseGridPointsAllPyromes.rds')
WFgrds_shp.DT <- readRDS(gridfile)
WFgrds_shp.DT[, grid_x := sprintf('%03d', as.numeric(as.character(grid_x)))]
WFgrds_shp.DT[, grid_y := sprintf('%03d', as.numeric(as.character(grid_y)))]
WFgrds_shp.DT[, GRIDID := paste0(grid_x,grid_y)]

# Read station list of all pyromes 
stationlistfile <- paste0('D:/ERCTimeSeriesAnalysis/Data/Pyromes_2019/','Top5CloseStationsAllPyromes_2019.rds')
stationlistfile.dt <- readRDS(stationlistfile)

# Extract data variables from station list
stationlistfile.dt.in <- stationlistfile.dt[LONGITUDE < 0.0]
stationlistfile.dt.in <- stationlistfile.dt.in[!is.na(LATITUDE)]
seltvar <- c('STATID','LONGITUDE','LATITUDE')
stationlistfile.dt.in <- stationlistfile.dt.in[, ..seltvar]
sapply(stationlistfile.dt.in, class)

colVar <- 'LATITUDE'
setorderv(stationlistfile.dt, colVar)[]

# Extract data variables from grid point list
seltvar <- c('GRIDID','longitude','latitude')
WFgrds_shp.DT.in <- WFgrds_shp.DT[, ..seltvar]
sapply(WFgrds_shp.DT.in, class)

WFgrds_shp.DT.in <- WFgrds_shp.DT.in[rep(seq(1, nrow(WFgrds_shp.DT.in)), nrow(stationlistfile.dt))]
colVar <- 'GRIDID'
setorderv(WFgrds_shp.DT.in, colVar)[]
stationlistfile.dt.in <- stationlistfile.dt.in[rep(seq(1, nrow(stationlistfile.dt.in)), nrow(WFgrds_shp.DT))]

# create Origin-destination matrix
orig <- WFgrds_shp.DT.in[]
dest <- stationlistfile.dt.in[]
odmatrix <- bind_cols(orig,dest)
colnames(odmatrix) <- c('GRIDID', 'long_orig', 'lat_orig', 'STATID', 'long_dest', 'lat_dest')
odmatrix[, dist := dt.haversine(lat_orig, long_orig, lat_dest, long_dest)/1000.0]

#odmatrix <- odmatrix[origi_id == '068199']
odmatrix.frst <- setorder(odmatrix, GRIDID, dist)[, indx := seq_len(.N), by = GRIDID][indx <= 1]
seltvar <- c('GRIDID', 'STATID')
odmatrix.frst <- odmatrix.frst[, ..seltvar]

# Append closest station to grid data
setkey(WFgrds_shp.DT,'GRIDID')
setkey(odmatrix.frst,'GRIDID')
WFgrds_shp.DT <- WFgrds_shp.DT[odmatrix.frst, allow=TRUE]
colVar <- c('PYROME', 'distrank')
setorderv(WFgrds_shp.DT, colVar)[]
sapply(WFgrds_shp.DT, class)

# Write concatenated data.table of all pyromes and closest grid locales
#gridfile <- paste0('D:/ERCTimeSeriesAnalysis/Data/Pyromes_2019/','Top4CloseGridPointsAllPyromes.rds')
#saveRDS(WFgrds_shp.DT, gridfile)

##################################################################################################################
## Section: Calculate and append elevation attributes to grid locations
##################################################################################################################

# Read grid point list of all pyromes 
gridfile <- paste0('D:/ERCTimeSeriesAnalysis/Data/Pyromes_2019/','Top4CloseGridPointsAllPyromes.rds')
WFgrds_shp.DT <- readRDS(gridfile) 
WFgrds_shp.DT.DF <- as.data.frame(WFgrds_shp.DT)

# Read 1km Elevation raster covering North America, data in meters
DEMProj <- '+proj=laea +lat_0=-100 +lon_0=6370997 +x_0=45 +y_0=0 +datum=WGS84 +units=m +no_defs'
laea_proj4 <- '+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
DEMProj <- laea_proj4
filename <- paste0(ROOT_DATA_DIR,'/',DATA_DIR,'/Elevation_GRID/','USIKMElev.tif')
elevdat <- raster::raster(filename)
plot(elevdat) 

# Correct elevation data for extreme value coding
correlev <- function(x) {
  z <- ifelse(x >= 32768, x - 65536 , x)
}
elevdat.CORR <- calc(elevdat, fun=correlev)
crs(elevdat.CORR) <- DEMProj
plot(elevdat.CORR)

# Create derivatives - slope terrain(x, opt=c('slope'), unit='degrees', neighbors=8)
slope_perc <- terrain(elevdat.CORR, 'slope', neighbors=8)
plot(slope_perc)
slope_perc <- round(slope_perc * 100.0, 2) 
title('1KM DEM or Derivatives')

# Reclassify slope percentages
m <- c(-0.1, 25.999, 1, 25.999, 40.999, 2, 40.999, 55.999, 3, 55.999, 75.999, 4, 75.999, 100, 5)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
slope_perc.rc <- reclassify(slope_perc, rclmat)
plot(slope_perc.rc)

# Save reprojected raster
#if (require(rgdal)) {
#  rf <- writeRaster(landform, filename=paste0('C:/ERCTimeSeriesAnalysis/Data/landform1.tif'), 
#                   format='GTiff', overwrite=TRUE)
#}

# Create derivatives - aspect
aspect <- terrain(elevdat.CORR, 'aspect', unit='degrees', neighbors=8)
plot(aspect) 
title('1KM DEM or Derivatives')

# Reclassify aspect percentages
m <- c(0, 1, 0, 1, 45, 1, 45, 90, 2, 90, 135, 3, 135, 180, 4, 
          180, 225, 5, 225, 270, 6, 270, 315, 7, 315, 360, 8)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
aspect.rc <- reclassify(aspect, rclmat)
plot(aspect.rc)

# Write grid locations shape out
#writeOGR(obj=WFgrds.spdf.trnsfrmd, dsn='C:/ERCTimeSeriesAnalysis/Data', 
#         layer='gridlocations', driver='ESRI Shapefile', overwrite=TRUE) 

# TPI for different neighborhood size:
# first step: define customized function
tpiw <- function(x, w=5) {
  m <- matrix(1/(w^2-1), nc=w, nr=w)
  m[ceiling(0.5 * length(m))] <- 0
  f <- focal(x, m) 
  x - f
}
# second step: apply the function
tpi5 <- tpiw(elevdat.CORR, w=5)
tpi5n <- setMinMax(tpi5)
col <- rainbow(20)
plot(tpi5n, col=col, main='Topographic Position Index in Death Valley')
quantile(tpi5[],na.rm=T)
hist(tpi5[])

# Get the standard deviation of the TPI
SD <- sd(tpi5[],na.rm=T)

# Make landform classes
#Morphologic class De Reu et al. 2013;  Weiss (2001)
landform <- reclassify(tpi5, matrix(c(-Inf, -SD, 1,
                                      -SD, -SD/2, 2,
                                      -SD/2, 0, 3,
                                      0, SD/2, 4,
                                      SD/2, SD, 5,
                                      SD, Inf, 6),
                                    ncol = 3, byrow = T),
                       right = T)

# Turn it into categorical raster
landform <- as.factor(landform) 
rat <- levels(landform)[[1]]
#rat[['landform']] <- c('Valley', 'Lower Slope', 
#                       'Flat Area','Middle Slope', 
#                       'Upper Slope', 'Ridge')

rat[['landform']] <- c('Valley bottom or flat ', 'Valley bottom or flat ', 
                       'Valley bottom or flat ', 'Valley bottom or flat ', 
                       'Mid-slope', 'Ridge or peak')
levels(landform) <- rat 
# Plot the classification
x11(12,12)

levelplot(landform, col.regions = rev(brewer.pal(6,'RdYlBu')),
          labels = rat$landform,
          main = 'Landform Classification',
          colorkey=list(labels=list(at=1:6, labels=rat[['landform']])))

# Create spatial points data frame from grid point locations data.table
options(stringsAsFactors=FALSE)
coords <- WFgrds_shp.DT.DF[ , c('longitude', 'latitude')] # coordinates
data   <- WFgrds_shp.DT.DF[ , 8:10] # data
crs    <- CRS('+init=epsg:4326') # proj4string of coords

# make the spatial points data frame object
WFgrds.spdf <- SpatialPointsDataFrame(coords = coords,
                                      data = data, 
                                      proj4string = crs)

# Project data points to DEM projection
WFgrds.spdf.trnsfrmd <- spTransform(WFgrds.spdf,DEMProj)
#plot(WFgrds.spdf.trnsfrmd, col = 'Blue',  add = TRUE, cex=0.5)
summary(WFgrds.spdf.trnsfrmd)
compareCRS(crs(WFgrds.spdf.trnsfrmd), crs(elevdat.CORR))

# Extract raster attributes
WFgrds.spdf.trnsfrmd.elv <- extract(elevdat.CORR, WFgrds.spdf.trnsfrmd)
WFgrds.spdf.trnsfrmd.elv[ is.na(WFgrds.spdf.trnsfrmd.elv) ] <- 0
WFgrds.spdf.trnsfrmd.elv <- as.integer(WFgrds.spdf.trnsfrmd.elv * 3.2808399) # Convert to feet
WFgrds.spdf.trnsfrmd.slp <- extract(slope_perc.rc, WFgrds.spdf.trnsfrmd)
WFgrds.spdf.trnsfrmd.slp[ is.na(WFgrds.spdf.trnsfrmd.slp) ] <- 1
WFgrds.spdf.trnsfrmd.apt <- extract(aspect.rc, WFgrds.spdf.trnsfrmd)
WFgrds.spdf.trnsfrmd.apt[ is.na(WFgrds.spdf.trnsfrmd.apt) ] <- 0
WFgrds.spdf.trnsfrmd.ldf <- extract(landform, WFgrds.spdf.trnsfrmd)
WFgrds.spdf.trnsfrmd.ldf[ is.na(WFgrds.spdf.trnsfrmd.ldf) ] <- 1
results.df <- cbind.data.frame(WFgrds.spdf.trnsfrmd.elv,WFgrds.spdf.trnsfrmd.slp,WFgrds.spdf.trnsfrmd.apt,WFgrds.spdf.trnsfrmd.ldf)
colnames(results.df) <- c('ELEVTN_ft', 'SLOPECLS', 'ASPCTCLS', 'SITE')
results.df <- cbind.data.frame(WFgrds_shp.DT.DF,results.df)
WFgrds_shp.DT.OUT <- setDT(results.df)

# Write concatenated data.table of all pyromes and closest stations
gridfile <- paste0('D:/ERCTimeSeriesAnalysis/Data/Pyromes_2019/','Top4CloseGridPointsAllPyromesCWAtrbs.rds')
saveRDS(WFgrds_shp.DT.OUT, gridfile)

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
##            Program end section to Process weather inputs to carry out trend analysis                         ##
##                                                                                                              ##
##################################################################################################################


