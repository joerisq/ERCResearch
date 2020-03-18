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

# Check directory
initial_path <- file.path(paste0(ROOT_DATA_DIR,'/',DATA_DIR))
setwd(initial_path)
wd <- getwd()
print(paste0('Current working dir: ', wd))

# Read Pyrome shapefile
Pyrms_shp <- readOGR(paste0(initial_path,'/','Pyrome_20150605.shp')) 
Pyrms_shp <- spTransform(Pyrms_shp, laea_proj4)
Pyrms.AOI <- Pyrms_shp[Pyrms_shp@data$PYROME == 26,]

# Read original bp data
DEMProj <- '+proj=laea +lat_0=-100 +lon_0=6370997 +x_0=45 +y_0=0 +datum=WGS84 +units=m +no_defs'
laea_proj4 <- '+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
DEMProj <- laea_proj4

filename <- paste0('D:/WildfireOriginalData/OriginalBPData/iBP.tif')
origbp <- raster::raster(filename)
plot(origbp) 
crs(origbp) <- laea_proj4

filename <- paste0('C:/WildfireFSIMRelated/Outputs_Rerun26CorCMDXFile/PY026r1_BurnProb.asc')
rerunbp <- raster::raster(filename)
plot(rerunbp) 
crs(rerunbp) <- laea_proj4
rerunbp[rerunbp <= 0] <- NA

rerunext <- extent(rerunbp)

origbp.cr <- crop(origbp, rerunext)
plot(origbp.cr) 

origbp.cr.msk <- mask(origbp.cr, rerunbp)
plot(origbp.cr.msk) 

diff <- origbp.cr.msk/rerunbp
diff[diff <= 0] <- NA
diff[diff >= 50] <- NA

plot.new()
par(mfrow=c(2,2))
plot(rerunbp, main='P26 Rerun by risQ, 0314/2020') 
plot(origbp.cr, main='P26 Original BP Cropped') 
plot(origbp.cr.msk, main='P26 Original BP Cropped/masked') 
plot(diff, main='P26 difference (Orig/new)') 


DSMhist<-hist(diff,
              breaks=10,
              main="",
              col="wheat3",  # changes bin color
              xlab= "")



r <- raster(volcano)
plot(diff, col=topo.colors(100), legend=FALSE, axes=FALSE)
r.range <- c(minValue(diff), maxValue(diff))


plot(diff, legend.only=TRUE, col=topo.colors(100),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(at=seq(r.range[1], r.range[2], 25),
                    labels=seq(r.range[1], r.range[2], 25), 
                    cex.axis=0.6),
     legend.args=list(text='Elevation (m)', side=4, font=2, line=2.5, cex=0.8))


plot(diff, legend.only=TRUE, col=topo.colors(100), legend.width=1, legend.shrink=0.75,
     smallplot=c(0,.09, .3,.75)); par(mar = par("mar"))

filename <- paste0('C:/WildfireFSIMRelated/Outputs_Rerun26CorCMDXFile/P26Diff.tif')

# Save reprojected raster
if (require(rgdal)) {
  rf <- writeRaster(diff, filename, format='GTiff', overwrite=TRUE)
}


plot.new()
par(mfrow=c(2,4))
filename <- paste0('C:/WildfireFSIMRelated/Pyrome_026_02/PY26_INPUTS/PY_026_2014b.lcp')
p26lcp <- raster::stack(filename) 
crs(p26lcp) <- laea_proj4
plot(p26lcp[[1]]) 
plot(p26lcp[[2]]) 
plot(p26lcp[[3]]) 
plot(p26lcp[[4]]) 
plot(p26lcp[[5]]) 
plot(p26lcp[[6]]) 
plot(p26lcp[[7]]) 
plot(p26lcp[[8]]) 


filename <- paste0('C:/WildfireFSIMRelated/Outputs_Rerun26CorCMDXFile/P26FuelsMod.tif')
# Save reprojected raster
if (require(rgdal)) {
  rf <- writeRaster(p26lcp[[4]], filename, format='GTiff', overwrite=TRUE)
}






