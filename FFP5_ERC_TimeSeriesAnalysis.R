##################################################################################################################
## Project: Investigate creating time series extropolations of ERC data sets 
##
## Script purpose: Input ERC steam created with Fire Family Plus Version 5.0
##
## Date: 8th November 2019
## Author: JR
##################################################################################################################

# Load Libraries
library(fUnitRoots)
library(lmtest)
library(forecast)
library(FitAR)
library(data.table)
library(stringr)
library(lubridate)
library(TTR)
library(tseries)

setDTthreads(6)

##################################################################################################################
## Section: Set global directory locations
##################################################################################################################

# Get the data directory from readtext
ROOT_DATA_DIR <- ('C:/ERCTimeSeriesAnalysis')
DATA_DIR <- ('Data') 
OUTPUT_DIR <- ('Output')
#DATA_TYPE <- ('Station')
DATA_TYPE <- ('Station')

# specify working directory; you will need to change this line if run elsewhere
setwd(ROOT_DATA_DIR)

# Check directory
getwd()
wd <- getwd()
print(paste0('Current working dir: ', wd))

##################################################################################################################
## Section: Read FireFamilyPlus 5 ERC stream DAILY output data to carry out trend analytics
##################################################################################################################

# Check directory
initial_path <- file.path(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/',DATA_TYPE))
setwd(initial_path)
wd <- getwd()
print(paste0('Current working dir: ', wd))

# Read input ERC data
#ERCFile <- paste0(initial_path,'/','DailyERC_FFPV5.txt') 
ERCFile <- paste0(initial_path,'/','DailyERC_FFPV5_FAMWEBData.txt')

# Read data to data.table line by line
ERC.DT <- fread(ERCFile, colClasses = "character",
                sep = "\n", header = FALSE, verbose = TRUE)

# Delete meta data and header 
ERC.DT <- ERC.DT[-(1:31), ]

# Column data entered by hand
if (DATA_TYPE == 'Station') {
  end_col <- 22 #19 #cumsum(cols)
  start_col <- c(1,11,17) #c(1,9,14) #end_col - cols + 1
  end_col <- c(10,16,22) #c(8,13,19)
} else {
  end_col <- 22 #cumsum(cols)
  start_col <- c(1,11,17) #end_col - cols + 1
  end_col <- c(10,16,22)
}
start_end <- cbind(start_col, end_col) # matrix of start and end positions

# Function to pass data by column widths
text <- lapply(ERC.DT, function(x) {
  apply(start_end, 1, function(y) substr(x, y[1], y[2])) 
})

# Pass data from text DT to final DT
ERC.DT.final <- data.table(text$V1)
Headernames <- c('DATE','Time','ERC')
setnames(ERC.DT.final, old = 1:ncol(ERC.DT.final), new = Headernames)

# Delete embedded header row
ERC.DT.final[, c(1:3) := lapply(.SD, function(x) trimws(x))]
sapply(ERC.DT.final, class)
collength <- dim(ERC.DT.final)[2] 
Cols <- colnames(ERC.DT.final)
convert_to_number <- colnames(ERC.DT.final)[collength:collength]
ERC.DT.final[, c(convert_to_number) := lapply(.SD, as.numeric), .SDcols=convert_to_number]
sapply(ERC.DT.final, class)

# convert to time series _ Edit start date by hand for current test version
data_ts <- ts(ERC.DT.final[, ERC],start = c(1963, 1),frequency=365.25)
attributes(data_ts)

# decompose into time series components
timeseriescomponents <- decompose(data_ts, type=c("additive"))
plot(timeseriescomponents)

# detemine stationarity of data
urkpssTest(data_ts, type = c("tau"), lags = c("short"),use.lag = NULL, doplot = TRUE)
tsstationary<-diff(data_ts, differences=1)
plot(tsstationary)
acf(data_ts,lag.max=34)

# remove seasonality
timeseriesseasonallyadjusted <- data_ts - timeseriescomponents$seasonal
plot(timeseriesseasonallyadjusted)
tsstationary <- diff(timeseriesseasonallyadjusted, differences=1)
plot(tsstationary)
par(mfrow=c(2,1))
acf(tsstationary, lag.max=34) 
pacf(tsstationary, lag.max=34)

# fit the model
fitARIMA<-arima(timeseriesseasonallyadjusted, order=c(1,1,1),seasonal = list(order = c(1,0,0), period = 12),method="ML")
fitARIMA

# significance of coefficients
coeftest(fitARIMA)
par(mfrow=c(1,1))
acf(fitARIMA$residuals)

# residual diagnostics
boxresult<-LjungBoxTest (fitARIMA$residuals,k=2,StartLag=1) 
par(mfrow=c(2,1))
plot(boxresult[,3],main="Ljung-Box Q Test", ylab="P-values", xlab="Lag")
qqnorm(fitARIMA$residuals)
qqline(fitARIMA$residuals)

# Improve fit
fitARIMA <- auto.arima(timeseriesseasonallyadjusted, trace=TRUE)

# forcast future values
par(mfrow=c(1,1))
predict(fitARIMA,n.ahead = 365)
futurVal <- forecast(fitARIMA,h=365*10, level=c(99.5))
plot(futurVal)

# Test simulation for 100 years
plot(arima.sim(list(fitARIMA, ma=.9),  n=365*100), col=4, ylab="x", main= 'Test simulation 365 days')

##################################################################################################################
## Section: Read FireFamilyPlus 5 ERC stream HOURLY output data to investigate parameters
##################################################################################################################

# Check directory
initial_path <- file.path(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/',DATA_TYPE))
setwd(initial_path)
wd <- getwd()
print(paste0('Current working dir: ', wd))

# Read input ERC data
ERCFile <- paste0(initial_path,'/','HourlyDailyERC_FFPV5.txt.txt')

# Read data to data.table line by line
ERC.DT <- fread(ERCFile, colClasses = "character",
                sep = "\n", header = FALSE, verbose = TRUE)

# Delete meta data and header 
ERC.DT <- ERC.DT[-(1:31), ]

# Column data entered by hand
end_col <- 22 #19 #cumsum(cols)
start_col <- c(1,11,17) #c(1,9,14) #end_col - cols + 1
end_col <- c(10,16,22) #c(8,13,19)
start_end <- cbind(start_col, end_col) # matrix of start and end positions

# Function to pass data by column widths
text <- lapply(ERC.DT, function(x) {
  apply(start_end, 1, function(y) substr(x, y[1], y[2])) 
})

# Pass data from text DT to final DT
ERC.DT.final <- data.table(text$V1)
Headernames <- c('DATE','Time','ERC')
setnames(ERC.DT.final, old = 1:ncol(ERC.DT.final), new = Headernames)

# Delete embedded header row
ERC.DT.final[, c(1:3) := lapply(.SD, function(x) trimws(x))]
sapply(ERC.DT.final, class)
collength <- dim(ERC.DT.final)[2] 
Cols <- colnames(ERC.DT.final)
convert_to_number <- colnames(ERC.DT.final)[collength:collength]
ERC.DT.final[, c(convert_to_number) := lapply(.SD, as.numeric), .SDcols=convert_to_number]
sapply(ERC.DT.final, class)

# Reformat date column
cols <- c('DATE')
#out_cols = paste('CORRTD', cols, sep = ".")
#ERC.DT.final[, c(out_cols) := lapply(.SD, function(x){gsub('/', '-', x)}), .SDcols = cols]
ERC.DT.final[, DATE := lapply(.SD, function(x){gsub('/', '-', x)}), .SDcols = cols]

# Create day/month and year columns
ERC.DT.final.out <- ERC.DT.final[, as.list(str_split_fixed(DATE, '-', 3)), by=DATE] 
setnames(ERC.DT.final.out, old = c('V1', 'V2', 'V3'), new = c('MONTH', 'DAY', 'YEAR'))
setkey(ERC.DT.final,'DATE')
setkey(ERC.DT.final.out,'DATE')
ERC.DT.final <- ERC.DT.final[ERC.DT.final.out, allow=TRUE]
ERC.DT.final[, CORDDATE := paste0(YEAR,'-',MONTH,'-',DAY)]
# Remove leap year days
ERC.DT.final <- ERC.DT.final[DAY != '29']
ERC.DT.final[, DAYOFYR := yday(CORDDATE)]
ERC.DT.final[, DAYOFYR := DAYOFYR][(DAYOFYR >= 28) & is.leapyear(as.integer(YEAR)) == TRUE, DAYOFYR := DAYOFYR - 1][]

# Create statistics by day of ERC
ERC.daily.stats <- ERC.DT.final[,list(NUMOFYRS = .N,
                                           mean=mean(ERC),
                                           sd=sd(ERC),                                      
                                           min=min(ERC),
                                           median = median(ERC),
                                           Q80=quantile(ERC, .80, na.rm=TRUE),
                                           Q90=quantile(ERC, .90, na.rm=TRUE),
                                           Q97=quantile(ERC, .97, na.rm=TRUE),
                                           max=max(ERC)),
                                     by=DAYOFYR]

# Plot summary statistics
tsDataMean<- ts(ERC.daily.stats[, mean])
tsDataMin<- ts(ERC.daily.stats[, min])
tsDataMax<- ts(ERC.daily.stats[, max])
x <- ERC.daily.stats[, DAYOFYR]
matplot(x, cbind(tsDataMean,tsDataMin,tsDataMax),type="l",col=c('red','green','blue'),lty=c(1,1,1))
tsm <- cbind(tsDataMean, tsDataMin, tsDataMax)
plot(tsm)

# Clean ERC Values - Fill in missing timeslots
lstr <- nrow(ERC.DT.final)
#strtdate <- paste0(ERC.DT.final[1,CORDDATE],' ',ERC.DT.final[1,Time])
#enddate <- paste0(ERC.DT.final[lstr,CORDDATE],' ',ERC.DT.final[lstr,Time])
#ts <- seq.POSIXt(as.POSIXlt(strtdate), as.POSIXlt(enddate), by="day", tz="MDT")
#ts <- format.POSIXct(ts,'%m/%d/%y %H:%M')
ts <- seq.POSIXt(as.POSIXlt(ERC.DT.final[1,CORDDATE]), as.POSIXlt(ERC.DT.final[lstr,CORDDATE]), by="day", tz="MDT")
ts <- format.POSIXct(ts,'%Y-%m-%d')
length(ts)
ts.dt <- setDT(as.data.frame(ts))
setkey(ERC.DT.final,'CORDDATE')
setkey(ts.dt,'ts')
ERC.DT.final <- ERC.DT.final[ts, allow=TRUE]
for (j in names(ERC.DT.final)) set(ERC.DT.final,which(is.na(ERC.DT.final[[j]])),j,0)
ERC.DT.final[, DAYOFYR := yday(CORDDATE)]
ERC.DT.final <- ERC.DT.final[DAYOFYR != 29]
ERC.DT.final[, DAYOFYR := DAYOFYR][(DAYOFYR >= 28) & is.leapyear(as.integer(YEAR)) == TRUE, DAYOFYR := DAYOFYR - 1][]

# convert to time series
tsdata <- ts(ERC.DT.final[, ERC], start = c(as.integer(ERC.DT.final[1,YEAR]), 1),frequency=365)
attributes(tsdata)
plot(tsdata) 

# Test for stationarity
#ERC.DT.final[, ERC := ERC][(ERC <= 0) , ERC := 0.01][]
#adickeytest <- adf.test(diff(log(tsdata)), alternative="stationary", k=0)
#cycle(tsdata)
#plot(aggregate(tsdata,FUN=mean))
#boxplot(tsdata~cycle(tsdata))

# decompose into time series components
timeseriescomponents <- decompose(tsdata, type=c("additive"))
plot(timeseriescomponents)

##################################################################################################################

# Test fit annual mean trend
library(MASS)
y_ERC = (ERC.daily.stats[, mean])
x_DYR = (ERC.daily.stats[, DAYOFYR])
w_SD = ERC.daily.stats[, sd]
Y <- data.frame(X = x_DYR, Y = y_ERC+0.1)
which(complete.cases(w_SD) == FALSE)
w_SD[is.na(w_SD)] <- mean(w_SD, na.rm = TRUE)
w_SD[w_SD < 1.0] <- 1.0

model <- lm(y_ERC ~ poly(x_DYR,26))
Y_pred <- predict(model, data = new_X_data)
plot(y_ERC)

matplot(x_DYR, cbind(y_ERC,Y_pred),type="l",col=c('red','green'),lty=c(1,1))
summary(model)

##################################################################################################################
## Section: Define miscellaneous functions
##################################################################################################################

# Test if year is leap or not
is.leapyear <- function(year){
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
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
##            Program end section to create time series extropolations of ERC data sets                         ##
##                                                                                                              ##
##################################################################################################################
