library(sf)
library(sp)
library(ggplot2)
library(lubridate)
library(gridExtra)
library(MASS)
library(dplyr)
library(mgcv)
library(Metrics)
library(tictoc)
library(stats)

ROOT_DATA_DIR <- "~/Documents/risQ/FinalERCProjectData/ERCTimeSeriesAnalysis"
DATA_DIR <- "Data"
OUTPUT_DIR <- "Output"

#data <- readRDS(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/P008_ERCStats.rds'))

data_base <- readRDS(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/USFSERCCmpCWWethr_P008BaseERCFSIMInput_92-19.rds'))
data_cc <- readRDS(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/USFSERCCmpCWWethr_P008CCDWERCFSIMInput_92-19.rds'))

ERC <- data_base %>%
  dplyr::group_by(SRCE,DAYOFYR) %>%
  dplyr::summarise(max_ERC = max(ERCYR,na.rm = TRUE),
            min_ERC = min(ERCYR,na.rm = TRUE),
            mean_ERC = mean(ERCYR,na.rm = TRUE)) %>%
  dplyr::ungroup()

ERC_stats <- data_base %>%
  dplyr::group_by(SRCE) %>%
  dplyr::summarise(base_Q80 = quantile(ERCYR,probs = 0.8),
                   base_Q90 = quantile(ERCYR,probs = 0.9),
                   base_Q97 = quantile(ERCYR,probs = 0.97)) %>%
  dplyr::ungroup()

risq_ERC <- data_base %>%
  filter(SRCE == 'RISQ') %>%
  dplyr::select(SRCE,DAYOFYR,ERCYR, DATE) %>%
  mutate(SRCE = 'Historical') %>%
  rbind(data_cc %>% 
          dplyr::select(SRCE,DAYOFYR,ERCYR,DATE) %>% 
          filter(SRCE == 'RISQ') %>%
          mutate(SRCE = '2050')) %>%
  dplyr::rename(Climate = SRCE)

risQ_ERC_mean <- risq_ERC %>%
  dplyr::group_by(Climate,DAYOFYR) %>%
  dplyr::summarise(mean_ERC = mean(ERCYR,na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  mutate(DATE = as.Date("2020-01-01","%Y-%m-%d") + DAYOFYR-1)

## plot the base risQ ERC and the climate conditioned ERC, with historical 80th percentile plotted
Q80 <- unlist(ERC_stats %>% filter(SRCE=='RISQ') %>% dplyr::select(Q80))

#ggplot() +
#  geom_hline(yintercept=Q80, linetype="dashed", color = "black") +
#  geom_line(data=risQ_ERC_mean,aes(x=DAYOFYR,y=mean_ERC,group=Climate,color=Climate)) +
#  geom_text(aes(0,Q80),label = "historical 80th percentile", hjust= 0,vjust = -1) +
#  theme_bw()+ xlab("Day of Year") + ylab("Mean ERC")

ggplot() +
  geom_hline(yintercept=Q80, linetype="dashed", color = "black") +
  geom_line(data=risQ_ERC_mean,aes(x=DATE,y=mean_ERC,group=Climate,color=Climate)) +
  scale_x_date(limits = as.Date(c("2020-01-01","2020-12-30")), date_breaks = "1 month",date_labels = "%b") +
  annotate(geom="text", x=as.Date("2020-01-01"),y=Q80,label = "historical 80th percentile", hjust= 0,vjust = -1) +
  theme_bw()+ xlab("Day of Year") + ylab("Mean ERC")
