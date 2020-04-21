###############################################################################################################################
##                                           condition_pyrome_3.1_V01.R
##
## last modified: 03-31-2020
## 
## Modified from condition_pyrome_v3.R
##
## Function to modify time series for a single weather station/grid point.
##
##  - modifies temperature, RH, Rain, and RnDr variables
##  - can be run for rcp 4.5, 8.5 and any decade from 2020 through 2100, including 2045
##  - all scenario_obj using mediam precip projecyion
##
##  - computes new_TEMP, MaxT and MinT by shifting all by the mean temperature change
##  - computes new_RH, MxRH and MnRH by shifting all by the mean RH change
##
##  - computes new_Rain two ways: 
##        new_Rain - uses absolute changes, 
##        new_Rain2 - uses percent changes
##  - computes new_RnDr 3 ways: 
##        new_RnDr - keeps intensity the same, new_RnDr = new_Rain2/intensity
##        new_RnDr2 - uses relationship b/w historical Rain and RnDr, then esimates RnDr using new_Rain2
##        new_RnDr3 - uses relationship b/w historical Rain and RnDr, then esimates RnDr using new_Rain
##
##  Inputs to Function
##    plat: latitude do station or grid point
##    plon: longitude do station or grid point
##    Year: the year you want to condition for (2020, 2030, 2040, 2045, 2050, ... 2100)
##    scenario: "rcp45" or "rcp85" 
##    scnenario_obj: "lwr", "med" , or "upr" for best case, medium case, or worst case scenario,
##    pyrome_data: dataframe containing the data to be modified
##    change_factors_tas: list of temperature climate conditioning change factors
##    change_factors_rh: list of RH climate conditioning change factors
##    change_factors_pr: list of precipitation climate conditioning change factors, must have percent changes
##
##  Outputs:
##    A table with column names matching the input, but values modified under climate conditioning. Output Rain is 
##    new_Rain2 and outut RnDr is new_RnDr2
##
##
###############################################################################################################################


condition_pyrome_3.1_V01 <- function(plat = 45, plon = -110, Year = 2040, scenario = "rcp85", scenario_obj = "med",
                             pyrome_data = data_tb, change_factors_tas = tas_change_factors_file, 
                             change_factors_rh = rh_change_factors_file, change_factors_pr = pr_change_factors_file){
  
  # FUNCTION TO FIND GRID POINT CLOSEST TO PYROME CENTRE
  find_gridpt <- function(plat = 45, plon = -115, risq_grid = grid_pts){
    
    dist <- (plat-risq_grid$lat)^2 + (plon-risq_grid$lon)^2
    output <- risq_grid$risq_2degree_id[which.min(dist)]
    return(output)
  }
  
  # get the risq grid ID number associated with the centroid of the pyrome
  grid_pts <- change_factors_tas$risQ_grid
  grid_pt <- find_gridpt(plat, plon, risq_grid = grid_pts)
  
  
  # extract T, Precip, and RH factors based on scenario -->  all scenarios use med proj changes
  if (scenario_obj=="upr"){
    RH_change_factors <- change_factors_rh[[paste0("change_factors_",Year,"_Q_0.025")]][[grid_pt]]
    T_change_factors <- change_factors_tas[[paste0("change_factors_",Year,"_Q_0.975")]][[grid_pt]]  
    P_change_factors <- change_factors_pr[[paste0("change_factors_",Year,"_Q_0.5")]][[grid_pt]]
  } else if (scenario_obj=="lwr"){
    RH_change_factors <- change_factors_rh[[paste0("change_factors_",Year,"_Q_0.975")]][[grid_pt]]
    T_change_factors <- change_factors_tas[[paste0("change_factors_",Year,"_Q_0.025")]][[grid_pt]]
    P_change_factors <- change_factors_pr[[paste0("change_factors_",Year,"_Q_0.5")]][[grid_pt]]
  } else if (scenario_obj=="med"){
    RH_change_factors <- change_factors_rh[[paste0("change_factors_",Year,"_Q_0.5")]][[grid_pt]]
    T_change_factors <- change_factors_tas[[paste0("change_factors_",Year,"_Q_0.5")]][[grid_pt]]
    P_change_factors <- change_factors_pr[[paste0("change_factors_",Year,"_Q_0.5")]][[grid_pt]]
  }
  
  
  # convert the pyrome data to numeric
  setDF(pyrome_data)
  cols.num <- c("TEMP","MinT","MaxT","RH","MxRH","MnRH","Rain","RnDr")
  pyrome_data[cols.num] <- sapply(pyrome_data[cols.num],as.numeric)
  
  # add column that contain the month as a numeric 
  pyrome_data <- pyrome_data %>% 
    dplyr::mutate(MONTH_int = as.numeric(MONTH))
  
  ## create a dataframe with the station data, the precip. intensity, and the change factors for all the variables 
  pyrome_df <- pyrome_data %>% 
    dplyr::mutate(RnIntensity = Rain/RnDr) %>%
    dplyr::inner_join(T_change_factors[[paste0("delta_tas_mean_",scenario)]] %>%
                        dplyr::select(Season, tas_mean_delta = delta), 
                      by = c("MONTH_int"="Season")) %>%
    dplyr::inner_join(P_change_factors[[paste0("delta_pr_q_75_",scenario)]] %>%
                        dplyr::select(Season, pr_q_75_delta = delta, pr_q_75_percent_change = percent_change),
                      by = c("MONTH_int"="Season")) %>%
    dplyr::inner_join(P_change_factors[[paste0("delta_pr_q_80_",scenario)]] %>%
                        dplyr::select(Season, pr_q_80_delta = delta, pr_q_80_percent_change = percent_change),
                      by = c("MONTH_int"="Season")) %>%
    dplyr::inner_join(P_change_factors[[paste0("delta_pr_q_90_",scenario)]] %>%
                        dplyr::select(Season, pr_q_90_delta = delta, pr_q_90_percent_change = percent_change),
                      by = c("MONTH_int"="Season")) %>%
    dplyr::inner_join(P_change_factors[[paste0("delta_pr_q_95_",scenario)]] %>%
                        dplyr::select(Season, pr_q_95_delta = delta, pr_q_95_percent_change = percent_change),
                      by = c("MONTH_int"="Season")) %>%
    dplyr::inner_join(P_change_factors[[paste0("delta_pr_q_99_",scenario)]] %>%
                        dplyr::select(Season, pr_q_99_delta = delta, pr_q_99_percent_change = percent_change),
                      by = c("MONTH_int"="Season")) %>%
    dplyr::inner_join(RH_change_factors[[paste0("delta_rhs_mean_",scenario)]] %>%
                        dplyr::select(Season, rhs_mean_delta = delta), 
                      by = c("MONTH_int"="Season")) %>%
    dplyr::mutate(index = 1:nrow(.)) # create variable "index" which numbers each row --> order of time series
  
  ## insure that the rain intensity makes sense, if RnDr = 0, then set intensity to 0
  pyrome_df$RnIntensity[pyrome_df$RnDr==0] <- 0
  
  ##################################################################################################################
  ##           modify temperature and RH using mean projected change in the variable
  ##################################################################################################################
  
  ## add average temperature and RH change to the original 1pm temperature and RH field (there are no avg values)
  pyrome_df$new_Temp <- pyrome_df$TEMP + (9/5)*pyrome_df$tas_mean_delta # convert change to deg F
  pyrome_df$new_RH <- pyrome_df$RH + pyrome_df$rhs_mean_delta 
  
  ## modify 1pm, min, and max temperautre and RH based on originL difference from daily average
  pyrome_df$new_MinT <- pyrome_df$MinT + (9/5)*pyrome_df$tas_mean_delta # convert change to deg F
  pyrome_df$new_MaxT <- pyrome_df$MaxT + (9/5)*pyrome_df$tas_mean_delta # convert change to deg F
  
  pyrome_df$new_MnRH <- pyrome_df$MnRH + pyrome_df$rhs_mean_delta
  pyrome_df$new_MxRH <- pyrome_df$MxRH + pyrome_df$rhs_mean_delta
  
  ##################################################################################################################
  ##       modify Rain 
  ##        method 1 (new_Rain): use absolute change of precip quantiles
  ##        method 2 (new_Rain2): use percent change of precip quantiles
  ##################################################################################################################
  
  ## modify precipitation using the quantiles
  pyrome_df_quantile <- pyrome_df %>%
    dplyr::group_by(YEAR, MONTH_int) %>%
    # get the index (to know the timing) of the maximum precip over 24, 48 72, and 96 hr in historical data
    dplyr::summarise(q_75 = quantile(Rain, probs=0.75),
                     q_80 = quantile(Rain, probs=0.80),
                     q_90 = quantile(Rain, probs=0.90),
                     q_95 = quantile(Rain, probs=0.95),
                     q_99 = quantile(Rain, probs=0.99)) %>%
    dplyr::ungroup()
  
  # function to modify precipitation
  monthly_processing <- function(tmp, modif){
    
    ## add the new rain amount for each quantile, but first convert from mm to inches
    #  new_Rain: uses absolutation change in precip quantile
    #  new_Rain2:  uses the percent change in the precip quantile
    tmpnew <- tmp %>% 
      mutate(new_Rain = case_when(Rain >= modif$q_99 ~ Rain + (0.0393701)*tmp$pr_q_99_delta, 
                                  Rain >= modif$q_95 ~ Rain + (0.0393701)*tmp$pr_q_95_delta, 
                                  Rain >= modif$q_90 ~ Rain + (0.0393701)*tmp$pr_q_90_delta,
                                  Rain >= modif$q_80 ~ Rain + (0.0393701)*tmp$pr_q_80_delta,
                                  Rain >= modif$q_75 ~ Rain + (0.0393701)*tmp$pr_q_75_delta,
                                  Rain < modif$q_75 ~ Rain), # if less than 75th percentile, keep the same
             new_Rain2 = case_when(Rain >= modif$q_99 ~ Rain*tmp$pr_q_99_percent_change, 
                                  Rain >= modif$q_95 ~ Rain*tmp$pr_q_95_percent_change, 
                                  Rain >= modif$q_90 ~ Rain*tmp$pr_q_90_percent_change,
                                  Rain >= modif$q_80 ~ Rain*tmp$pr_q_80_percent_change,
                                  Rain >= modif$q_75 ~ Rain*tmp$pr_q_75_percent_change,
                                  Rain < modif$q_75 ~ Rain)
      )
    
    # ensure no negative precip
    tmpnew$new_Rain[tmpnew$new_Rain < 0] <- 0
    tmpnew$new_Rain2[tmpnew$new_Rain2 < 0] <- 0
    
    # if there was no precip oroiginally, set to zero after climate conditioning
    #  - a problem b/c some quantiles may be zero is there were few rainy days historically
    tmpnew$new_Rain[tmpnew$Rain == 0] <- 0
    tmpnew$new_Rain2[tmpnew$Rain2 == 0] <- 0
    
    return(tmpnew %>% 
             dplyr::select(YEAR,MONTH,DAY,new_Rain, new_Rain2))
  }
  
  ## split the data so that each season-year is its own list
  pyrome_df_quantile_split <- pyrome_df_quantile %>%
    mutate(grp = paste(YEAR, MONTH_int, sep="_")) %>% 
    split(f = .$grp)
  
  pyrome_df_split <- pyrome_df %>%
    mutate(grp = paste(YEAR, MONTH_int, sep="_")) %>% 
    split(f = .$grp)
  
  ## perform the daily seasonal adjustment for each season-year list
  pyrome_df <- pyrome_df %>%
    inner_join(purrr::map2(.x = pyrome_df_split, 
                           .y = pyrome_df_quantile_split, 
                           .f = ~monthly_processing(.x, .y)) %>%
                 data.table::rbindlist(), by = c("YEAR","MONTH","DAY"))
  
  
  ##################################################################################################################
  ##                  modify RnDr 
  ## new_RnDr: keeps the intensity the same as before, new_RnDr = new_Rain2/RnIntensity
  ## new_RnDr2: uses the linar relationship b/w historical Rain and RnDr to adjust new_RnDr based on the new_Rain2 (percent change)
  ## new_RnDr3: uses the linar relationship b/w historical Rain and RnDr to adjust new_RnDr based on the new_Rain (absolute change)
  ##
  ##################################################################################################################
  
  # keep the original intensity, compute new duration as amount of precip/itensity
  pyrome_df <- pyrome_df %>%
    mutate(new_RnDr = new_Rain2/(RnIntensity))
  
  # if the climate conditioned RnDr is less than 30min, set to zero
  pyrome_df$new_RnDr[pyrome_df$new_RnDr<=0.5] <- 0
  
  # keep all days with original RnDr = 0 as zero
  pyrome_df$new_RnDr[pyrome_df$RnDr==0] <- 0
  
  # the RnDr can't be greater than 24 hr
  pyrome_df$new_RnDr[pyrome_df$new_RnDr>24] <- 24
  
  # check for any NaN's, set to zero
  pyrome_df$new_RnDr[is.na(pyrome_df$new_RnDr)] <- 0
  
  # round all RnDr values up to the nearest integer
  pyrome_df$new_RnDr <- ceiling(pyrome_df$new_RnDr)
  
  ## function to determine RnDr vs Rain trends for the second method, force it to have no intercept (i.. RnDr=0 when Rain=0)
  get_coef <- function(x,y){
    MASS::rlm(y ~ x +0, maxit=200, data = data.frame(x,y)) %>%  # increase maxint from default of 20 to ensure convergence
      summary(.) %>%
      coef(.) %>% 
      .[1,1]
  }
  
  # function to modify new_RnDr2 and new_RnDr3
  monthly_RnDr_processing <- function(tmp, modif){
    
    # add the change in Rain Duration based on monthly relationship b/w RnDr and Rain
    tmpnew <- tmp %>% 
      mutate(new_RnDr2 = RnDr + (new_Rain2 - Rain)*modif$trend,
             new_RnDr3 = RnDr + (new_Rain - Rain)*modif$trend)
    
    # keep all days with original RnDr = 0 as zero
    tmpnew$new_RnDr2[tmpnew$RnDr2==0] <- 0
    tmpnew$new_RnDr3[tmpnew$RnDr3==0] <- 0
    
    # the RnDr can't be greater than 24 hr
    tmpnew$new_RnDr2[tmpnew$new_RnDr2>24] <- 24
    tmpnew$new_RnDr3[tmpnew$new_RnDr3>24] <- 24
    
    # set any RnDr that is 30min or less to zero
    tmpnew$new_RnDr2[tmpnew$new_RnDr2<=0.5] <- 0
    tmpnew$new_RnDr3[tmpnew$new_RnDr3<=0.5] <- 0
    
    # round all RnDr to the closest hour
    tmpnew$new_RnDr2 <- ceiling(tmpnew$new_RnDr2)
    tmpnew$new_RnDr3 <- ceiling(tmpnew$new_RnDr3)
    
    return(tmpnew %>% 
             dplyr::select(YEAR,MONTH,DAY,new_RnDr2,new_RnDr3))
  }
 
   ## compute the trends RnDr vs Rain for every month
  pyrome_df_trends <- pyrome_df %>%
    dplyr::group_by(MONTH_int) %>%
    dplyr::summarise(trend = get_coef(Rain,RnDr)) %>%
    dplyr::ungroup()
  
  ## split the data so that each Month is its own list
  pyrome_df_trends_split <- pyrome_df_trends %>%
    mutate(grp = MONTH_int) %>% 
    split(f = .$grp)
  
  pyrome_df_split <- pyrome_df %>%
    mutate(grp = MONTH_int) %>% 
    split(f = .$grp)
  
  ## perform the daily seasonal adjustment for each season-year list
  pyrome_df <- pyrome_df %>%
    inner_join(purrr::map2(.x = pyrome_df_split, 
                           .y = pyrome_df_trends_split, 
                           .f = ~monthly_RnDr_processing(.x, .y)) %>%
                 data.table::rbindlist(), by = c("YEAR","MONTH","DAY"))
  
  
  ##################################################################################################################
  ##                  make some plots if desired
  ##################################################################################################################
  
  if (FALSE){
    ## plots of Rain vs RnDr or Intensity
    p1 = ggplot(data=pyrome_df %>% filter(MONTH_int==1),aes(x=Rain,y=RnDr)) + geom_point() + geom_smooth(method = "lm") + 
      labs(title="January")
    p2 = ggplot(data=pyrome_df %>% filter(MONTH_int==2),aes(x=Rain,y=RnDr)) + geom_point() + geom_smooth(method = "lm") + 
      labs(title="February")
    p3 = ggplot(data=pyrome_df %>% filter(MONTH_int==3),aes(x=Rain,y=RnDr)) + geom_point() + geom_smooth(method = "lm") + 
      labs(title="March")
    p4 = ggplot(data=pyrome_df %>% filter(MONTH_int==4),aes(x=Rain,y=RnDr)) + geom_point() + geom_smooth(method = "lm") + 
      labs(title="April")
    p5 = ggplot(data=pyrome_df %>% filter(MONTH_int==5),aes(x=Rain,y=RnDr)) + geom_point() + geom_smooth(method = "lm") + 
      labs(title="May")
    p6 = ggplot(data=pyrome_df %>% filter(MONTH_int==6),aes(x=Rain,y=RnDr)) + geom_point() + geom_smooth(method = "lm") + 
      labs(title="June")
    
    
    g <- arrangeGrob(p1, p2, p3, p4, p5,p6, nrow = 3)
    ggsave(paste0(ROOT_DATA_DIR,'/',OUTPUT_DIR,'/',"JanuaryToJune_RainRnIntensity_pyrome26.png"),g)
    
     ## plot the change factors to compare with Riley and Loehman (2016) in deg C and % RH
    delta_tas <- T_change_factors$delta_tas_mean_rcp85
    delta_rhs <- RH_change_factors$delta_rhs_mean_rcp85
    delta_pr_mean <- P_change_factors$delta_pr_mean_rcp85
    delta_pr_q_99 <- P_change_factors$delta_pr_q_99_rcp85
    delta_pr_q_95 <- P_change_factors$delta_pr_q_95_rcp85
    delta_pr_q_90 <- P_change_factors$delta_pr_q_90_rcp85
    delta_pr_q_80 <- P_change_factors$delta_pr_q_80_rcp85
    delta_pr_q_75 <- P_change_factors$delta_pr_q_75_rcp85
    
    ## compute the monthly precipitation totals and the difference
    pyrome_rain <- pyrome_df %>%
      dplyr::group_by(YEAR, MONTH_int) %>%
      #dplyr::group_by(YEAR) %>%
      # get the index (to know the timing) of the maximum precip over 24, 48 72, and 96 hr in historical data
      dplyr::summarise(total_rain = sum(Rain),
                       new_total_rain = sum(new_Rain),
                       new_total_rain2 = sum(new_Rain2),
                       rain_diff = sum(new_Rain) - sum(Rain),
                       rain_diff2 = sum(new_Rain2) - sum(Rain)) %>%
      dplyr::ungroup() %>%
      filter(MONTH_int==12)
      #filter(YEAR=="90")
    
    
    # plot difference in momntly precip totals for a given year
    ggplot(data = pyrome_rain %>% filter(YEAR=="90"), aes(x = MONTH_int, y = rain_diff)) + 
      geom_point(color = "black", size=0.65) +
      geom_line() +
      scale_x_continuous(breaks=1:12, labels=c("1" = "Jan", "2" = "Feb", "3" = "March", "4" = "April", "5" = "May", "6" = "June", "7" = "July",
                                               "8" = "Aug", "9" = "Sept", "10" = "Oct", "11" = "Nov", "12" = "Dec")) +
      labs(x="Month",y=expression('Total Precip. Change (inches)')) +
      theme_bw()  +
      theme(text = element_text(size=20))
    
    ## plots!
    ggplot(data = delta_pr_q_90 %>% mutate(percent_change = 100*(percent_change-1)), aes(x = Season, y = percent_change)) +
      geom_point(color = "black", size=0.65) +
      geom_line() +
      scale_x_continuous(breaks=1:12, labels=c("1" = "Jan", "2" = "Feb", "3" = "March", "4" = "April", "5" = "May", "6" = "June", "7" = "July",
                                               "8" = "Aug", "9" = "Sept", "10" = "Oct", "11" = "Nov", "12" = "Dec")) +
      labs(x="Month",y=expression('Change in 0.90 Precip. Quantile')) +
      theme_bw()  +
      theme(text = element_text(size=20))
    #ggsave("/Users/carlinghay/Documents/risQ/ERCTimeSeriesAnalysis/Figures/pyrome8_pr_q_95_percent_change.png")
    
    ggplot(data = delta_rhs, aes(x = Season, y = delta)) +
      geom_point(color = "black", size=0.65) +
      geom_line() +
      scale_x_continuous(breaks=1:12, labels=c("1" = "Jan", "2" = "Feb", "3" = "March", "4" = "April", "5" = "May", "6" = "June", "7" = "July",
                                               "8" = "Aug", "9" = "Sept", "10" = "Oct", "11" = "Nov", "12" = "Dec")) +
      labs(x="Month",y=expression('RH Change (%)')) +
      theme_bw()  +
      theme(text = element_text(size=20))
    # ggsave("/Users/carlinghay/Documents/risQ/ERCTimeSeriesAnalysis/RH_Change_IdahoStation.png")
    
    ggplot(data = delta_pr_q_75, aes(x = Season, y = delta)) +
      geom_point(color = "black", size=0.65) +
      geom_line() +
      scale_x_continuous(breaks=1:12, labels=c("1" = "Jan", "2" = "Feb", "3" = "March", "4" = "April", "5" = "May", "6" = "June", "7" = "July",
                                               "8" = "Aug", "9" = "Sept", "10" = "Oct", "11" = "Nov", "12" = "Dec")) +
      labs(x="Month",y=expression('Precip. Q 0.95 Change')) +
      theme_bw()  +
      theme(text = element_text(size=20))
    
    
    ## plot the original and modified time series
    
    RH_data <- pyrome_df %>% dplyr::select(YEAR,MONTH,DAY,RH,new_RH,rhs_mean_delta) %>%
      mutate(DATE = ymd(paste(YEAR, MONTH, DAY, sep = ' ')),
             delta = new_RH - RH) %>%
      filter(MONTH=='07' | MONTH=='08' | MONTH=='09', YEAR=='14')
    
    tas_data <- pyrome_df %>% dplyr::select(YEAR,MONTH,DAY,TEMP,new_Temp,tas_mean_delta) %>%
      mutate(DATE = ymd(paste(YEAR, MONTH, DAY, sep = ' ')),
             delta = new_Temp - TEMP)  %>%
      filter(MONTH=='07' | MONTH=='08' | MONTH=='09', YEAR=='14')
    
    pr_data <- pyrome_df %>% dplyr::select(YEAR,MONTH,DAY,Rain,new_Rain, new_Rain2, RnDr, new_RnDr, new_RnDr2,new_RnDr3) %>%
      mutate(DATE = ymd(paste(YEAR, MONTH, DAY, sep = ' '))) %>%
      #filter(YEAR=='12' | YEAR=='13' | YEAR=='14')
      filter(MONTH=='07' | MONTH=='08' | MONTH=='09', YEAR=='14')
    
    # Temperatue
    ggplot() +
      # plot original value
      geom_line(data = tas_data, aes(x = DATE, y = TEMP), color = "gray45", size = 0.65) +
      
      # plot modified value
      #geom_line(data = pr_data, aes(x = DATE, y = new_Rain), color = "blue", size=0.65) +
      geom_line(data = tas_data, aes(x = DATE, y = new_Temp), color = "red", size=0.65) +
      
      labs(x="Date",y=expression('Temperature (F)')) +
      theme_bw()  +
      theme(text = element_text(size=14))
    
    ggsave("/Users/carlinghay/Documents/risQ/ERCTimeSeriesAnalysis/Figures/pyrome8_tas_JulytoSept.png")
    
    # RH
    ggplot() +
      # plot original value
      geom_line(data = RH_data, aes(x = DATE, y = RH), color = "gray45", size = 0.65) +
      
      # plot modified value
      #geom_line(data = pr_data, aes(x = DATE, y = new_Rain), color = "blue", size=0.65) +
      geom_line(data = RH_data, aes(x = DATE, y = new_RH), color = "red", size=0.65) +
      
      labs(x="Date",y=expression('RH (%)')) +
      theme_bw()  +
      theme(text = element_text(size=14))
    
    ggsave("/Users/carlinghay/Documents/risQ/ERCTimeSeriesAnalysis/Figures/pyrome8_RH_JulytoSept2014.png")
    
    # Rain
    ggplot() +
      # plot original value
      geom_line(data = pr_data, aes(x = DATE, y = Rain), color = "gray45", size = 0.65) +
      
      # plot modified value
      #geom_line(data = pr_data, aes(x = DATE, y = new_Rain), color = "blue", size=0.65) +
      geom_line(data = pr_data, aes(x = DATE, y = new_Rain2), color = "red", size=0.65) +
      
      labs(x="Date",y=expression('Precipitation (inches)')) +
      theme_bw()  +
      theme(text = element_text(size=14))
   
    ggsave("/Users/carlinghay/Documents/risQ/ERCTimeSeriesAnalysis/Figures/pyrome8_Rain_JulytoSept2014.png") 
    
    # RnDr 
    ggplot() +
      # plot original value
      geom_line(data = pr_data, aes(x = DATE, y = RnDr), color = "gray45", size = 0.65) +
      
      # plot modified value
      #geom_line(data = pr_data, aes(x = DATE, y = new_RnDr), color = "darkgreen", size=0.65) +
      geom_line(data = pr_data, aes(x = DATE, y = new_RnDr2), color = "red", size=0.65) +
      #geom_line(data = pr_data, aes(x = DATE, y = new_RnDr3), color = "darkgreen", size=0.65) +
      
      labs(x="Date",y=expression('Precipitation  Duration (hr)')) +
      theme_bw()  +
      theme(text = element_text(size=14))
    
    ggsave("/Users/carlinghay/Documents/risQ/ERCTimeSeriesAnalysis/Figures/pyrome8_RnDr_JulytoSept2014.png") 
  }
  
  
  ##################################################################################################################
  ##                  put the conditioned variables into a data frame in the correct form for outputting
  ##################################################################################################################
  
  ## convert the variables to integers 
  cols.num <- c("new_Temp","new_MinT","new_MaxT","new_RH","new_MxRH","new_MnRH","new_Rain","new_RnDr")
  pyrome_df[cols.num] <- sapply(pyrome_df[cols.num],as.integer)
  
  ## convert the columns to characters
  pyrome_df[cols.num] <- sapply(pyrome_df[cols.num],as.character)
  
  ## select only the modified columns to output
  out <- pyrome_df %>%
    dplyr::select(-tas_mean_delta,-rhs_mean_delta, -pr_q_75_delta, -pr_q_80_delta, -pr_q_90_delta, pr_q_95_delta, -pr_q_99_delta,
                  -TEMP, -RH, -MinT, -MaxT, -MnRH,-MxRH, -Rain,-new_Rain, -RnDr,-new_RnDr3, -RnIntensity, -MONTH_int, -index) %>%
    dplyr::rename(TEMP = new_Temp,
                  RH = new_RH,
                  MinT = new_MinT,
                  MaxT = new_MaxT,
                  MnRH = new_MnRH,
                  MxRH = new_MxRH,
                  Rain = new_Rain2,  # use new_Rain2 (modified using % change)
                  RnDr = new_RnDr2)  # use new_RnDr2 (modified using % change new_Rain and Rain vs RnDr relationship)
  
  # drop the column names we don't want
  pyrome_data$index <- NULL
  pyrome_data$MONTH_int <- NULL
  
  require(data.table) 
  out <- setDT(out[colnames(pyrome_data)])
  
  return(out)
}