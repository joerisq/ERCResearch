###############################################################################################################################
##                                           condition_pyrome.R
##
## last modified: 02-04-2020
## 
## Modified from Evan's pilot_input_modifictions.R
##
## Function to modify time series for a single weather station/grid point.
##
##  - currently only modifying temperature and RH datasets based on tas_mean and rhs_mean simulations   
##  - currently only option is to modify based on rcp85, simulations haven't been run yet for rcp45
##  - fixed units --> change factors for temperature now converted to def F
## 
##  Inputs to Function
##    plat: latitude do station or grid point
##    plon: longitude do station or grid point
##    Year: the year you want to condition for (at the moment 2030, 2040, 2045, 2050)
##    scenario: "rcp45" or "rcp85" (currenlty only works for rcp85)
##    scnenario_obj: "lwr", "med" , or "upr" for best case, medium case, or worst case scenario
##    pyrome_data: dataframe containing the data to be modified
##    change_factors: list read in from the climate conditioning filefile
##
###############################################################################################################################


condition_pyrome <- function(plat = 45, plon = -110, Year = 2040, scenario = "rcp85", scenario_obj = "med",
                             pyrome_data = data_tb, change_factors = change_factors_file){
  
  # Specify plat and plon if running a test case
  #plon <- -116.9086111
  #plat <- 46.92972222
  
  # FUNCTION TO FIND GRID POINT CLOSEST TO PYROME CENTRE
  find_gridpt <- function(plat = 45, plon = -115, risq_grid = grid_pts){
    
    dist <- (plat-risq_grid$lat)^2 + (plon-risq_grid$lon)^2
    output <- risq_grid$risq_2degree_id[which.min(dist)]
    return(output)
  }
  
  # get the risq grid ID number associated with the centroid of the pyrome
  grid_pts <- change_factors$risQ_grid
  grid_pt <- find_gridpt(plat, plon, risq_grid = grid_pts)
  
  
  # extract T and RH factos based on scenario
  if (scenario_obj=="upr"){
    RH_change_factors <- change_factors[[paste0("change_factors_",Year,"_Q_0.025")]][[grid_pt]]
    T_change_factors <- change_factors[[paste0("change_factors_",Year,"_Q_0.975")]][[grid_pt]]  
  } else if (scenario_obj=="lwr"){
    RH_change_factors <- change_factors[[paste0("change_factors_",Year,"_Q_0.975")]][[grid_pt]]
    T_change_factors <- change_factors[[paste0("change_factors_",Year,"_Q_0.025")]][[grid_pt]]
  } else if (scenario_obj=="med"){
    RH_change_factors <- change_factors[[paste0("change_factors_",Year,"_Q_0.5")]][[grid_pt]]
    T_change_factors <- change_factors[[paste0("change_factors_",Year,"_Q_0.5")]][[grid_pt]]
  }
  
  
  # convert the pyrome data to numeric
  setDF(pyrome_data)
  cols.num <- c("TEMP","MinT","MaxT","RH","MxRH","MnRH")
  pyrome_data[cols.num] <- sapply(pyrome_data[cols.num],as.numeric)
  
  # add column that contain the month as a numeric
  pyrome_data$MONTH_int <- as.numeric(pyrome_data$MONTH)
  
  ## create a dataframe with the station data, and the change factors for all the variables 
  pyrome_df <- pyrome_data %>% 
    dplyr::inner_join(T_change_factors[[paste0("delta_tas_mean_",scenario)]] %>%
                        dplyr::select(Season, tas_mean_delta = delta), 
                      by = c("MONTH_int"="Season")) %>%
    #dplyr::inner_join(grid_pt_scenario[[paste0("delta_pr_mean_",scenario)]] %>%
    #                    dplyr::select(Season, pr_mean_delta = delta), 
    #                  by = c("MONTH_int"="Season")) %>%
    dplyr::inner_join(RH_change_factors[[paste0("delta_rhs_mean_",scenario)]] %>%
                        dplyr::select(Season, rhs_mean_delta = delta), 
                      by = c("MONTH_int"="Season")) %>%
    dplyr::mutate(index = 1:nrow(.)) # create variable "index" which numbers each row --> order of time series
  
  ## add average temperature and RH change to the original 1pm temperature and RH field (there are no avg values)
  pyrome_df$new_Temp <- pyrome_df$TEMP + (9/5)*pyrome_df$tas_mean_delta # convert change to deg F
  pyrome_df$new_RH <- pyrome_df$RH + pyrome_df$rhs_mean_delta 
  
  ## modify 1pm, min, and max temperautre and RH based on originL difference from daily average
  pyrome_df$new_MinT <- pyrome_df$MinT + (9/5)*pyrome_df$tas_mean_delta # convert change to deg F
  pyrome_df$new_MaxT <- pyrome_df$MaxT + (9/5)*pyrome_df$tas_mean_delta # convert change to deg F
  
  pyrome_df$new_MnRH <- pyrome_df$MnRH + pyrome_df$rhs_mean_delta
  pyrome_df$new_MxRH <- pyrome_df$MxRH + pyrome_df$rhs_mean_delta
  
  ## plot the change factors to compare with Riley and Loehman (2016) in deg C and % RH
  # delta_tas <- T_change_factors$delta_tas_mean_rcp85 
  # elta_rhs <- RH_change_factos$delta_rhs_mean_rcp85
  # 
  # ggplot(data = delta_tas, aes(x = Season, y = delta)) + 
  #   geom_point(color = "black", size=0.65) +
  #   geom_line() +
  #   scale_x_continuous(breaks=1:12, labels=c("1" = "Jan", "2" = "Feb", "3" = "March", "4" = "April", "5" = "May", "6" = "June", "7" = "July",
  #                             "8" = "Aug", "9" = "Sept", "10" = "Oct", "11" = "Nov", "12" = "Dec")) +
  #   labs(x="Month",y=expression('Temperature Change (deg C)')) +
  #   theme_bw()  +
  #   theme(text = element_text(size=20))
  # ggsave("/Users/carlinghay/Documents/risQ/ERCTimeSeriesAnalysis/Tas_Change_IdahoStation.png")
  # 
  # ggplot(data = delta_rhs, aes(x = Season, y = delta)) + 
  #   geom_point(color = "black", size=0.65) +
  #   geom_line() +
  #   scale_x_continuous(breaks=1:12, labels=c("1" = "Jan", "2" = "Feb", "3" = "March", "4" = "April", "5" = "May", "6" = "June", "7" = "July",
  #                                            "8" = "Aug", "9" = "Sept", "10" = "Oct", "11" = "Nov", "12" = "Dec")) +
  #   labs(x="Month",y=expression('RH Change (%)')) +
  #   theme_bw()  +
  #   theme(text = element_text(size=20))
  # ggsave("/Users/carlinghay/Documents/risQ/ERCTimeSeriesAnalysis/RH_Change_IdahoStation.png")
  
  ## plot the original and modified time series
  # library(lubridate)
  # RH_data <- pyrome_df %>% dplyr::select(YEAR,MONTH,DAY,RH,new_RH,rhs_mean_delta) %>%
  #   mutate(DATE = ymd(paste(YEAR, MONTH, DAY, sep = ' ')),
  #          delta = new_RH - RH)
  # 
  # tas_data <- pyrome_df %>% dplyr::select(YEAR,MONTH,DAY,TEMP,new_Temp,tas_mean_delta) %>%
  #   mutate(DATE = ymd(paste(YEAR, MONTH, DAY, sep = ' ')),
  #          delta = new_Temp - TEMP)
  # 
  # ggplot() +
  #   # plot original value
  #   #geom_line(data = RH_data[1:365,], aes(x = DATE, y = RH), color = "gray45", size = 0.65) +
  # 
  #   # plot the change in temperature in deg F
  #   geom_line(data = RH_data[1:365,], aes(x = DATE, y = delta), color = "gray45", size = 0.65) +
  # 
  #   # plot modified value
  #  #geom_line(data = RH_data[1:365,], aes(x = DATE, y = new_RH), color = "red", size=0.65) +
  # 
  #   #labs(x="Date",y=expression('Temperature Change (F)')) +
  #   labs(x="Date",y=expression('RH Change')) +
  #   #labs(x="Date",y=expression('RH (%)')) +
  #   theme_bw()  +
  #   theme(text = element_text(size=14))
  
  ## convert the variables to integers 
  cols.num <- c("new_Temp","new_MinT","new_MaxT","new_RH","new_MxRH","new_MnRH")
  pyrome_df[cols.num] <- sapply(pyrome_df[cols.num],as.integer)
  
  ## convert the columns to characters
  pyrome_df[cols.num] <- sapply(pyrome_df[cols.num],as.character)
  
  ## select only the modified columns to output
  out <- pyrome_df %>%
    dplyr::select(-tas_mean_delta,-rhs_mean_delta, -TEMP, -RH, -MinT, -MaxT, -MnRH,-MxRH, -MONTH_int, -index) %>%
    dplyr::rename(TEMP = new_Temp,
                  RH = new_RH,
                  MinT = new_MinT,
                  MaxT = new_MaxT,
                  MnRH = new_MnRH,
                  MxRH = new_MxRH)
  
  # drop the column names we don't want
  pyrome_data$index <- NULL
  pyrome_data$MONTH_int <- NULL
  
  require(data.table) 
  out <- setDT(out[colnames(pyrome_data)])
  
  return(out)
  
  # ## 
  # pyrome_tbl_mod_max <- zoo_tbl_mod %>%
  #   dplyr::group_by(Year, Season) %>%
  #   # get the index (to know the timing) of the maximum precip over 24, 48 72, and 96 hr in historical data
  #   dplyr::summarise(max_24_index = index[which.max(RcppRoll::roll_sum(P, n = 1))], # compute rolling sum of time series
  #                    max_48_index = index[which.max(RcppRoll::roll_sum(P, n = 2))],
  #                    max_72_index = index[which.max(RcppRoll::roll_sum(P, n = 3))],
  #                    max_96_index = index[which.max(RcppRoll::roll_sum(P, n = 4))], 
  #                    
  #                    # get the maximum of the 24, 48, 72, 96 hr precip values in historical data
  #                    max_24 = max(RcppRoll::roll_sum(P, n = 1)),
  #                    max_48 = max(RcppRoll::roll_sum(P, n = 2)),
  #                    max_72 = max(RcppRoll::roll_sum(P, n = 3)),
  #                    max_96 = max(RcppRoll::roll_sum(P, n = 4)), 
  #                    
  #                    # compute the quantiles of the daily precipitation for each season/year historically
  #                    q_75 = quantile(P, probs=0.75), 
  #                    q_90 = quantile(P, probs=0.90),
  #                    q_95 = quantile(P, probs=0.95),
  #                    q_99 = quantile(P, probs=0.99)) %>%
  #   dplyr::ungroup()
  
  # seasonal_processing <- function(tmp, modif){
  # 
  # 
  # 
  #   ## determine which statistic the historical precipation falls into, then modify it based on that statistics projected change
  #   ### THIS HAS BEEN MODIFIED FROM ORIGINAL CODE, NEED MOST RESTRICTIVE CONDItion FIRST
  #   tmpnew <- tmp %>% 
  #     mutate(P = case_when(P >= modif$q_99 ~ P*tmp$pr_q_99_delta[1], # new_P = hist_P * delta_raio_of_change
  #                          P >= modif$q_95 ~ P*tmp$pr_q_95_delta[1], 
  #                          P >= modif$q_90 ~ P*tmp$pr_q_90_delta[1],
  #                          P >= modif$q_75 ~ P*tmp$pr_q_75_delta[1],
  #                          P == modif$max_24 ~ P*tmp$pr_max_24_delta,   # is ths more important than the abvove conditions?
  #                          P < modif$q_75 ~ P) # if less than 75th percentile, keep the same
  #            )
  #   # now figure out how to adjust mean precip
  #   target_total <- tmp$pr_mean_delta[1]*sum(tmp$P) # seasonal mean change ratio x total historical precip in the season
  #   current_total <- sum(tmpnew$P) # projected total by adjusting extremes only
  #   
  #   # diff & ratio b/w what the total should be based on projecting using mean, and total by projecting extremes
  #   margin <- target_total - current_total 
  #   adj_ratio <- target_total/current_total
  #   
  # 
  #   adj_ind <- which(tmpnew$P < modif$q_75 & tmpnew$P > 0) # find the location where the precip is less that the 75th percentil
  #   adj <- margin/length(adj_ind) # find the precip amount per day (i.e. how much more precip is needed per day)
  #   tmpnew$P[adj_ind] <- tmpnew$P[adj_ind] + adj # adjust the non-extreme days to bring the sum of  projected precip values into agreement
  # 
  #   
  #   tmpnew$P[tmpnew$P < 0] <- 0
  #   
  #   # add the temperature change to the original temperature field
  #   tmpnew$E <- tmp$E + tmp$tas_mean_delta[1] # in degrees C
  # 
  #   #cat(tmp$Year[1], tmp$Season[1],  min(tmpnew$P-tmp$P), "\n")
  #   
  #   return(tmpnew %>% 
  #            dplyr::select(Date, P, E))
  # }
  # 
  # ## split the data so that each season-year is its own list
  # zoo_tbl_mod_max_split <- zoo_tbl_mod_max %>%
  #   mutate(grp = paste(Year, Season, sep="_")) %>% 
  #   split(f = .$grp)
  #   
  # zoo_tbl_mod_split <- zoo_tbl_mod %>%
  #   mutate(grp = paste(Year, Season, sep="_")) %>% 
  #   split(f = .$grp)
  # 
  # ## perform the daily seasonal adjustment for each season-year list
  # zoo_out <- purrr::map2(.x = zoo_tbl_mod_split, 
  #             .y = zoo_tbl_mod_max_split, 
  #             .f = ~seasonal_processing(.x, .y)) %>%
  #   data.table::rbindlist() %>%
  #   inner_join(zoo_tbl %>% dplyr::rename(P_old = P, E_old = E), by = "Date")
  # 
  # #print(max(zoo_out$P - zoo_out$P_old))
}
