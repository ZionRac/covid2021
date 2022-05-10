# This R code aims to create a function for Alberta in covid-19 death prediction
# Feb 20 2022
## tvlm led_death ~ cumulative cases + cvaccine
## data source: github UofT data, Summary data from www.Alberta.ca
## vaccination data, cases and death
# Created by Juxin liu, Dawo Zhang
#--------------------------------------
install.packages("devtools")
library(tidyverse)
library(devtools)
#library(covid19nytimes)
library(timetk)
library(lubridate)
library(broom)
library(knitr)
library(gt)
library(tvReg)


knitr::opts_chunk$set(echo = TRUE, eval = TRUE)

table_style <- list(cell_text(font = "Arial", size = "small"))
#open file

Adata <- "/Users/dawo/data/covid-19-alberta-statistics.csv"
ABdata = read.csv(Adata)
ABdata$Date = as.Date(ABdata$Date.reported.to.Alberta.Health, format="%Y-%m-%d")

urlvac <- 'https://raw.githubusercontent.com/ccodwg/Covid19Canada/master/timeseries_prov/vaccine_completion_timeseries_prov.csv'
com_vaccine = read.csv(urlvac)
com_vaccine$Date = as.Date(com_vaccine$date_vaccine_completed, format="%d-%m-%Y")

save(ABdata,file="/Users/dawo/data/ABdata.RData")
save(com_vaccine,file="/Users/dawo/data/ABvac.RData")

### loading data
load("/Users/dawo/data/ABdata.RData")
load("/Users/dawo/data/ABvac.RData")

#head(ABdata)
#select model data
cutoff_start <- "2021-10-21" #alberta
cutoff_end <- max(com_vaccine$Date)  # data included for model fitting
death_long <- ABdata %>% filter(Date >= cutoff_start & Date <= cutoff_end)
vac_long <- com_vaccine %>% filter (Date >= cutoff_start & Date <= cutoff_end & province=="Alberta")
cases_long <- ABdata %>% filter(Date >= cutoff_start & Date <= cutoff_end)
# display the table
death_long %>%
  head() %>%
  gt() %>%
  tab_options(table.width = "100%") %>%
  tab_style(style = table_style, 
            locations = cells_body()) %>% 
  opt_all_caps()

#######################################################
# projection
######################################################

### combining death, case, and vaccination into one dataframe
#Alberta <-  merge (death_long, vac_long , cases_long, by = "Date")

Alberta <- merge ( merge (death_long[,c( "Number.of.deaths", "Cumulative.number.of.deaths","Date" )], 
                       cases_long[,c( "Number.of.cases", "Cumulative.number.of.cases", "Date" ) ], 
                       by = "Date") ,
                vac_long[,c("cvaccine","cumulative_cvaccine","Date")],
                by = "Date")
names(Alberta)[5] <- "cumulative_cases"
names(Alberta)[3] <- "cumulative_deaths"

  ###### varying lags = 14 days ahead(tested)
###tested by the best fit model giving the largest R^2 -> by GK's function n_ahead = 21
  n_ahead = 21
  Alberta <- Alberta %>%
    # create lags by day
    tk_augment_lags(cumulative_deaths, .lags = 0:-max_lead, .names = "auto") 
  # fix names to remove minus sign
  names(Alberta) <- names(Alberta) %>% str_replace_all("lag-|lag", "lead")

  # use only case dates where we have complete future knowledge of deaths for all lead times.
  Alberta$led_deaths <- lead (Alberta$cumulative_deaths, n_ahead)
 # Alberta <- Alberta %>% mutate (cum_led_deaths = cumsum(led_deaths))
  all_data_raw <- na.omit (Alberta)
  
  #data sample
  all_data <- all_data_raw [1:(dim(all_data_raw)[1] - n_ahead),]
  
  complete_date <- all_data %>% complete(Date = seq.Date(min(Date), max(Date), by = "day")) 
  id <- which(all_data$Date %in% complete_date$Date)
  
  #### lag = n_ahead
  nfit = dim(all_data)[1]
  #new data set for prediction out of sample
  newdataraw <- Alberta %>% filter (Date > max (all_data$Date) & Date <= (max (all_data$Date) + n_ahead)) %>% 
    select (c("Date", "cumulative_cases", "cumulative_cvaccine"))
  
  new_complete_date <- newdataraw %>% complete(Date = seq.Date(min(Date), max(Date), by = "day")) 
  
  new_id <- which(newdataraw$Date %in% new_complete_date$Date)
  
    ####################
  ## Run a regression on lagged cases and date vs deaths
  # tvLM: time varying coefficient in linear regression model
  # 
  tvLMdf <- tvLM(led_deaths ~ cumulative_cases + cumulative_cvaccine , data = all_data)
  
  fit_tvreg <- tvLMdf$fitted 
# Create predicted deaths based on most recent case counts and vaccination counts
  cbind(fit_tvreg, all_data$led_deaths)
  
  predicted_tvreg <- forecast (tvLMdf, newdata = as.matrix(newdataraw[,2:3]), n.ahead = n_ahead)

  #deaths observations
  death_cum_obs <- all_data_raw$led_deaths
  
  #tvreg = fitted tvreg + predict tvreg
  tvreg <- c(fit_tvreg, predicted_tvreg)

 
  predict <- data.frame(tvreg, death_cum_obs)
  predict$Date <- all_data_raw $ Date
 
  ##################### ggplot in sample prediction############

  gg <- predict %>%
    pivot_longer(cols = where(is.numeric)) %>%
    filter(name %in% c("death_cum_obs", "tvreg")) %>%
    ggplot(aes(x = Date, y = value, color = name, linetype = name)) +
    geom_line() + ylim(range(c(tvreg, death_cum_obs))) +
    geom_vline(xintercept = max(all_data$Date), col = "black", linetype = 1) +
    labs(
      title = "Actual vs. Predicted Deaths for Alberta provincial data",
      subtitle = "lead death ~ cumulative cases + cumulative Vaccinations (lead days = 21)",
      x = "Date",
      y = "Count",
      caption ="Source from U of T data & www.alberta.ca: covid-19 data summary ")
    gg
 ##########ggplot out of sample prediction############
  pred_date = max(Alberta$Date) - n_ahead + 1
  tvLM_out <- tvLM(led_deaths ~ cumulative_cases + cumulative_cvaccine, data = Alberta)
  #new data sample prediction
  newdata_out <- Alberta %>% filter (Date >= pred_date) %>% select (c("cumulative_cases", "cumulative_cvaccine", "led_deaths"))
  tvLM_out_prediction <- forecast (tvLM_out, newdata = as.matrix(newdata_out[,1:2]), n.ahead = n_ahead)
  pred_tvreg <- c(fit_tvreg, predicted_tvreg, tvLM_out_prediction)
  #death observations
  death_obs <- Alberta$led_deaths
  prediction <- data.frame(pred_tvreg, death_obs)
  prediction$date <- Alberta$Date + n_ahead
    
  gg_out <- prediction %>%
    pivot_longer(cols = where(is.numeric)) %>%
    filter(name %in% c("pred_tvreg", "death_obs")) %>%
    ggplot(aes(date, value, color = name, linetype = name)) +
    geom_line() + ylim(range(c(pred_tvreg, death_obs_pred))) +
    geom_vline(xintercept = max(all_data$Date), color = "Black") +
    geom_vline(xintercept = pred_date + n_ahead, color = "Black") + 

    labs(
      title = "Actual vs. Predicted Deaths for Alberta provincial data",
      subtitle = "lead death ~ cumulative cases + cumulative Vaccinations (lead days = 21)",
      x = "Date",
      y = "Count",
      caption = paste("Source from U of T data & www.alberta.ca: covid-19 data summary ")
    )
  gg_out
  
    
    
    