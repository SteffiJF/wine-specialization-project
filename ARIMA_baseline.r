#!/usr/bin/env Rscript

#libraries:
library(forecast)
library(quantmod)
library(sarima)
library(tidyverse)
library(RcppRoll)
library(DataCombine)
library(roll)


#Reading data
df <- read.csv(paste0("C:/Users/Stephanie/OneDrive - NTNU/Skole/Master/final_agg_data_56.csv"), sep=';', header = TRUE, stringsAsFactors = FALSE)

#Empty matrices for Y-values
MY <- matrix(NA, nrow = 154, ncol = 56)
MY2 <- matrix(NA, nrow = 154, ncol = 56)
MY3 <- matrix(NA, nrow = 154, ncol = 56)

pred.length <- 12


#Analysis with lag=12

for (area in 1:56){
  print('Area')
  print(area)
  data <- df[area+1]
  ts <- ts(data, start=c(2007,1), end=c(2019,10), frequency=12)
  for (date in 100:142){
    print(date)
    
    #Shortened time series
    short.ts <- ts[-(date+1):-length(ts)]
    
    #Using box-cox transformation
    lambda <- BoxCox.lambda(short.ts)
    arima.model <- auto.arima(diff(short.ts, lag=12), lambda=lambda, seasonal = TRUE)
    pred <- forecast(arima.model,h=pred.length)
    prediction <- pred$mean
    
    #Last 12 dates
    last.data <- short.ts[(length(short.ts)-11):length(short.ts)]
    
    #Using last 12 values to find inverse diference
    pred.inv <- diffinv(prediction, lag=12, xi=last.data)
    pred.inv <- pred.inv[(length(pred.inv)-pred.length+1):length(pred.inv)]
    
    #Adding prediction to time series for analysis
    prediction <- as.vector(pred.inv)
    full.prediction <- c(short.ts,prediction)
    actual.data <- ts(data, start=1, end=154, frequency=1)
    
    #Calculating dfY-values
    rs <- shift(roll_sum(full.prediction,12),12)
    vY <- Delt(rs,k=12)
    #vY <- vY[-1]
    vY[is.na(vY)] <- 0
    for (i in 1:length(vY)){
      if (vY[i]>=0.3){vY[i]=2}
      else if (vY[i]<=-0.3){vY[i]=0}
      else {vY[i]=1}
    }
    MY[date,area] <- vY[date]
    
    #Calculating dfY2-values
    vY2 <- as.vector(diff(zoo(full.prediction), lag = 12, na.pad=TRUE))
    vY2 <- vY2/shift(roll_sd(vY2, 12, min_obs=12),-1)
    vY2 <- shift(rollmean(vY2,12,align='left', fill=NA),1)
    vY2[is.na(vY2)] <- 0    
    for (i in 1:length(vY2)){
      if (vY2[i]>=1){vY2[i]=2}
      else if (vY2[i]<=-1){vY2[i]=0}
      else {vY2[i]=1}
    }
    
    MY2[date,area] <- vY2[date]
    
    #Calculating dfY3-values
    pct.diff <- Delt(full.prediction,k=12)
    vY3 <- pct.diff/shift(roll_sd(pct.diff, 12, min_obs=12),-1)
    vY3[is.na(vY3)] <- 0
    for (i in 1:length(vY3)){
      if (vY3[i]>=2){vY3[i]=1}
      else if (vY3[i]<=-2){vY3[i]=-1}
      else {vY3[i]=0}
    }
    vY3 <- shift(roll_sum(vY3,12),12)
    vY3[is.na(vY3)] <- 0
    for (i in 1:length(vY3)){
      if (vY3[i]>=3){vY3[i]=2}
      else if (vY3[i]<=-3){vY3[i]=0}
      else {vY3[i]=1}
    }
    
    MY3[date,area] <- vY3[date]
    
    print(pred$model)
  }
}

forecasts <- as.data.frame(MY)
forecasts2 <- as.data.frame(MY2)
forecasts3 <- as.data.frame(MY3)


path <- "C:/Users/Stephanie/OneDrive - NTNU/Skole/Master/Sarima models\\"
name <- paste(path,"MYdiff",".csv", sep = "")
write.csv(forecasts,name)
name <- paste(path,"MY2diff",".csv", sep = "")
write.csv(forecasts2,name)
name <- paste(path,"MY3diff",".csv", sep = "")
write.csv(forecasts3,name)

