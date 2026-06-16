#load data 
#library(tidyverse)
greenhouse2 <- read_csv("GreenhouseEnvironmentTimeSeries.csv")

#calculate NeweyWest
library(sandwich)

greenhouseNW <- as.matrix(greenhouse2[,-1])

#NW gives the covariance of betahat
#Multiply by X'X, here, n, to get covariance of sum Gamma(j)
#This is not Gamma(0)
NeweyWest( lm(greenhouseNW ~ 1)) * dim(greenhouseNW)[1]

