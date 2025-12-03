#load data 
source("greenhouseData.R")

#create figure 1: time series plot of CO_2, Humidity, and Soil Temp
source("greenhouseTSplot.R")

#create figure 2: ACFs and CCFs
source("greenhouseACFCCF.R")

#create figures 3 and 4 (prelim): variance and covariance traces
source("greenhouseTraces.R")

#create table 1: estimates against T
source("greenhouseTTable.R")

#choose T and create figure 4 (final)
source("greenhouseChooseT.R")