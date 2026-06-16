#load data 
source("greenhouseData.R")

#create figure 1: time series plot of CO_2, Humidity, and Soil Temp
source("greenhouseTSplot.R")

#create scatterplot matrix---not included in paper
source("greenhouseScatterplotMat.R")

#create figure 2: ACFs and CCFs
source("greenhouseACFCCF.R")

#create figures 3 and 4 (prelim): variance and covariance traces
source("greenhouseTraces.R")

#create table 1: estimates against T
source("greenhouseTTable.R")

#choose T and create figure 4 (final)
source("greenhouseChooseT.R")

#create table 2 (including calculation of Newey-West covariance matrix)
source("greenhouseCovTable.R")

#create table 3
source("greenhouseTable3.R")


