##----------------------------------------------------------
# Author:  Kennedy Lushasi
# This script produces a table of the annual vaccination coverage 
# for each district of Pemba region, with the mean, standard deviation, 
# median and range of each campaign
##---------------------------------------------------------

# import data
rm(list=ls())

library(rgdal)
require(tidyverse)
library(data.table)
library(dplyr)
options(stringsAsFactors=F) 

## Vaccination coverage data
vacc <- read.csv("output/vcYearPembaDist.csv", row.names = 1,header=F)
colnames(vacc) <- c("2010","2011", "2012","2013","2014","2015","2016","2017","2018","2019","2020")
setDT(vacc, keep.rownames = "District") # add column name for the districts
View(vacc)

vacc[,c("2015")] <- list(NULL)
colnames(vacc)
names(vacc)
vacc$`2010`<- vacc$`2010`*100
vacc$`2011`<- vacc$`2011`*100
vacc$`2012`<- vacc$`2012`*100
vacc$`2013` <- vacc$`2013`*100
vacc$`2014` <- vacc$`2014`*100
vacc$`2015` <- vacc$`2015`*100
vacc$`2016` <- vacc$`2016`*100
vacc$`2017` <- vacc$`2017`*100
vacc$`2018` <- vacc$`2018`*100
vacc$`2019` <- vacc$`2019`*100

# calculate the mean and Sd, range and median

# for 2011
m10 <- mean(vacc$`2010`);m10
s10 <- sd(vacc$`2010`);s10
r10 <- range(vacc$`2010`);r10
med <- median(vacc$`2010`);med

# for 2011
m11 <- mean(vacc$`2011`);m11
s11 <- sd(vacc$`2011`);s11
r11 <- range(vacc$`2011`);r11
med <- median(vacc$`2011`);med

# for 2012
m12 <- mean(vacc$`2012`);m12
s12 <- sd(vacc$`2012`);s12
r12 <- range(vacc$`2012`);r12
med12 <- median(vacc$`2012`);med12

# for 2013
m13 <- mean(vacc$`2013`);m13
s13 <- sd(vacc$`2013`);s13
r13 <- range(vacc$`2013`);r13
med13 <- median(vacc$`2013`);med13

# for 2014
m14 <- mean(vacc$`2014`);m14
s14 <- sd(vacc$`2014`);s14
r14 <- range(vacc$`2014`);r14
med14 <- median(vacc$`2014`);med14

# for 2016
m16 <- mean(vacc$`2016`); m16
s16 <- sd(vacc$`2016`); s16
r16 <- range(vacc$`2016`); r16
med16 <- median(vacc$`2016`); med16 

# for 2017
m17 <- mean(vacc$`2017`); m17
s17 <- sd(vacc$`2017`); s17
r17 <- range(vacc$`2017`); r17
med17 <- median(vacc$`2017`); med17 

# for 2018
m18 <- mean(vacc$`2018`); m18
s18 <- sd(vacc$`2018`); s18
r18 <- range(vacc$`2018`); r18
med18 <- median(vacc$`2018`); med18 

# for 2019
m19 <- mean(vacc$`2019`); m19
s19 <- sd(vacc$`2019`); s19
r19 <- range(vacc$`2019`); r19
med19 <- median(vacc$`2019`); med19 


# Dog population data 2020

Dogpop <- read.csv("output/dogPopMatPembaDist.csv",row.names = 1,header=F)
setDT(Dogpop, keep.rownames = "District") # add column name for the districts
View(Dogpop)

# Months are arranged from Jan 2010 to Dec 2020 corresponding to 132 months
# The last month on the column represents the dog populatiom for Dec 2020

totalDogs2020 <- sum(Dogpop$V133);totalDogs2020


