#---------------------------------------------------------
# Author: Kennedy Lushasi
# Clean and processes animal and human rabies contact tracing data for further analysis
#---------------------------------------------------------
rm(list=ls())
## load all libraries
library(tidyverse)
library(lubridate)

# ----- Import & process contact tracing data ------------------------------------------------------------
# Humans:
datestamp = "20220227040144" # "20211201071412" 
humans <- read.csv(paste("data/Tanzania_Human_Contact_Tracing_", datestamp, ".csv", sep=""), stringsAsFactors = FALSE)

# Animals:
datestamp = "20220210225548" # "20210506171615"  
animals <- read.csv(paste("data/Tanzania_Animal_Contact_Tracing_", datestamp, ".csv", sep=""), stringsAsFactors = FALSE)

#----- Process CT data ---------------------------------------------------
# Subset for pemba data only
study_regions <- c("Kusini Pemba", "Kaskazini Pemba")
humans <- humans[which(humans$Region %in% study_regions),]
animals <- animals[which(animals$Region %in% study_regions),]

#----- FIX DATES ---------------------------------------------------
startDate <- as.Date("2010-01-01")

humans <- humans %>%
  mutate( # Format dates
    Date.bitten = as.Date(Date.bitten, format="%Y-%m-%d"),
    Date.reported = as.Date(Date.reported, format="%Y-%m-%d"),
    PEP.1.Date = as.Date(PEP.1.Date, format="%Y-%m-%d")) %>%
  mutate(Date_proxy = as.Date(ifelse(!is.na(Date.bitten), Date.bitten, # DATE PROXY for missing dates (bitten or reported)
                                 ifelse(!is.na(Date.reported), Date.reported,
                                        ifelse(!is.na(PEP.1.Date), PEP.1.Date, NA))), origin="1970-01-01")) %>%
  mutate(Year.bitten = year(Date_proxy), # Create year bitten and month
         month_bitten = month(Date_proxy) + (Year.bitten-2010)*12) %>%
  filter(Date_proxy >= startDate)
nrow(humans) #

animals <- animals %>%
  mutate( # Format dates
    Symptoms.started = as.Date(Symptoms.started), # format="%d-%b-%Y"),
    Date.bitten = as.Date(Date.bitten), # format="%d-%b-%Y"),
    Year = year(Symptoms.started),
    month_symp = month(Symptoms.started) + (Year-2010)*12) %>%
  filter(Symptoms.started >= startDate | Date.bitten >= startDate)  %>%
  mutate(Loc_ID = paste(District, Ward, sep="_")) # Create location IDs
nrow(animals)

#-----  SELECT RABID ANIMAL DATA -------------------------------------------------
unknown <- subset(animals, Suspect=="Yes"| Suspect=="Unknown"| is.na(Symptoms.started)); dim(unknown) 
which(is.na(animals$Symptoms.started)) # SHOULD WE REMOVE NAs for symptoms started?
rabid <- subset(animals, Suspect=="Yes"); dim(rabid)

# ---------------------------------------------------------------------------
# Write a csv files
 write.csv(animals, "Output/animal_data.csv", row.names=FALSE) # Contact tracing data - bitten and biters
 write.csv(humans, "Output/human_bites.csv", row.names=FALSE) # For exposures
 write.csv(rabid, "Output/animal_cases.csv", row.names=FALSE) # For rabid animals
# # ---------------------------------------------------------------------------
