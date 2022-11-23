#---------------------------------------------------------
# Author: Kennedy Lushasi
# Explores exposure and rabies data - compare contact tracing vs hospital records, PEP use & adherance
#---------------------------------------------------------
rm(list=ls())

## load all libraries
library(lubridate)
library(tidyverse)
library(Hmisc) # binomial probability of death if no PEP
options(stringsAsFactors=F)

#---------------------------------------------------------
#  Set start dates
startDate <- as.Date("2010-01-01")

## Read in data - human exposures and animal cases
humans <- read.csv("Output/human_bites.csv",row.names = NULL) # People bitten by both healthy & rabid animals, irrespective of health seeking
cases <- read.csv("Output/animal_cases.csv", na.strings=c("", "NA")) # only probable rabid animal cases
animals <- read.csv("Output/animal_data.csv", na.strings=c("", "NA")) # All suspect and non suspect animal data
humanpop <- 472958 # from 2012 census

#---- Process data -------------------------------------------------------------

# Identify if received PEP
humans$no.PEP <- (humans$If.no.PEP..why.not.Sought.but.not.available=="true" |
                    humans$If.no.PEP..why.not.Not.aware=="true" |
                    humans$If.no.PEP..why.not.Too.expensive=="true" |
                    humans$If.no.PEP..why.not.Not.advised.by.health.worker=="true" |
                    humans$If.no.PEP..why.not.Not.advised.by.dog.owner=="true" |
                    humans$If.no.PEP..why.not.Didn.t.think.was.needed=="true" |
                    humans$If.no.PEP..why.not.Currently.seeking....outcome.unknown=="true")

# or got late PEP or were searching
humans$late.PEP <- (humans$PEP.Status.Raising.funds=="true"| humans$PEP.Status.Seeking=="true")

# Create columns stating pre-/post-outbreak
table(humans$Year.bitten) # None in 2015
humans$outbreak <- ifelse(humans$Year.bitten < 2015, "Pre-outbreak",
                          ifelse(humans$Year.bitten > 2015, "Post-outbreak", NA))
table(humans$outbreak, useNA="always")

# Create subset of exposures (known rabid)
exposures <- subset(humans, Rabid == "Yes"); nrow(exposures)

#---------- Calculate estimates of incidence per year -----------------------------------
# Create table of exposures, deaths and patients (that went for PEP) for bites by healthy animals by YEAR
bites_y <- humans %>%
  mutate(Year = factor(Year.bitten, levels = 2010:2021)) %>%
  group_by(Year, .drop=FALSE)  %>%
  dplyr::summarize(
    exposures = length(which(Rabid=="Yes")), # Number of likely exposures
    poss_exposures = length(which(Rabid == "Unknown")), # Number of possible exposures - removed Rabid=="Yes" | as other "poss" columns exclude this
    healthy = length(which(Rabid=="No" & PEP.1=="true")), # Number of Healthy bites that received PEP
    deaths = length(which(Patient.outcome=="Died")), # Number of rabid exposures that died
    exp_inc = exposures*100000/humanpop, # Exposure incidence
    death_inc = deaths*100000/humanpop, # Death incidence
    PEP1 = length(which(PEP.1=="true")), # Number of patients (rabid & healthy) that received PEP dose 1
    PEP2 = length(which(PEP.2=="true")), # Number of patients (rabid & healthy) that received PEP dose 2
    PEP3 = length(which(PEP.3=="true")), # Number of patients (rabid & healthy) that received PEP dose 3
    noPEP_rabid = length(which(Rabid=="Yes" & no.PEP == TRUE)), # Number of rabid exposures that did not receive PEP
    noPEP_rabid_poss = length(which(Rabid=="Unknown" & no.PEP == TRUE)), # Number of possible exposures that did not receive PEP
    latePEP_rabid = length(which(Rabid=="Yes" & late.PEP == TRUE)), # Number of rabid exposures that received PEP late
    latePEP_rabid_poss = length(which(Rabid=="Unknown" & late.PEP == TRUE)), # Number of possible exposures that received PEP late
    rabid_noPEPsought = length(which(Rabid=="Yes" & PEP.Status.Not.sought=="true")), # Number of rabid exposures that did not seek PEP
    possrabid_noPEPsought = length(which(Rabid=="Unknown" & PEP.Status.Not.sought=="true")), # Number of rabid exposures that did not seek PEP
    rabid_PEPcomplete = length(which(Rabid=="Yes" & PEP.Status.Completed=="true")), # Number of rabid exposures that completed their PEP treatment
    rabid_PEPIncomplete =length(which(Rabid=="Yes" & PEP.Status.Incomplete=="true"))) # Number of rabid exposures that did not complete their PEP treatment
# View table
View(bites_y)

# Create alternative table that groups by pre-outbreak and post-outbreak
## CAN WE INTEGRATE THE SUMMARY VARIABLES BELOW INTO THE ABOVE TABLE AND THEN JUST AGGREGATE PRE/POST 2015 MORE 
bites_outbreak_gr <- humans %>%
  mutate(outbreak = factor(outbreak, levels = c("Pre-outbreak", "Post-outbreak"))) %>%
  group_by(outbreak)  %>%
  dplyr::summarize(
    exposures = length(which(Rabid=="Yes")), # Number of likely exposures
    poss_exposures = length(which(Rabid == "Unknown")), # Number of possible exposures - removed Rabid=="Yes" | as other "poss" columns exclude this
    healthy = length(which(Rabid=="No" & PEP.1=="true")), # Number of Healthy bites that received PEP
    deaths = length(which(Patient.outcome=="Died")), # Number of rabid exposures that died
    exp_inc = round(exposures*100000/humanpop, digits=3), # Exposure incidence
    death_inc = round(deaths*100000/humanpop, digits=3), # Death incidence
    PEP1 = length(which(PEP.1=="true")), # Number of patients (rabid & healthy) that received PEP dose 1
    PEP2 = length(which(PEP.2=="true")), # dose 2
    PEP3 = length(which(PEP.3=="true")), # dose 3
    noPEP_rabid = length(which(Rabid=="Yes" & no.PEP == TRUE)), # Number of rabid exposures that did not receive PEP
    noPEP_rabid_poss = length(which(Rabid=="Unknown" & no.PEP == TRUE)), # Number of possible exposures that did not receive PEP
    noPEP_healthy = length(which(Rabid=="No" & no.PEP == TRUE)), # Number of Healthy bites that did not receive PEP
    latePEP_rabid = length(which(Rabid=="Yes" & late.PEP == TRUE)), # Number of rabid exposures that received PEP late
    latePEP_rabid_poss = length(which(Rabid=="Unknown" & late.PEP == TRUE)), # Number of possible exposures that received PEP late
    latePEP_healthy = length(which(Rabid=="No" & late.PEP == TRUE)), # Number of Healthy bites that received PEP late
    noPEPsought_rabid = length(which(Rabid=="Yes" & PEP.Status.Not.sought=="true")), # Number of rabid exposures that did not seek PEP
    noPEPsought_possrabid = length(which(Rabid=="Unknown" & PEP.Status.Not.sought=="true")), # Number of possible exposures that did not seek PEP
    noPEPsought_healthy = length(which(Rabid=="No" & PEP.Status.Not.sought=="true")), # Number of Healthy bites that did not seek PEP
    PEPcomplete_rabid = length(which(Rabid=="Yes" & PEP.Status.Completed=="true")), # Number of rabid exposures that completed their PEP treatment
    PEPcomplete_possrabid = length(which(Rabid=="Unknown" & PEP.Status.Completed=="true")), # Number of possible exposures that completed their PEP treatment
    PEPcomplete_healthy = length(which(Rabid=="No" & PEP.Status.Completed=="true")), # Number of Healthy bites that completed their PEP treatment
    PEPIncomplete_rabid =length(which(Rabid=="Yes" & PEP.Status.Incomplete=="true")), # Number of rabid exposures that did not complete their PEP treatment
    PEPIncomplete_possrabid =length(which(Rabid=="Unknown" & PEP.Status.Incomplete=="true")), # Number of possible exposures that did not complete their PEP treatment
    PEPIncomplete_healthy =length(which(Rabid=="No" & PEP.Status.Incomplete=="true"))) # Number of Healthy bites that did not complete their PEP treatment
# View table
View(bites_outbreak_gr)
# Transform to long dataframe
bites_outbreak_gr_t <- t(bites_outbreak_gr)
colnames(bites_outbreak_gr_t) <- bites_outbreak_gr_t[1,]
bites_outbreak_gr_t <- bites_outbreak_gr_t[-1,]
# Save as output
write.csv(bites_outbreak_gr_t, "Output/Table_1.csv")


# Examine PEP access and deaths
n <- sum(bites_y$exposures); n # Total rabid exposures
n2 <- sum(bites_y$poss_exposures, bites_y$exposures); n2 # Total possible exposures
deaths <- subset(humans, Patient.outcome == "Died");
nDeaths <- nrow(deaths); nDeaths; sum(bites_y$deaths) # 6 deaths
noPEP <- sum(bites_y$noPEP_rabid); noPEP2 <- sum(bites_y$noPEP_rabid_poss) # 60 rabid/1 possible did not complete timely PEP

# THINK SOME MISUNDERSTANDING OF THE CALCULATIONS HERE - not sure who corrected but needs revisiting and cleaning up
pPEP <- (n-noPEP)/n; round(pPEP, 2); # probability of PEP for definite  exposures
pPEP2 <- (n2-noPEP2-noPEP)/n2; round(pPEP2, 2) # 0.75-0.78 probability of getting good PEP? (removed "not" - this calculated probability of receiving PEP)
# RS: In the line above, I changed (n2-noPEP)/n2 to (n2-noPEP2)/n2 - you were subtracting the wrong number

# Deaths
deaths$Date.bitten; deaths$When.died;
d1 <- which(deaths$Year.bitten == 2010)
d2 <- which(deaths$Year.bitten == 2017)
as.Date(deaths$When.died)-as.Date(deaths$Date.bitten) # 1 person died on same day as bite - is this correct? KENNEDY PLEASE CHECK!
# RS: This is easier to see with formatting
message("Head: ", paste0(deaths$Bite.site.Head...neck, collapse=", ")); 
message("Arms: ", paste0(deaths$Bite.site.Arms, collapse=", ")); 
message("Hands: ", paste0(deaths$Bite.site.Hands, collapse=", ")); 
message("Trunk: ", paste0(deaths$Bite.site.Trunk, collapse=", ")); 
message("Legs: ", paste0(deaths$Bite.site.Legs, collapse=", ")); 
message("Feet: ", paste0(deaths$Bite.site.Feet, collapse=", "))
message("PEP 1: ", paste0(deaths$PEP.1.Date, collapse=", ")); 
message("PEP 2: ", paste0(deaths$PEP.2.Date, collapse=", ")); 
message("PEP 3: ", paste0(deaths$PEP.3.Date, collapse=", "))
deaths$PEP.Status.Not.sought[d1]; deaths$PEP.Status.Completed[d1]; deaths$PEP.Status.Incomplete[d1]

# Examine proportion of patients due to healthy vs rabid dogs
bites_y$Year <- as.numeric(levels(bites_y$Year))[bites_y$Year]
p1 <- which(bites_y$Year < 2015)
p2 <- which(bites_y$Year > 2015)

# Note that rabid bite victims that did not seek care are NOT patients!
total_patients_p1 <- sum(bites_y$healthy[p1]) + sum(bites_y$exposures[p1]) + sum(bites_y$poss_exposures[p1]) - sum(bites_y$rabid_noPEPsought[p1]); total_patients_p1
total_patients_p2 <- sum(bites_y$healthy[p2]) + sum(bites_y$exposures[p2]) + sum(bites_y$poss_exposures[p2]) - sum(bites_y$rabid_noPEPsought[p2]); total_patients_p2
pExposures <- (sum(bites_y$exposures[p1]) - sum(bites_y$rabid_noPEPsought[p1]))/total_patients_p1; pExposures
# RS: Below, why do you subtract rabid exposures from possible (unknown) exposures? I've added a new column in bites_y for possible exposures (unknown) that did not seek PEP
# The exposures and possible exposures need combining - they have been seperated here which doesn't make logical sense - can we fix!
pExposures_uk <- (sum(bites_y$poss_exposures[p1]) - sum(bites_y$possrabid_noPEPsought[p1]))/total_patients_p1; pExposures_uk + pExposures
pExposures2 <- (sum(bites_y$exposures[p2]) - sum(bites_y$rabid_noPEPsought[p2]))/total_patients_p2; pExposures2

# Incidence of bite by healthy dogs (includes possible exposures that add uncertainty)
mean(bites_y$healthy[p1]); mean(bites_y$healthy[p1] + bites_y$poss_exposures[p1]) # bites_y$unknowns[p1])
mean(bites_y$healthy[p1])*100000/humanpop; mean(bites_y$healthy[p1] + bites_y$poss_exposures[p1])*100000/humanpop #bites_y$unknowns[p1])*100000/humanpop
mean(bites_y$healthy[p2]); mean(bites_y$healthy[p2])*100000/humanpop

# Numbers of patients that likely did not get PEP
pRabies <- 0.17 # from Changalucha et al probability
nDeaths/pRabies # would expect ~35 to not obtain PEP
pDeath_range <- binconf(sum(bites_y$exposures[p1])*pRabies, sum(bites_y$exposures[p1]))
pDeath_range_uk <- binconf(sum(bites_y$poss_exposures[p1])*pRabies, sum(bites_y$poss_exposures[p1])) # This now correctly only uses unknowns, previously included rabid exposures
pDeath_range; pDeath_range_uk
# RS: you don't add these to the bites_y table you created

sum(bites_y$noPEP_rabid[p1]) 
sum(bites_y$deaths[p1])/c(pDeath_range,pDeath_range_uk) # expect this many people to not seek treatment (prior to 2015)!
sum(bites_y$rabid_noPEPsought[p2]) # during the 2016 outbreak all sought PEP
# RS: surely this gives 39 that did not obtain ANY PEP? To get those that did not get adequate PEP, surely they need to have received at least 1 PEP dose?
sum(bites_y$noPEP_rabid[p2]) # but 39 didn't obtain (adequate) PEP - SO WHAT DID THEY GET? either late or incomplete PEP

#' _Exposures who did not not get PEP with reasons_
n1 <- sum(bites_y$exposures[p1]);n1
n2 <- sum(bites_y$exposures[p2]);n2

pre <- which(exposures$Year.bitten<2015)
post <- which(exposures$Year.bitten>2015)
  
shortages1 <- table(exposures$If.no.PEP..why.not.Sought.but.not.available[pre])
shortages2 <- table(exposures$If.no.PEP..why.not.Sought.but.not.available[post])
round(shortages1 *100/n1, 2)
round(shortages2*100/n2, 2)
# No pep shortage before the outbreak, but later 0.56% ran into shortages

unaware1 <- table(exposures$If.no.PEP..why.not.Not.aware[pre])
unaware2 <- table(exposures$If.no.PEP..why.not.Not.aware[post])
round(unaware1*100/n1,2)
round(unaware2*100/n2,2)
# 7.94% of exposures were not aware of PEP before outbreak, but after the outbreak it decreased to 1.67%

expense1 <- table(exposures$If.no.PEP..why.not.Too.expensive[pre])
expense2 <- table(exposures$If.no.PEP..why.not.Too.expensive[post])
round(expense1*100/n1,2)
round(expense2*100/n2,2)
# Before the outbreak cost was not a problem as PEP was provided at no cost, but during the outbreak, 0.56% failed to get PEP due to cost

illadvised1 <- table(exposures$If.no.PEP..why.not.Not.advised.by.health.worker[pre])
illadvised2 <- table(exposures$If.no.PEP..why.not.Not.advised.by.health.worker[post])
round(illadvised1*100/n1,2)
round(illadvised2*100/n2,2)
# 19.05% of exposures were not advised by HW before the outbreak vs 2.78 during the outbreak

dogowner1 <- table(exposures$If.no.PEP..why.not.Not.advised.by.dog.owner[pre])
dogowner2 <- table(exposures$If.no.PEP..why.not.Not.advised.by.dog.owner[post])
round(dogowner1*100/n1,2)
round(dogowner2*100/n2,2)
# During endemic period, no exposures failed to get PEP due not being advised by dog owners, but during the outbreak 0.56 did not get pep as were not advised by dog owners

#### IS THIS NOT THE SAME AS BEING NOT AWARE!????
table(exposures$If.no.PEP..why.not.Didn.t.think.was.needed[pre])
round(table(exposures$If.no.PEP..why.not.Didn.t.think.was.needed[pre])*100/n1,2)
table(exposures$If.no.PEP..why.not.Didn.t.think.was.needed[post])
round(table(exposures$If.no.PEP..why.not.Didn.t.think.was.needed[post])*100/n2,2)
# 12.7% did not think PEP was need during endemic period Vs 2.78 during the outbreak

seeking1 <- table(exposures$If.no.PEP..why.not.Currently.seeking....outcome.unknown[pre])
seeking2 <- table(exposures$If.no.PEP..why.not.Currently.seeking....outcome.unknown[post])
round(seeking1*100/n1,2)
round(seeking2*100/n2,2)
# 14.44% were still searching for pep at the time of investigation during the outbreak, and none during endemic period

table("PEP complete before 2015:"=exposures$PEP.Status.Complete[pre])
table("PEP incomplete before 2015:"=exposures$PEP.Status.Incomplete[pre])
table("PEP complete after 2015:"=exposures$PEP.Status.Complete[post])
table("PEP not sought after 2015:"=exposures$PEP.Status.Not.sought[post])
table("PEP late after 2015:"=exposures$late.PEP[post])

#----- RS: Create table of individual human deaths, ages, and reason for not obtaining PEP -----

# Create table of required columns
deaths_table <- data.frame("Year"=substr(deaths$When.died, 1, 4), # year(as.Date(deaths$When.died)) # better way to pull out year? (as.Date)
                           "Age"=deaths$Age..in.years.,
                           "Head"=deaths$Bite.site.Head...neck,
                           "Arms"=deaths$Bite.site.Arms,
                           "Hands"=deaths$Bite.site.Hands,
                           "Trunk"=deaths$Bite.site.Trunk,
                           "Legs"=deaths$Bite.site.Legs,
                           "Feet"=deaths$Bite.site.Feet,
                           "Severe_broken_bones"=deaths$Bite.details.Severe..broken.bones.,
                           "Severe_hospitalisation"=deaths$Bite.details.Severe..hospitalization.,
                           "Large_wound"=deaths$Bite.details.Large.wound.s.,
                           "Minor_wound"=deaths$Bite.details.Minor.wound.s.,
                           "Scratch"=deaths$Bite.details.Scratch,
                           "Currently_seeking"=deaths$If.no.PEP..why.not.Currently.seeking....outcome.unknown,
                           "Thought_unecessary"=deaths$If.no.PEP..why.not.Didn.t.think.was.needed,
                           "Not_advised_by_owner"=deaths$If.no.PEP..why.not.Not.advised.by.dog.owner,
                           "Not_advised_by_HW"=deaths$If.no.PEP..why.not.Not.advised.by.health.worker,
                           "Not_aware"=deaths$If.no.PEP..why.not.Not.aware,
                           "Not_available"=deaths$If.no.PEP..why.not.Sought.but.not.available,
                           "Too_expensive"=deaths$If.no.PEP..why.not.Too.expensive)

# Transform "true" to the column heading
deaths_table$Head <- ifelse(deaths_table$Head=="true", "Head", NA)
deaths_table$Arms <- ifelse(deaths_table$Arms=="true", "Arm", NA)
deaths_table$Hands <- ifelse(deaths_table$Hands=="true", "Hand", NA)
deaths_table$Trunk <- ifelse(deaths_table$Trunk=="true", "Trunk", NA)
deaths_table$Legs <- ifelse(deaths_table$Legs=="true", "Leg", NA)
deaths_table$Feet <- ifelse(deaths_table$Feet=="true", "Foot", NA)
deaths_table$Severe_broken_bones <- ifelse(deaths_table$Severe_broken_bones=="true", "Severe (broken bones)", NA)
deaths_table$Severe_hospitalisation <- ifelse(deaths_table$Severe_hospitalisation=="true", "Severe (hospitalisation)", NA)
deaths_table$Large_wound <- ifelse(deaths_table$Large_wound=="true", "Large wound", NA)
deaths_table$Minor_wound <- ifelse(deaths_table$Minor_wound=="true", "Minor wound", NA)
deaths_table$Scratch <- ifelse(deaths_table$Scratch=="true", "Scratch", NA)
deaths_table$Currently_seeking <- ifelse(deaths_table$Currently_seeking=="true", "Currently seeking", NA)
deaths_table$Thought_unecessary <- ifelse(deaths_table$Thought_unecessary=="true", "Thought unecessary", NA)
deaths_table$Not_advised_by_owner <- ifelse(deaths_table$Not_advised_by_owner=="true", "Not advised by dog owner", NA)
deaths_table$Not_advised_by_HW <- ifelse(deaths_table$Not_advised_by_HW=="true", "Not advised by Healthworker", NA)
deaths_table$Not_aware <- ifelse(deaths_table$Not_aware=="true", "Not aware", NA)
deaths_table$Not_available <- ifelse(deaths_table$Not_available=="true", "Not available", NA)
deaths_table$Too_expensive <- ifelse(deaths_table$Too_expensive=="true", "Too expensive", NA)

# Paste results together in new column, ignoring NAs
deaths_table <- deaths_table %>%
  unite("Bite_location", Head:Feet, sep=", ", remove = F, na.rm = T)
deaths_table <- deaths_table %>%
  unite("Bite_injury", Severe_broken_bones:Scratch, sep=", ", remove = F, na.rm = T)
deaths_table <- deaths_table %>%
  unite("Reason_for_no_PEP", Currently_seeking:Too_expensive, sep=", ", remove = F, na.rm = T)

# View table to check
View(deaths_table)

# Remove now surplus columns, and arrange by year then age
deaths_table <- deaths_table %>%
  dplyr::select(Year, Age, Bite_location, Bite_injury, Reason_for_no_PEP) %>%
  arrange(Year, Age)

# Save output
write.csv(deaths_table, "Output/Table_2.csv", row.names = F)

#-----Disease incidence in dogs per year--------

# Cases per year and between different species
rabid <- subset(animals, Suspect=="Yes") # only rabid animals
table(Year = rabid$Year, Species=rabid$Species) # does not list all years (missing 2015, 2019-2021 because no cases)
table(rabid$Species)*100/nrow(rabid) # proportion of each rabid species

dogpopulation <- 3156 # dog population has been adjusted from 4000 to the current estimates
cases_y <- rabid %>%
  mutate(Year = factor(Year, levels = 2010:2020)) %>%
  group_by(Year, .drop=FALSE)  %>%
  dplyr::summarize(n=n(),
                   dogs = length(which(Species=="Domestic dog")),
                   cows = length(which(Species=="Livestock: Cow")),
                   goats = length(which(Species=="Livestock: Goat")),
                   confirmed = length(na.omit(Lateral.flow.test=="Positive")),
                   case_inc = dogs*100/dogpopulation,# disease incidence in all animals - actually tractable as %
                   posivity_inc = confirmed*100/dogpopulation) # sample positivity incidence - actually tractable as %
cases_y
sum(cases_y$n); sum(cases_y$dogs);

# Look at dog rabies cases per year and case detection extrapolation
cases_y$Year <- as.numeric(levels(cases_y$Year))[cases_y$Year]
pDetect <- c(0.5383871, 0.6935484) # detection probability of cases in Pemba before 2015 and after # 
pBite <- 0.38 # mean bites per rabid dog
pBite_range <- binconf(sum(bites_y$exposures[p1])*pBite, sum(bites_y$exposures[p1])) # expected proportions
pBite_range2 <- binconf(sum(bites_y$exposures[p2])*pBite, sum(bites_y$exposures[p2])) # expected proportions
p1_exposures <- c(sum(bites_y$exposures[p1]), sum(bites_y$exposures[p1])+sum(bites_y$poss_exposures[p1]))
p2_exposures <- sum(bites_y$exposures[p2])

rabid_dogs <- c(sum(cases_y$dogs[p1]), sum(cases_y$dogs[p2], na.rm=TRUE))
est_rabid_dogs <- rabid_dogs/pDetect; est_rabid_dogs # Adjust rabid dog estimates given case detection
est_exposures_p1 <- est_rabid_dogs[1] * pBite_range; est_exposures_p1 # estimated exposures based on dog behaviour
est_exposures_p2 <- est_rabid_dogs[2] * pBite_range2; est_exposures_p2 # estimated exposures based on dog behaviour

# Not sure that the following lines - 325-331 are needed?
pPEP_p1 <- c(p1_exposures-sum(bites_y$rabid_noPEPsought[p1]),
             p1_exposures-sum(bites_y$noPEP_rabid[p1]),
             p1_exposures-sum(bites_y$rabid_PEPcomplete[p1]))/p1_exposures # This last row is incorrect!
pPEP_p1
PEP_incomplete <- c(p2_exposures-sum(bites_y$rabid_PEPcomplete[p2]), sum(bites_y$noPEP_rabid[p2]))
pPEP_p2 <- (p2_exposures-PEP_incomplete)/p2_exposures # pPEP after 2015
pPEP_p2

# Ratio of exposures per rabid dog (based on adjusted rabid dog cases according to case detection)
p1_exposures/est_rabid_dogs[1]
p2_exposures/est_rabid_dogs[2]

p2_exposures
#----- Case confirmation --------
# Number and % of confirmed animal cases
confirmed = sum(na.omit(cases$Lateral.flow.test=="Positive")) # are these LFT positive? There were laboratory confirmed by the vets department in Pemba
confirmed*100/nrow(cases) # report the proportion of cases confirmed

# Check if case confirmation (i.e. cases confirmed with the RDT) varies by species.
results_tab <- as.data.frame(table(results = cases$Lateral.flow.test, spp = cases$Species, useNA = "always"))
subset(results_tab, spp == "Domestic dog")$Freq
Total_dogs <- sum(subset(results_tab, spp == "Domestic dog")$Freq); Total_dogs
dogs_confirmed <- sum(na.omit(cases$Lateral.flow.test=="Positive" & cases$Species=="Domestic dog"))
Prop_dogs_Confirmed <- dogs_confirmed/Total_dogs*100; Prop_dogs_Confirmed # 18.6%
dogs_confirmed/confirmed

# Livestock confirmed
Livestock_confirmed <- confirmed - dogs_confirmed; Livestock_confirmed
cows <- sum(subset(results_tab, spp == "Livestock: Cow")$Freq)
goats <- sum(subset(results_tab, spp == "Livestock: Goat")$Freq)
Total_livestock <- sum(cows, goats)
prop_livestock_confirmed <- Livestock_confirmed/Total_livestock*100
prop_livestock_confirmed #23.1 %


#------------------------------------------------------------------
#----- Cases/ bites/ per rabid dog -----------
# How many human exposures per rabid dogs on average
n.exposures <- sum(bites_y$exposures); n.exposures
n.cases <- sum(cases_y$n) ; n.cases  
n.exposures/n.cases # 1.3 - BUT THIS IS NOT ADJUSTED FOR CASE DETECTION  - think this can be deleted!
length(unique(exposures$Biter.ID)) # How many dogs are responsible for all the bites?
n.exposures/length(unique(exposures$Biter.ID)) # Biting dogs bite ~how many people each? 
which(is.na(exposures$Biter.ID)) # No exposed persons with NA Biter.ID 

# Function that takes the rabid animals and the bite victims and calculates numbers bitten per rabid animal
bites <- function(biters, bite_victims){
  rmax = max(biters, na.rm=TRUE) + 1 # Find the range of the IDs of the bite victims
  bitten_index = hist(bite_victims, breaks = -1:rmax, plot=F)$counts[-1] # aggregate bite victims by biter ID
  NB = rep(NA, length(biters)) # create a vector to store the numbers bitten (NB)
  for (i in 1:length(biters)){
    NB[i] = bitten_index[biters[i]] # match the numbers bitten to each biter
  }
  NB
}

# Alternative function to check biting!
bites2 <- function(biters, bite_victims){
  bites_table = table(bite_victims) 
  NB = rep(0, length(biters))
  index = match(names(bites_table), biters)
  NB[index] =bites_table
  NB
}

n.exposed <- bites(cases$ID, exposures$Biter.ID) # number of people exposed per rabid ANIMAL (because cases are ALL spp!)
test = bites2(cases$ID, exposures$Biter.ID)
sum(n.exposed); sum(test)


# bite_stats <- 
#   biting_animals %>%
#   filter(Biter.ID > 0 & !is.na(Biter.ID)) %>%
#   group_by(Biter.ID) %>%
#   dplyr::summarise(Dogs.bitten = sum(Species %in% "Domestic dog"), 
#                    Animals.bitten = n())

##################################################################
n.exposed <- bites(cases$ID, exposures$Biter.ID) # number of people exposed per rabid ANIMAL (because cases are ALL spp!)
n.bitten <- bites2(animals$ID, !is.na(animals$Biter.ID)) # Number of dogs BITTEN per RABID DOG

max(animals$Biter.ID)
max(animals$ID)
test1 = hist(animals$Biter.ID, breaks = -1:607, plot=F)$counts[-1] 
test2 = rep(NA, length(rabid$ID))

n.rabid <- bites(rabid$ID, rabid$Biter.ID) # Number of dogs INFECTED (i.e. secondary cases) per RABID DOG

hist(n.exposed, -1:15) # THE DISTRIBUTION OF EXPOSED PEOPLE! 
hist(n.bitten, -1:15) # THE DISTRIBUTION OF BITTEN ANIMALS! 
hist(n.rabid, -1:15) # THE DISTRIBUTION OF SECONDARY CASES (RABID ANIMALS) - I think this will be quite biased but we should look at the transmission trees!

# Export animal cases (with IDs) and persons bitten:
rabid_exposed = data.frame(id_case = rabid$ID, species = rabid$Species, exposed = n.exposed)
write.csv(rabid_exposed, "Output/ID_exposed.csv", row.names = FALSE) # to be used for Fig 4




