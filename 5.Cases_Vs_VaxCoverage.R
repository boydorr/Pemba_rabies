#'---
#'Author: KL
#'Title: Modelling rabies cases against vaccination coverage
#'---
#'
# Explore to see whether there is any relationship between coverage and cases
# Need to do a regression with number of cases as the response variable and vaccination coverage as the explanatory variable
# Dealing with count data, so need to use poisson or negative binomial approach
# Need to account for the serious zero inflation!

rm(list=ls())

require(AER)
require(MASS)
require(pscl)
require(lme4)
library ("RVAideMemoire")
library(reshape2)
library("glmmTMB")
library(lattice)
library(ggplot2)
library(forcats)
library(tidyverse)

# Import data
casecov <- read.csv("output/VillageData.csv")
#names(casecov)[names(casecov)=="Var2"] <- "months"
names(casecov)

# to check for overdispersion of the data
attach(casecov)
mean(animal_cases)
var(animal_cases)
var(animal_cases)/mean(animal_cases)
# data is overdispersed as the mean is different from the variance.
# Negative binomial suits the data
# For poisson distribution, mean and variance should be equal or almost the same

# explore data
table(casecov$animal_cases)
barplot(table(factor(casecov$animal_cases, 0:max(casecov$animal_cases))))
nrow(casecov)
hist(casecov$dogs_vacc)
table(casecov$dogs_vacc == 0)
hist(casecov$vc_est)
table(casecov$vc_est == 0)
table(casecov$vc_est == 0, casecov$dogs_vacc == 0) # use waning-adjusted VC?

plot(jitter(vc_waning_est) ~ jitter(dogs_vacc), data = casecov)
plot(jitter(animal_cases) ~ vc_waning_est, data = casecov)

#-- Categorise vaccination coverage and add an extra column of vaccination data with the above categories
#  as either zero,poor,medium or good--

casecov$Factor_Coverage<-""
casecov$Factor_Coverage[casecov$vc_waning_est==0.00] <-"Zero"
casecov$Factor_Coverage[(casecov$vc_waning_est> 0.00 & casecov$vc_waning_est < 0.30)] <-"Low"
casecov$Factor_Coverage[(casecov$vc_waning_est>= 0.30 & casecov$vc_waning_est < 0.70)] <-"Medium"
casecov$Factor_Coverage[casecov$vc_waning_est>= 0.70] <-"High"

casecov$Factor_Coverage<-factor(casecov$Factor_Coverage)
table(casecov$Factor_Coverage, casecov$animal_cases)

# Compute the mean and standard  deviation and error for the vaccination coverage
meancom<-function(x, k)
{# computing the means
  Mn<-tapply(x, list(k), mean)
  # computing the standard error
  sD<-tapply(x, list(k), sd)
  lD<-tapply(x, list(k), length)
  se<-sD/sqrt(lD)
  print("95% CI")
  Lower = Mn - qnorm(0.975)*se
  Upper = Mn + qnorm(0.975)*se
  All<-cbind(Mn, Lower, Upper)
  print(All)
}

allmtr<-as.data.frame(meancom(casecov$animal_cases, casecov$Factor_Coverage))
allmtr$Coverage<-rownames(allmtr)

#'*: Order to Set the levels*
allmtr$Coverage <- factor(allmtr$Coverage,
                          levels=c("Zero", "Low", "Medium", "High"))

ggplot(allmtr, aes(fct_rev(fct_reorder(Coverage,-Mn)), y=Mn)) +
  geom_bar(aes(x=Coverage, y=Mn), stat="identity", fill="skyblue", alpha=0.5) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.4, colour="orange", alpha=0.9, size=1.3)+
  xlab("Coverage") + ylab("Mean Coverage")

# making everything else low, then recategorise vaccination coverage
casecov$Factor_Coverage<-"Low"
casecov$Factor_Coverage[casecov$vc_waning_est<=0.30] <-"Low"
casecov$Factor_Coverage[(casecov$vc_waning_est> 0.30 & casecov$vc_waning_est < 0.70)] <-"Medium"
casecov$Factor_Coverage[casecov$vc_waning_est>= 0.70] <-"High"

# using the mean function above, recalculate
allmtr<-as.data.frame(meancom(casecov$animal_cases, casecov$Factor_Coverage))
allmtr$Coverage<-rownames(allmtr)

#' _: To set the the orders, to set the levels:_
allmtr$Coverage <- factor(allmtr$Coverage,
                          levels=c("Low", "Medium", "High"))

  ggplot(allmtr, aes(fct_rev(fct_reorder(Coverage, -Mn)),y=Mn)) +
  geom_bar(aes(x=Coverage, y=Mn), stat="identity", fill="skyblue", alpha=0.5) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.4, colour="orange", alpha=0.9, size=1.3)+
  xlab("Vaccination Coverage") + ylab("Mean Rabies Cases per month")

pdf("figures/Figure_3.meancases_per_Vaxcoverage.pdf", width = 9, height = 6)
allmtr %>%
  mutate(name = fct_reorder(Coverage, desc(Mn))) %>%
  ggplot( aes(x=Coverage, y=Mn)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.8, width=.6) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.4, colour="red", alpha=0.9, size=1.3) +
  xlab("Vaccination coverage") + 
  ylab("Mean rabies cases per month")+
  #theme_bw()+
  theme_classic()+
  theme(axis.text.x  = element_text(size = 16, face="plain"),
        axis.title.x = element_text(size=16, face="plain"),
        axis.text.y = element_text(size = 16, face="plain"),
        axis.title.y = element_text(size=16, face="plain"))

dev.off()


#' Potential figure:Lagged coverage in shehia (i.e. at month prior) vs. cases in the 
#' next month







#'----------------------------------------------------------------------------
#'*------- Run glms to see how cases  vary with vaccination coverage*

m1<-glm(animal_cases~vc_est,family=poisson, data= casecov)
summary(m1)
exp(coef(m1)) # getting back the coefficients
anova(m1,test="Chisq") # data overdispersed

# Trying a negative binomial and checking for model fits
m2<-glm.nb(animal_cases~vc_est, data= casecov) # also highly overdispersed as the residual degrees of freedom are greater the the null deviences
summary(m2)
summary(m1)

# Checking for model fit
1-pchisq(deviance(m2),df.residual(m2)) # the model also fits the data as P > 0.05

# Opt for the standard poison model
# Accounting for the zero inflation values in the model
m3 <- zeroinfl(animal_cases ~ vc_est, data= casecov)
summary(m3) # zero inflation model statistically significant except for the count model

# We expect cases to decrease with increase in coverage
# Use waning coverages and include random effects into the model

# consider cases and coverages only
m4 = glm.nb(animal_cases ~ vc_waning_est, data = casecov) 
summary(m4) # not significant 

#'----
# Create a full model with all the possible combination of interactions and use the drop1 function
# to simplify the model using LRT. 
# Based on output of drop1 I will delete the least significant variable and will then run a simpler model to carry on the model selection.
m5 = glmer.nb(animal_cases ~ vc_waning_est  + (1|shehia) + factor(year),data = casecov) 
summary(m5)  # sigificant relationship btn cases and vacc coverage <0.0001 except in 2016, 0.871,2019-0.78


m6 = glmer(animal_cases~vc_waning_est + factor(year) + (1|factor(year))  + (1|shehia) 
           + dog_pop_est,family = poisson, data = casecov) # factor()random and fixed effect for time
summary(m6)


m6.1 = glmer(animal_cases~vc_waning_est + (1|factor(year)) + (1|shehia), 
             family = poisson, data = casecov) 
summary(m6.1) # significant  Pvalue = 0.0001

anova(m6.1,test="Chi")

#--------
# Explore how vaccination coverages grouped into either low, medium or high affects rabies cases

# run the model
m7 = glmer.nb(animal_cases~Factor_Coverage + (1|shehia), data = casecov) 
summary(m7)

# Simpler model
m8 = glmer(animal_cases ~ Factor_Coverage + (1|shehia), family = poisson, data = casecov) 
summary(m8)

m9 = glmer(animal_cases ~ Factor_Coverage + (1|shehia) +dog_pop_est, family = poisson, data = casecov) 
summary(m9)


# M4 seems most intuitive - but with other terms used significance is no longer clear!
m10 = glmer(animal_cases ~ Factor_Coverage*factor(year) + (1|factor(year)) + (1|shehia), family = poisson, data = casecov) # random and fixed effect for time
summary(m10) 

# Removed months as a random effect
m11 = glmer(animal_cases ~ Factor_Coverage*factor(year) + (1|shehia), family = poisson, data = casecov)  
summary(m11) # still not significant 

# Look over larger timescales (quarters to reduce n zeros, or years) and groups of shehias

# Include offset (might improve fit because already scaled by population)
# e.g. + offset(log(ndogs))

# Change cases to binomial level
casecov$C = as.integer(casecov$animal_cases>0.5)
mB1 = glmer(C ~ Factor_Coverage + (1|shehia), family = binomial, data = casecov) # 
summary(mB1)

mB2 = glmer(C~ factor(year) + (1|shehia), family = binomial, data = casecov) # 
summary(mB2)
anova(mB2,alpha=0.05) # There is  good correlation btn cases and time...cases declining with time


#----------------------------------------------------------
# Estimating the dog population to date
# Estimating the human dog ratio
#------------------------------------------------------------
data2018 <- subset(casecov, casecov$year =="2018")
dat18 <- as.data.frame(data2018)
vacc18 <- sum(dat18$dogs_vacc)/sum(dat18$dog_pop_est)*100;vacc18 # 22.12 % from HDR...
dogs18 <-  sum(dat18$dog_pop_est);dogs18 # 2916 dog population for 2018

data2019 <- subset(casecov, casecov$year =="2019")
dat19 <- as.data.frame(data2019)
vacc19 <- sum(dat19$dogs_vacc)/sum(dat19$dog_pop_est)*100;vacc19 # 50.9 % from HDR...
dogs19 <-  sum(dat19$dog_pop_est);dogs19 # 3131,dog population for 2019

# Total human population of Pemba Island
humanpop <- 472958

# Vaccination coverage  for 2016 and 17
# Using the dog population estimated from the the dog to human ratio greatly
# underestimates the vaccination coverage
# Better to use coverages estimated by LFOs using their dog estimates
data2017 <- subset(casecov, casecov$year =="2017")
dat17 <- as.data.frame(data2017)
vacc17 <- sum(dat17$dogs_vacc)/sum(dat17$dog_pop_est)*100;vacc17 # 49.2

data2016 <- subset(casecov, casecov$year =="2016")
dat16 <- as.data.frame(data2016)
vacc16 <- sum(dat16$dogs_vacc)/sum(dat16$dog_pop_est)*100;vacc16 # 66.72

# Dog to guman ratio
HDR <- humanpop/dogs18;HDR # 1:61

# Calculate mean and range in vaccinated shehias
meanvacc <- mean(casecov$vc_waning_est[which(casecov$dogs_vacc>0)], na.rm=T);meanvacc
range <- range(casecov$vc_waning_est[which(casecov$dogs_vacc>0)], na.rm=T);range



