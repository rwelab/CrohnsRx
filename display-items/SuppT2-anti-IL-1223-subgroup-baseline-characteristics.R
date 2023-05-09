library(dplyr)
library(table1)     # table 1

setwd('UCSF/ipd-ma-cd2-2')

source('helper/subgroup.R') # subgroup.cohorts()

# load preprocessed data
drug_ordering <- read.csv('data/treatment_preferences-n10000.csv')

# IL12 > (TNFI, INTG)
prefer_il12 <- drug_ordering %>% 
  filter(p12_ohe == 1 & drug1 %in% c('il12'))

#------------------------------------------------------------------------------#

tbl1 <- prefer_il12 %>% 
  # re-factor categorical variables
  mutate(Sex = ifelse(Sex_Male == 1, 'Male', 'Female')) %>% 
  mutate(HxOfTNFi   = factor(HxOfTNFi,
                             levels=c(1,0),
                             labels=c('Yes','No'))) %>% 
  mutate(SteroidUse = factor(SteroidUse,
                             levels=c(1,0),
                             labels=c('Yes','No'))) %>% 
  mutate(ImmUse     = factor(ImmUse,
                             levels=c(1,0),
                             labels=c('Yes','No'))) %>% 
  mutate(Ileal      = factor(Ileal,
                             levels=c(1,0),
                             labels=c('Yes','No')))

tbl1 %>% glimpse()

#------------------------------------------------------------------------------#

# Format continuous covariates Mean (sd)
my.render.cont <- function(x){
  with(stats.apply.rounding(stats.default(x), digits=2), 
       c("", "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}

# Format binary covariates N (%)
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y, sprintf("%d (%0.0f %%)", FREQ, PCT))))
}

# Format binary covariates - include positive (1; Yes) only
my.render <- function(x,...) {
  y <- render.default(x,...)
  if(is.factor(x)) y[2] else y
}

table1(~ Age + BMI + CDAI_baseline + CRP + HxOfTNFi + 
         SteroidUse + ImmUse + Ileal | Sex, 
       data=tbl1,
       render = my.render, 
       render.continuous = my.render.cont, 
       render.categorical= my.render.cat,
       droplevels = T)

#------------------------------------------------------------------------------#

# copy-paste tabel to word document and modify

#------------------------------------------------------------------------------#