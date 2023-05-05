library(dplyr)
library(table1)     # table 1

setwd('UCSF/ipd-ma-cd2-2')

## Characteristics Table 1

tbl1 <- crohnsData_wk8_imputed %>% 
  # re-factor categorical variables
  mutate(SEX.FEMALE = abs(SEX.MALE - 1)) %>% 
  mutate(SEX.FEMALE = factor(SEX.FEMALE, 
                             levels=c(1,0), 
                             labels=c('Female','Male'))) %>% 
  mutate(HxOfTNFi   = factor(HxOfTNFi,
                             levels=c(1,0),
                             labels=c('Yes','No'))) %>% 
  mutate(STEROID    = factor(STEROID,
                             levels=c(1,0),
                             labels=c('Yes','No'))) %>% 
  mutate(IMMUNOMOD  = factor(IMMUNOMOD,
                             levels=c(1,0),
                             labels=c('Yes','No'))) %>% 
  mutate(LOC.ILEAL  = factor(LOC.ILEAL,
                             levels=c(1,0),
                             labels=c('Yes','No'))) %>% 
  filter(DRUG != 'ADA')

my.render.cont <- function(x){
  with(stats.apply.rounding(stats.default(x), digits=2), 
       c("", "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}

my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y, sprintf("%d (%0.0f %%)", FREQ, PCT))))
}

my.render <- function(x,...) {
  y <- render.default(x,...)
  if(is.factor(x)) y[2] else y
}

table1(~ TRTGRP + AGE + SEX.FEMALE + BMI + CDAI_BASELINE + CRP..mgL + HxOfTNFi + 
         STEROID + IMMUNOMOD + LOC.ILEAL | DRUG*TRIAL, 
       data=tbl1,
       overall=FALSE,
       # render = my.render,
       render.continuous = my.render.cont, 
       render.categorical= my.render.cat,
       droplevels = T)