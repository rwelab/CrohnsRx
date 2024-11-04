#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# table-1 RCT Baseline Char.R
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(dplyr)
library(table1)     # table 1

#------------------------------------------------------------------------------#

# Load data
crohns <- read.csv("data/crohns_biologics_raw.csv") %>%
  dplyr::select(Group, Trial, DrugClass, Year, CDAI_baseline, Age, BMI, 
                CRP, HxOfTNFi, Sex_Male, SteroidUse, ImmUse, Ileal)

#------------------------------------------------------------------------------#

tbl1 <- crohns %>% 
  # re-factor categorical variables
  mutate(Sex_Female = abs(Sex_Male - 1)) %>% 
  mutate(Sex_Female = factor(Sex_Female, 
                             levels=c(1,0), 
                             labels=c('Female','Male'))) %>% 
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

table1(~ Group + Age + Sex_Female + BMI + CDAI_baseline + CRP + HxOfTNFi + 
         SteroidUse + ImmUse + Ileal | DrugClass*Trial, 
       data=tbl1,
       overall=FALSE,
       render.continuous = my.render.cont, 
       render.categorical= my.render.cat,
       droplevels = T)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
