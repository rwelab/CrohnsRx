#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# table-sup-2 Anti-IL12 Cohort Baseline Char.R
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(dplyr)
library(table1)     # table 1

#------------------------------------------------------------------------------#

# Load data
crohns <- read.csv("data/crohns_biologics_clean.csv")
load("data/formula.RData")

# Train placebo, drug attributable models
srs.obj <- srs(formula         = formula, 
               data            = crohns, 
               treatment.field = Treatment, 
               placebo.label   = "Placebo",
               exclude.cov     = c("Year_Cent"))

# Find Subgroups
srs.cf <- batch.srs.predict(object.list = srs.obj$DrugModels, 
                            newdata     = crohns, 
                            level       = 0.95, 
                            interval    = 'confidence', 
                            method      = 'bootstrap', 
                            nsim        = 10000, 
                            parallel    = 'multicore', 
                            ncpus       = 8, 
                            seed        = 1234)

srs.sbgrps <- srs.summarise( srs.subgroups(srs.cf, srs.obj$DrugModels) )
knitr::kable(srs.sbgrps)

#------------------------------------------------------------------------------#

# IL12 > (TNFI, INTG)
prefer_il12 <- srs.sbgrps %>% 
  filter(p1.ohe == '>' & drug1 == 'il12')

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
                             labels=c('Yes','No'))) %>% 
  mutate(Age        = Age_Cent + 35) %>%
  mutate(BMI        = BMI_Cent + 20) %>%
  mutate(CRP        = CRP_Cent + 10) %>% 
  mutate(CDAI_baseline = CDAI_baseline_Cent + 300)

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

#------------------------------------------------------------------------------#

table1(~ Age + BMI + CDAI_baseline + CRP + HxOfTNFi + 
         SteroidUse + ImmUse + Ileal | Sex, 
       data=tbl1,
       render = my.render, 
       render.continuous = my.render.cont, 
       render.categorical= my.render.cat,
       droplevels = T)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#