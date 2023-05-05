library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)  # get p-values for lmer covariates

setwd('UCSF/ipd-ma-cd2-2')

# data frame
crohns_data1 <- read.csv("data/crohns_data_raw.csv") %>% 
  # rename
  rename(
    Year_Cent = Year_Norm,                   # Year - 2000
    CDAI_baseline_Cent = CDAI_baseline_Norm, # CDAI_baseline - 300
    Age_Cent = Age_Norm,                     # Age - 35
    BMI_Cent = BMI_Norm,                     # BMI - 20
    CRP_Cent = CRP_Norm                      # CRP (mg/L) - 10
  ) %>% 
  dplyr::select(
    # Participant group: active (drug) or placebo
    Group, 
    # Random effect variables
    Trial, 
    # Independent variables
    Year_Cent, CDAI_baseline_Cent, Age_Cent, BMI_Cent, CRP_Cent, 
    HxOfTNFi, Sex_Male, SteroidUse, ImmUse, Ileal, 
    # Dependent variables
    CDAI_reduction)

crohns_data1 %>% glimpse()

################################################################################
# Sequential Regression & Simulation (SRS)
#   1. Train placebo model
#   2. Subtract placebo attributable effect from data (isolate cdai reduction due to drug only)
#   3. Train drug models
################################################################################

# a. placebo model - train placebo attributable effect model

## isolate placebo data
placebo_df <- crohns_data1 %>% 
  filter(Group == 'Placebo') %>% 
  dplyr::select(Trial, Year_Cent:Ileal, CDAI_reduction)
placebo_df %>% glimpse()

## build formula objects
dep_var = 'CDAI_reduction'
covariate_list <- names(placebo_df)[-which(names(placebo_df) %in% c('Trial',dep_var))]
f1 <- paste(dep_var, paste(covariate_list, collapse = '+'), sep='~') # lm formula
f2 <- paste0( f1 , '+(1|Trial)' ) # lmer formula

## compute models

### intercept only model
int_plac <- lme4::lmer('CDAI_reduction ~ (1|Trial)', data=placebo_df)

### fixed effect model
fm_plac <- lm(f1, data = placebo_df)
summary(fm_plac)

### random effect model 
rm_plac <- lme4::lmer(f2, data = placebo_df)
summary( lmerTest::lmer(f2, data = placebo_df) )
ranef(rm_plac)

## likelihood ratio test - linear mixed effect models > intercept only (averaged effects)
anova(rm_plac, int_plac)

#------------------------------------------------------------------------------#

## b. find 'drug reduction' = 'CDAI_reduction' - 'Placebo_attributable' 

fm_crohns_data1 <- crohns_data1 %>%
  
  # impute (predicted) placebo reduction using placebo model
  mutate(Placebo_attributable = predict(fm_plac, newdata=crohns_data1), 
         .after='CDAI_reduction') %>%
  
  # calculate drug reduction (attributable effect)
  mutate(Drug_reduction = CDAI_reduction - Placebo_attributable, 
         .after='Placebo_attributable')

rm_crohns_data1 <- crohns_data1 %>%
  
  # impute (predicted) placebo reduction using placebo model
  mutate(Placebo_attributable = predict(rm_plac, newdata=crohns_data1, re.form=~0), 
         .after='CDAI_reduction') %>%
  
  # calculate drug reduction (attributable effect)
  mutate(Drug_reduction = CDAI_reduction - Placebo_attributable, 
         .after='Placebo_attributable')

fm_crohns_data1 %>% glimpse()
rm_crohns_data1 %>% glimpse()

#------------------------------------------------------------------------------#

## c. drug class models

### isolate active data
fm_active_df <- fm_crohns_data1 %>% 
  filter(Group == 'Active')

rm_active_df <- rm_crohns_data1 %>%
  filter(Group == 'Active')

### Split into 3 training datasets by drug class type (TNFi, IL12, Integrin)
fm_h2h_tnfi <- fm_active_df %>% filter(TNFi_Active == 1)
fm_h2h_il12 <- fm_active_df %>% filter(Il12_Active == 1) 
fm_h2h_intg <- fm_active_df %>% filter(Integrin_Active == 1)

rm_h2h_tnfi <- rm_active_df %>% filter(TNFi_Active == 1)
rm_h2h_il12 <- rm_active_df %>% filter(Il12_Active == 1) 
rm_h2h_intg <- rm_active_df %>% filter(Integrin_Active == 1)

### build formula objects
covariate_list_2 <- covariate_list[-1] # remove Year_Cent
dep_var = 'Drug_reduction'
f3 <- paste(dep_var, paste(covariate_list_2, collapse = '+'), sep='~')
f4 <- paste0( paste(dep_var, paste(covariate_list_2, collapse = '+'), sep='~'), '+(1|Trial)' )
f4

# intercept only models
int_tnfi <- lme4::lmer('Drug_reduction ~ (1|Trial)', data=rm_h2h_tnfi)
int_il12 <- lme4::lmer('Drug_reduction ~ (1|Trial)', data=rm_h2h_il12)
int_intg <- lme4::lmer('Drug_reduction ~ (1|Trial)', data=rm_h2h_intg)

# fixed effect models 
fm_tnfi <- lm(f3, data=fm_h2h_tnfi)
fm_il12 <- lm(f3, data=fm_h2h_il12)
fm_intg <- lm(f3, data=fm_h2h_intg)

summary(fm_tnfi)
  
# random effect models
rm_tnfi <- lme4::lmer(f4, rm_h2h_tnfi)
rm_il12 <- lme4::lmer(f4, rm_h2h_il12)
rm_intg <- lme4::lmer(f4, rm_h2h_intg)

summary( lmerTest::lmer(f4, data=rm_h2h_tnfi) )
ranef(rm_tnfi)
summary( lmerTest::lmer(f4, data=rm_h2h_il12) )
ranef(rm_il12)
summary( lmerTest::lmer(f4, data=rm_h2h_intg) )
ranef(rm_intg)

### likelihood ratio test
anova(rm_tnfi, int_tnfi)
anova(rm_il12, int_il12)
anova(rm_intg, int_intg)

#------------------------------------------------------------------------------#

## Save
saveRDS( fm_plac ,  file = 'models/fm_plac.rds')
saveRDS( fm_il12 ,  file = 'models/fm_il12.rds')
saveRDS( fm_intg ,  file = 'models/fm_intg.rds')
saveRDS( fm_tnfi ,  file = 'models/fm_tnfi.rds')

saveRDS( rm_plac ,  file = 'models/rm_plac.rds')
saveRDS( rm_il12 ,  file = 'models/rm_il12.rds')
saveRDS( rm_intg ,  file = 'models/rm_intg.rds')
saveRDS( rm_tnfi ,  file = 'models/rm_tnfi.rds')
