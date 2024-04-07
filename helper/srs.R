################################################################################
# Sequential Regression & Simulation (SRS)
#   1. Train placebo model
#   2. Subtract placebo attributable effect from data (isolate cdai reduction due to drug only)
#   3. Train drug models
################################################################################

library(dplyr)
library(tidyr)
library(lme4)
source('helper/predict.lme.R') # predict.lme()
source('helper/subgroup.R')    # patient.drug.preferences(), subgroup.cohorts()

#------------------------------------------------------------------------------#

srs <- function(data) {
  
  # a. placebo model - train placebo attributable effect model
  placebo_df <- data %>% 
    filter(Group == 'Placebo') %>% 
    dplyr::select(Trial, Year_Cent:Ileal, CDAI_reduction)
  
  dep_var = 'CDAI_reduction'
  covariate_list <- names(placebo_df)[-which(names(placebo_df) %in% c('Trial',dep_var))]
  f1 <- paste(dep_var, paste(covariate_list, collapse = '+'), sep='~') # lm formula
  f2 <- paste0( f1 , '+(1|Trial)' ) # lmer formula
  
  int_plac <- lme4::lmer('CDAI_reduction ~ (1|Trial)', data=placebo_df)
  fm_plac <- lm(f1, data = placebo_df)
  rm_plac <- lme4::lmer(f2, data = placebo_df)

  #----------------------------------------------------------------------------#
  # b. find 'drug reduction' = 'CDAI_reduction' - 'Placebo_attributable' 
  fm_crohns_data1 <- data %>%
    mutate(Placebo_attributable = predict(fm_plac, newdata=data), 
           .after='CDAI_reduction') %>%
    mutate(Drug_reduction = CDAI_reduction - Placebo_attributable, 
           .after='Placebo_attributable')
  
  rm_crohns_data1 <- data %>%
    mutate(Placebo_attributable = predict(rm_plac, newdata=data, re.form=~0), 
           .after='CDAI_reduction') %>%
    mutate(Drug_reduction = CDAI_reduction - Placebo_attributable, 
           .after='Placebo_attributable')
  
  #----------------------------------------------------------------------------#
  ## c. drug class models
  fm_active_df <- fm_crohns_data1 %>% filter(Group == 'Active')
  fm_h2h_tnfi <- fm_active_df %>% filter(TNFi_Active == 1)
  fm_h2h_il12 <- fm_active_df %>% filter(Il12_Active == 1) 
  fm_h2h_intg <- fm_active_df %>% filter(Integrin_Active == 1)
  
  rm_active_df <- rm_crohns_data1 %>% filter(Group == 'Active')
  rm_h2h_tnfi <- rm_active_df %>% filter(TNFi_Active == 1)
  rm_h2h_il12 <- rm_active_df %>% filter(Il12_Active == 1) 
  rm_h2h_intg <- rm_active_df %>% filter(Integrin_Active == 1)
  
  covariate_list_2 <- covariate_list[-1] # remove Year_Cent
  dep_var = 'Drug_reduction'
  f3 <- paste(dep_var, paste(covariate_list_2, collapse = '+'), sep='~')
  f4 <- paste0( paste(dep_var, paste(covariate_list_2, collapse = '+'), sep='~'), '+(1|Trial)' )
  
  int_tnfi <- lme4::lmer('Drug_reduction ~ (1|Trial)', data=rm_h2h_tnfi)
  int_il12 <- lme4::lmer('Drug_reduction ~ (1|Trial)', data=rm_h2h_il12)
  int_intg <- lme4::lmer('Drug_reduction ~ (1|Trial)', data=rm_h2h_intg)
  fm_tnfi <- lm(f3, data=fm_h2h_tnfi)
  fm_il12 <- lm(f3, data=fm_h2h_il12)
  fm_intg <- lm(f3, data=fm_h2h_intg)
  rm_tnfi <- lme4::lmer(f4, rm_h2h_tnfi)
  rm_il12 <- lme4::lmer(f4, rm_h2h_il12)
  rm_intg <- lme4::lmer(f4, rm_h2h_intg)
  
  #----------------------------------------------------------------------------#
  # results <- list("df_placebo"=placebo_df,
  #                 "df_active_fm"=fm_crohns_data1,
  #                 "df_active_rm"=rm_crohns_data1,
  #                 "plac_00"=int_plac,
  #                 "plac_fm"=fm_plac,
  #                 "plac_rm"=rm_plac,
  #                 "il12_00"=int_il12,
  #                 "il12_fm"=fm_il12,
  #                 "il12_rm"=rm_il12,
  #                 "intg_00"=int_intg,
  #                 "intg_fm"=fm_intg,
  #                 "intg_rm"=rm_intg,
  #                 "tnfi_00"=int_tnfi,
  #                 "tnfi_fm"=fm_tnfi,
  #                 "tnfi_rm"=rm_tnfi,
  #                 "f1"=f1,
  #                 "f2"=f2,
  #                 "f3"=f3,
  #                 "f4"=f4)
  
  results <- list("plac_rm"=rm_plac,
                  "il12_rm"=rm_il12,
                  "intg_rm"=rm_intg,
                  "tnfi_rm"=rm_tnfi)
  return(results)
}

#------------------------------------------------------------------------------#

srs.placebo <- function(data) {
  placebo_df <- data %>% 
    filter(Group == 'Placebo') %>% 
    dplyr::select(Trial, Year_Cent:Ileal, CDAI_reduction)
  
  dep_var = 'CDAI_reduction'
  covariate_list <- names(placebo_df)[-which(names(placebo_df) %in% c('Trial',dep_var))]
  f1 <- paste(dep_var, paste(covariate_list, collapse = '+'), sep='~') # lm formula
  f2 <- paste0( f1 , '+(1|Trial)' ) # lmer formula
  
  rm_plac <- lme4::lmer(f2, data = placebo_df)
  
  results <- list("plac_rm"=rm_plac)
  return(results)
}

srs.active <- function(data, rm_plac) {
  rm_crohns_data1 <- data %>%
    mutate(Placebo_attributable = predict(rm_plac, newdata=data, re.form=~0), 
           .after='CDAI_reduction') %>%
    mutate(Drug_reduction = CDAI_reduction - Placebo_attributable, 
           .after='Placebo_attributable')
  
  rm_active_df <- rm_crohns_data1 %>% filter(Group == 'Active')
  rm_h2h_tnfi <- rm_active_df %>% filter(TNFi_Active == 1)
  rm_h2h_il12 <- rm_active_df %>% filter(Il12_Active == 1) 
  rm_h2h_intg <- rm_active_df %>% filter(Integrin_Active == 1)
  
  covariate_list <- names(data %>% dplyr::select(Year_Cent:Ileal))
  covariate_list_2 <- covariate_list[-1] # remove Year_Cent
  dep_var = 'Drug_reduction'
  f3 <- paste(dep_var, paste(covariate_list_2, collapse = '+'), sep='~')
  f4 <- paste0( paste(dep_var, paste(covariate_list_2, collapse = '+'), sep='~'), '+(1|Trial)' )
  
  rm_tnfi <- lme4::lmer(f4, rm_h2h_tnfi)
  rm_il12 <- lme4::lmer(f4, rm_h2h_il12)
  rm_intg <- lme4::lmer(f4, rm_h2h_intg)
  
  results <- list("il12_rm"=rm_il12,
                  "intg_rm"=rm_intg,
                  "tnfi_rm"=rm_tnfi)
  return(results)
}

#------------------------------------------------------------------------------#

subgroup.analysis <- function(data, rm_plac, rm_il12, rm_intg, rm_tnfi) {
  
  data <- predict.lme(rm_plac, data, model.name='plac', 
                             interval='confidence', level=0.95, method='bootstrap', 
                             nsim=10000, parallel='multicore', ncpus=8, seed=1234)
  data <- predict.lme(rm_il12, data, model.name='il12', 
                             interval='confidence', level=0.95, method='bootstrap', 
                             nsim=10000, parallel='multicore', ncpus=8, seed=1234)
  data <- predict.lme(rm_intg, data, model.name='intg', 
                             interval='confidence', level=0.95, method='bootstrap',  
                             nsim=10000, parallel='multicore', ncpus=8, seed=1234)
  data <- predict.lme(rm_tnfi, data, model.name='tnfi', 
                             interval='confidence', level=0.95, method='bootstrap',  
                             nsim=10000, parallel='multicore', ncpus=8, seed=1234)
  
  df_list <- list("il12" = nrow(rm_il12@frame) - 10,
                  "intg" = nrow(rm_intg@frame) - 10,
                  "tnfi" = nrow(rm_tnfi@frame) - 10)
  
  treatment_preferences <- patient.drug.preferences(data, df_list)
  
  return(treatment_preferences)
} 

#------------------------------------------------------------------------------#
