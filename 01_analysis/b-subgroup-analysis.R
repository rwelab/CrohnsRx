library(dplyr)
library(tidyr)
library(lme4)     # bootMer() for bootstrapping  
library(parallel) # parallel processing

#------------------------------------------------------------------------------#

setwd('UCSF/ipd-ma-cd2-2')

rm_plac <- readRDS('models/rm_plac.rds')
rm_il12 <- readRDS('models/rm_il12.rds')
rm_intg <- readRDS('models/rm_intg.rds')
rm_tnfi <- readRDS('models/rm_tnfi.rds')

# see helper files for more details
source('helper/predict.lme.R') # predict.lme()
source('helper/subgroup.R')    # patient.drug.preferences(), subgroup.cohorts()

#------------------------------------------------------------------------------#

## data frame
cd_patients <- read.csv("data/crohns_data_raw.csv") %>% 
  # rename
  rename(
    Year_Cent = Year_Norm,                   # Year - 2000
    CDAI_baseline_Cent = CDAI_baseline_Norm, # CDAI_baseline - 300
    Age_Cent = Age_Norm,                     # Age - 35
    BMI_Cent = BMI_Norm,                     # BMI - 20
    CRP_Cent = CRP_Norm                      # CRP (mg/L) - 10
  ) %>% 
  # remove columns
  dplyr::select(Year, CDAI_baseline, Age, BMI, CRP, 
                Year_Cent, CDAI_baseline_Cent, Age_Cent, BMI_Cent, CRP_Cent, 
                HxOfTNFi, Sex_Male, SteroidUse, ImmUse, Ileal)

cd_patients %>% glimpse()

#------------------------------------------------------------------------------#

# determine number of cores available
# only linux and macOS can perform multicore processing in R 
# specify parallel to 'multicore' and ncpus > 1
# if on windows, specify parallel to 'no' and ncpus = 1
parallel::detectCores(logical = FALSE)

# predict drug class attributable (.attrib) and standard error of mean response (.se)
# for each trial-based patient for each model. Performing a 10,000 simulation 
# parametric bootstrap (use.u = TRUE). 
# 8 cores, 10,000 simulations ~ 10-20 seconds per model prediction
cd_patients <- predict.lme(rm_plac, cd_patients, model.name='plac', 
                           interval='confidence', level=0.95, method='bootstrap', 
                           nsim=10000, parallel='multicore', ncpus=8, seed=1234)

cd_patients <- predict.lme(rm_il12, cd_patients, model.name='il12', 
                           interval='confidence', level=0.95, method='bootstrap', 
                           nsim=10000, parallel='multicore', ncpus=8, seed=1234)

cd_patients <- predict.lme(rm_intg, cd_patients, model.name='intg', 
                           interval='confidence', level=0.95, method='bootstrap',  
                           nsim=10000, parallel='multicore', ncpus=8, seed=1234)

cd_patients <- predict.lme(rm_tnfi, cd_patients, model.name='tnfi', 
                           interval='confidence', level=0.95, method='bootstrap',  
                           nsim=10000, parallel='multicore', ncpus=8, seed=1234)

cd_patients %>% glimpse()

#------------------------------------------------------------------------------#

# determine drug class ranking and preferences

## model degrees of freedom (for two-sample t-test evaluation : check if two 
## drug classes are significantly different from one another for a patient)
df_list <- list("il12" = nrow(rm_il12@frame) - 10, #  587 - 10
                "intg" = nrow(rm_intg@frame) - 10, # 1818 - 10
                "tnfi" = nrow(rm_tnfi@frame) - 10) # 1677 - 10

treatment_preferences <- patient.drug.preferences(cd_patients, df_list)

treatment_preferences %>% glimpse()

#------------------------------------------------------------------------------#

# aggregate patients into cohorts based on similar treatment preferences

subgroup.cohorts(treatment_preferences)
#   drug1 drug2 drug3 p12_ohe p23_ohe     n
# 1 drug1 drug2 drug3       0       0  3142
# 2 drug1 drug2 il12        0       1     4
# 3 drug1 drug2 intg        0       1   354
# 4 il12  drug2 drug3       1       0   138
# 5 il12  tnfi  intg        1       1     1
# 6 tnfi  drug2 drug3       1       0  2021
# 7 tnfi  intg  il12        1       1    43

#------------------------------------------------------------------------------#

# save
write.csv(treatment_preferences, 'data/treatment_preferences-n10000.csv')

#------------------------------------------------------------------------------#

