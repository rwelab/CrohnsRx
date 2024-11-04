#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# table-2 Model Summaries.R
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(dplyr)
library(lme4)
library(sjPlot)    # Publication ready tables
source("R/srs.R")

#------------------------------------------------------------------------------#

crohns <- read.csv("data/crohns_biologics_clean.csv")
load("data/formula.RData")
srs.obj <- srs(formula         = formula, 
               data            = crohns, 
               treatment.field = Treatment, 
               placebo.label   = "Placebo",
               exclude.cov     = c("Year_Cent"))

#------------------------------------------------------------------------------#

# T2 Linear Mixed Effect Regression Models
tab_model(srs.obj$PlaceboModel, srs.obj$DrugModels$Il12, srs.obj$DrugModels$Integrin, srs.obj$DrugModels$TNFi,
          show.ci=F,
          show.se=T,
          pred.labels = c('Intercept','Year (Centered)','Baseline CDAI (Centered)','Age (Centered)','BMI (Centered)',
                          'CRP (mg/L) (Centered)','Sex: Male','HxOfTNFi','Steroid Use','Immunomodulator Use','Ileal Disease'),
          dv.labels = c('Placebo','Anti-Il-12/23','Anti-Integrin','Anti-TNF'),
          file='images/table-2.html')

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#