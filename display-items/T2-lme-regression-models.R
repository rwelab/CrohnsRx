library(dplyr)
library(lme4)
library(sjPlot)    # Publication ready tables

setwd('UCSF/ipd-ma-cd2-2')

rm_plac <- readRDS('models/rm_plac.rds')
rm_il12 <- readRDS('models/rm_il12.rds')
rm_intg <- readRDS('models/rm_intg.rds')
rm_tnfi <- readRDS('models/rm_tnfi.rds')

# T2 Linear Mixed Effect Regression Models
tab_model(rm_plac, rm_il12, rm_intg, rm_tnfi,
          show.ci=F,
          show.se=T,
          pred.labels = c('Intercept','Year (Centered)','Baseline CDAI (Centered)','Age (Centered)','BMI (Centered)',
                          'CRP (mg/L) (Centered)','Sex: Male','HxOfTNFi','Steroid Use','Immunomodulator Use','Ileal Disease'),
          dv.labels = c('Placebo','Anti-Il-12/23','Anti-Integrin','Anti-TNF'),
          file='display-items/T2-lme-regression-models.html')