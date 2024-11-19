#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# table-sup-5 Ablation Analysis.R
#
# Testing subgroup accuracy when certain input parameters are ablated (missing). 
# We chose to ablate CDAI baseline, CRP, and ileal involvement, as those three 
# variables in the RCT data are not often recorded or are difficult to calculate.
#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

setwd('UCSF-R/ipd-ma-cd2-3')

# Load functions
source("R/srs.R")
source("R/predict.R")
source("R/subgroups.R")
load("data/formula.RData")

# SRS wrapper 
ablation.analysis <- function(data, f) {
  srs.obj <- srs(formula         = f, 
                 data            = data, 
                 treatment.field = Treatment, 
                 placebo.label   = "Placebo",
                 exclude.cov     = c("Year_Cent"))
  
  srs.cf <- batch.srs.predict(object.list = srs.obj$DrugModels, 
                              newdata     = data, 
                              level       = 0.95,
                              interval    = 'confidence', 
                              method      = 'bootstrap', 
                              nsim        = 10000,
                              parallel    = "multicore",
                              ncpus       = 8,
                              seed        = 1234)
  srs.sbg <- srs.subgroups(srs.cf, srs.obj$DrugModels)
  return(srs.sbg)
}

#------------------------------------------------------------------------------#
# Ablation Analysis 1: % is Equal
#------------------------------------------------------------------------------#
# Ground Truth Subgroups 

crohns <- read.csv("data/crohns_biologics_clean.csv")
abl.truth <- ablation.analysis(crohns, formula)
(res.truth <- isEqual.pct(abl.truth$Subgroup, abl.truth$Subgroup) )
# 1

#------------------------------------------------------------------------------#
# Ablation 1 : -CDAI_baseline

for.cdai <- updateFormula(formula, exclude.cov = c("CDAI_baseline_Cent"))
abl.cdai <- ablation.analysis(crohns, for.cdai)
(res.cdai <- isEqual.pct(abl.truth$Subgroup, abl.cdai$Subgroup) )
# 1
(npr.cdai <- noPref.pct(abl.cdai))
# 0.63

#------------------------------------------------------------------------------#
# Ablation 2 : -CRP

for.crp <- updateFormula(formula, exclude.cov = c("CRP_Cent"))
abl.crp <- ablation.analysis(crohns, for.crp)
(res.crp <- isEqual.pct(abl.truth$Subgroup, abl.crp$Subgroup) )
# 1
(npr.crp <- noPref.pct(abl.crp))
# 0.53

#------------------------------------------------------------------------------#
# Ablation 3 : -Ileal

for.ileal <- updateFormula(formula, exclude.cov = c("Ileal"))
abl.ileal <- ablation.analysis(crohns, for.ileal)
(res.ileal <- isEqual.pct(abl.truth$Subgroup, abl.ileal$Subgroup) )
# 1
(npr.ileal <- noPref.pct(abl.ileal))
# 0.50

#------------------------------------------------------------------------------#
# Ablation 4 : -CDAI_baseline, CRP, Ileal

for.all <- updateFormula(formula, exclude.cov = c("CDAI_baseline","CRP_Cent","Ileal"))
abl.all <- ablation.analysis(crohns, for.all)
(res.all <- isEqual.pct(abl.truth$Subgroup, abl.all$Subgroup) )
# 1
(npr.all <- noPref.pct(abl.all))
# 0.48

#------------------------------------------------------------------------------#

knitr::kable( 
  data.frame(
    `Ablation Column` = c("Ground Truth","CDAI Baseline","CRP","Ileal","CDAI Baseline, CRP, Ileal"),
    Accuracy          = c(res.truth, res.cdai, res.crp, res.ileal, res.all) 
  )
)
# |Ablation.Column           | Accuracy|
# |:-------------------------|--------:|
# |Ground Truth              |        1|
# |CDAI Baseline             |        1|
# |CRP                       |        1|
# |Ileal                     |        1|
# |CDAI Baseline, CRP, Ileal |        1|

#------------------------------------------------------------------------------#
# Ablation Analysis 2: No Preference % 
#------------------------------------------------------------------------------#
# Ground Truth 
(npr.truth <- noPref.pct(abl.truth))
# 0.55

#------------------------------------------------------------------------------#
# Ablation 1 : -CDAI_baseline
(npr.cdai <- noPref.pct(abl.cdai))
# 0.63

#------------------------------------------------------------------------------#
# Ablation 2 : -CRP
(npr.crp <- noPref.pct(abl.crp))
# 0.53

#------------------------------------------------------------------------------#
# Ablation 3 : -Ileal
(npr.ileal <- noPref.pct(abl.ileal))
# 0.50

#------------------------------------------------------------------------------#
# Ablation 4 : -CDAI_baseline, CRP, Ileal
(npr.all <- noPref.pct(abl.all))
# 0.48

#------------------------------------------------------------------------------#

knitr::kable( 
  data.frame(
    `Ablation Column` = c("Ground Truth","CDAI Baseline","CRP","Ileal","CDAI Baseline, CRP, Ileal"),
    `No Preference`   = c(npr.truth, npr.cdai, npr.crp, npr.ileal, npr.all) 
  )
)
# |Ablation.Column           | No.Preference|
# |:-------------------------|-------------:|
# |Ground Truth              |          0.55|
# |CDAI Baseline             |          0.63|
# |CRP                       |          0.53|
# |Ileal                     |          0.50|
# |CDAI Baseline, CRP, Ileal |          0.48|

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
