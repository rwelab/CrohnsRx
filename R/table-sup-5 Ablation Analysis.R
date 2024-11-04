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

# Function to calculate sensitivity analysis subgroup accuracy
sa.acc <- function(res) { mean(res$isEqual) }

#------------------------------------------------------------------------------#
# Ground Truth Subgroups 

crohns <- read.csv("data/crohns_biologics_clean.csv")
abl.truth <- ablation.analysis(crohns, formula)
(res.truth <- sa.acc( srs.isEqual(abl.truth$Subgroup, abl.truth$Subgroup) ))
# 1

#------------------------------------------------------------------------------#
# Ablation 1 : -CDAI_baseline

for.cdai <- updateFormula(formula, exclude.cov = c("CDAI_baseline_Cent"))
abl.cdai <- ablation.analysis(crohns, for.cdai)
(res.cdai <- sa.acc( srs.isEqual(abl.truth$Subgroup, abl.cdai$Subgroup) ))
# 1

#------------------------------------------------------------------------------#
# Ablation 2 : -CRP

for.crp <- updateFormula(formula, exclude.cov = c("CRP_Cent"))
abl.crp <- ablation.analysis(crohns, for.crp)
(res.crp <- sa.acc( srs.isEqual(abl.truth$Subgroup, abl.crp$Subgroup) ))
# 1

#------------------------------------------------------------------------------#
# Ablation 3 : -Ileal

for.ileal <- updateFormula(formula, exclude.cov = c("Ileal"))
abl.ileal <- ablation.analysis(crohns, for.ileal)
(res.ileal <- sa.acc( srs.isEqual(abl.truth$Subgroup, abl.ileal$Subgroup) ))
# 1

#------------------------------------------------------------------------------#
# Ablation 4 : -CDAI_baseline, CRP, Ileal

for.all <- updateFormula(formula, exclude.cov = c("CDAI_baseline","CRP_Cent","Ileal"))
abl.all <- ablation.analysis(crohns, for.all)
(res.all <- sa.acc( srs.isEqual(abl.truth$Subgroup, abl.all$Subgroup) ))
# 1

#------------------------------------------------------------------------------#
# Print results

knitr::kable( 
  data.frame(
    `Ablation Column` = c("Ground Truth","CDAI Baseline","CRP","Ileal","CDAI Baseline, CRP, Ileal"),
    Accuracy          = c(res.truth, res.cdai, res.crp, res.ileal, res.all) 
  )
)

# |Ablation.Column           |  Accuracy|
# |:-------------------------|---------:|
# |Ground Truth              | 1.0000000|
# |CDAI Baseline             | 0.7976504|
# |CRP                       | 0.8593723|
# |Ileal                     | 0.9535332|
# |CDAI Baseline, CRP, Ileal | 0.8369279|

# |Ablation.Column           | Accuracy|
# |:-------------------------|--------:|
# |Ground Truth              |        1|
# |CDAI Baseline             |        1|
# |CRP                       |        1|
# |Ileal                     |        1|
# |CDAI Baseline, CRP, Ileal |        1|

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
