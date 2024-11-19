\#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# table-sup-6 Model Benchmarking.R
#
# Based on Supplementary Table 2: Model Evaluation from `A method for estimating 
# causal effects from heterogeneous clinical trials without a common control group 
# using sequential regression and simulation: an individual participant data 
# meta-analysis and validation study`, a Random Forest model resulted in a 
# superior mean 5-fold cross-validative RMSE compared to a linear mixed-effect
# model when conducting the sequential regression and simulation (SRS) method 
# on biologic treatments for Crohn's Disease.
# 
# As a model benchmark, we compare the subgroup similarity between a mixed-effect
# and a random forest based SRS analysis. 
#
# Below, we have modified the srs() function to be modular, in that it can accept
# custom functions for model training, prediction and 95% CI. 
#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(dplyr)
library(lme4)
library(ranger)
source("R/srs.R")
source("R/subgroups.R")

crohns <- read.csv("data/crohns_biologics_clean.csv")
load("data/formula.RData")

#------------------------------------------------------------------------------#

# Modular SRS - input train, predict functions specific to model
srs.mod <- function(formula, data, train.func, predict.func, 
                    treatment.field=Treatment, placebo.label="Placebo", exclude.cov=c('Year_Cent')) {
  
  if (missing(formula))         { stop("Argument 'formula' was not provided.\n") } 
  if (missing(data))            { stop("Argument 'data' was not provided.\n") } 
  
  srs.obj <- list()
  
  # Constants
  ## Y = CDAI_reduction ~ X + (1|RV)
  PLACEBO.RESP  <- as.character(attr(terms(formula), "variables")[[2]]) 
  ## Treatment = c('Placebo','TNFi','Il12','Integrin') - c('Placebo') = c('TNFi','Il12','Integrin')
  TRT.LABELS    <- setdiff(data %>% distinct( {{ treatment.field }} ) %>% pull( {{ treatment.field }} ) , placebo.label)
  
  # Data
  placebo_df <- data %>% dplyr::filter( {{ treatment.field }} == placebo.label )
  active_df  <- data %>% dplyr::filter( {{ treatment.field }} != placebo.label )
  
  # Train Placebo Attributable Model
  srs.plac <- train.func(formula, placebo_df)
  srs.obj$PlaceboModel <- srs.plac
  
  # Infer Drug Attributable Effects
  active_df$ObsTotalEffect     <- active_df[[ PLACEBO.RESP ]]
  active_df$PredPlaceboEffect  <- predict.func(srs.plac, active_df, formula)
  active_df$InferredDrugEffect <- active_df[[ PLACEBO.RESP ]] - active_df$PredPlaceboEffect
  
  # Train Drug Attributable Models
  formula2 <- updateFormula(formula, "InferredDrugEffect", exclude.cov)
  for (trt in TRT.LABELS) {
    trt_df <- active_df %>% dplyr::filter( {{ treatment.field }} == trt )
    m <- train.func(formula2, trt_df)
    srs.obj$DrugModels[[trt]] <- m
  }
  
  return(srs.obj)
}

# Modular CI - input CI function specific to model 
srs.ci <- function(srs.obj, data, ci.func) {
  # Calculate 95% CI
  name.list <- names(srs.obj)
  for(i in 1:length(srs.obj)) {
    data <- ci.func(srs.obj[[i]], data, name.list[[i]])
    print(paste0("Completed batch predict ",i,": ",name.list[[i]]))
  }
  return(data)
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Define SRS Modular Functions: Train, Predict, 95%CI
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Mixed Effects (lmer) Functions

LMER.train <- function(f, d) { 
  lme4::lmer(formula = f, data = d)  
}

LMER.predict <- function(m, d, ...) {
  predict(m, d, re.form = NA)
}

LMER.ci <- function(m, d, model.name) {
  # return predicted values from bootstrap (random effects set to 0)
  myPred <- function(.){ predict(., newdata = d, re.form = NA) }
  boot.pred <- lme4::bootMer(m, myPred, nsim = 10000, use.u = TRUE, type = "parametric", 
                             seed = 1234, parallel = 'multicore', ncpus = 8)
  # find median (drug class attrib) and standard deviation (standard error) for bootstrap
  sumBoot <- function(merBoot){
    return(data.frame(attrib = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))), 
                      se = apply(merBoot$t, 2, sd)))
  }
  boot.fit <- sumBoot( boot.pred )
  d[[ paste0(model.name, '.attrib') ]] = boot.fit$attrib 
  d[[ paste0(model.name, '.se') ]]     = boot.fit$se
  return(d)
}

#------------------------------------------------------------------------------#
# Random Forest Functions

RF.train <- function(f, d) {
  ranger(formula = f, data = d, num.trees=10^3, keep.inbag=TRUE)
}

RF.predict <- function(m, d, ...) {
  pred <- predict(m, data = d)
  pred$predictions
}

RF.ci <- function(m, d, model.name) {
  pred <- predict(m, data = d, type = "se", se.method = "infjack")
  d[[ paste0(model.name, '.attrib') ]] = pred$predictions
  d[[ paste0(model.name, '.se') ]] = pred$se
  return(d)
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Model Benchmarking
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Random Effects (lmer)
LMER.obj <- srs.mod( formula, crohns, LMER.train, LMER.predict )
LMER.cf <- srs.ci( LMER.obj$DrugModels, crohns, LMER.ci )
LMER.sbgrps <- srs.subgroups( LMER.cf, LMER.obj$DrugModels )
(res.lmer <- isEqual.pct(LMER.sbgrps$Subgroup, LMER.sbgrps$Subgroup) )
# 1

#------------------------------------------------------------------------------#

formula2 <- updateFormula(formula, exclude.cov=c('(1 | Trial)')) # remove rand eff

# Random Forest
RF.obj <- srs.mod( formula2, crohns, RF.train, RF.predict )
RF.cf <- srs.ci( RF.obj$DrugModels, crohns, RF.ci )
RF.sbgrps <- srs.subgroups( RF.cf, RF.obj$DrugModels, df=srs.df(LMER.obj$DrugModels))
(res.rf <- isEqual.pct(RF.sbgrps$Subgroup, LMER.sbgrps$Subgroup) )
# 0.9238997

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
