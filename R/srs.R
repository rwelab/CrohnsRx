#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Sequential Regression & Simulation (SRS)
#   1. Train placebo model
#   2. Subtract placebo attributable effect from data (isolate cdai reduction due to drug only)
#   3. Train drug models
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(dplyr)
library(lme4)
library(rlang)
# library(lmerTest)

#------------------------------------------------------------------------------#
# more testing

srs <- function(formula, data=NULL, treatment.field=NULL, placebo.label="Placebo", exclude.cov=list()) {
  
  if (missing(formula))         { stop("Argument 'formula' was not provided.\n") } 
  if (missing(data))            { stop("Argument 'data' was not provided.\n") } 
  if (missing(treatment.field)) { stop("Argument 'treatment.field' was not provided.\n") } 
  
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
  srs.plac <- lme4::lmer(formula, data=placebo_df)
  srs.obj$PlaceboModel <- srs.plac

  # Infer Drug Attributable Effects
  active_df$ObsTotalEffect     <- active_df[[ PLACEBO.RESP ]]
  active_df$PredPlaceboEffect  <- predict(srs.plac, newdata=active_df, re.form=~0)
  active_df$InferredDrugEffect <- active_df[[ PLACEBO.RESP ]] - active_df$PredPlaceboEffect

  # Train Drug Attributable Models
  formula2 <- updateFormula(formula, "InferredDrugEffect", exclude.cov)
  for (trt in TRT.LABELS) {
    trt_df <- active_df %>% dplyr::filter( {{ treatment.field }} == trt )
    m <- lme4::lmer(formula2, data=trt_df)
    srs.obj$DrugModels[[trt]] <- m
  }

  return(srs.obj)
}

updateFormula <- function(orig.formula, new.response=NULL, exclude.cov=list()) {
  
  formula2 <- orig.formula
  
  # Update new response variable Y1 ~. -> Y2 ~.
  formula2 <- update(formula2, as.formula(paste(new.response, "~ .")))
  
  # Remove dependent variables
  remove.cov.str <- paste("- ", exclude.cov, collapse = " ")
  formula2 <- update(formula2, paste(". ~ .", remove.cov.str))
  
  return(formula2)
}

rhs.vars <- function(f) { as.character(attr(terms(f), "term.labels")) }
lhs.vars <- function(f) { as.character(attr(terms(formula), "variables")[[2]])  }

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
