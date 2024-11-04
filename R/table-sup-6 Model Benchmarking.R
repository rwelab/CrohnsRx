crohns <- read.csv("data/crohns_biologics_clean.csv")
# crohns %>% glimpse()
load("data/formula.RData")
formula2 <- updateFormula(formula, exclude.cov=c('(1 | Trial)'))

apply_predict_with_ci <- function(data, model, algorithm, n_boot = 1000, ci_level = 0.95) {
  library(boot)  # Ensure the boot library is loaded for custom bootstrapping if needed
  
  # Internal function to calculate predictions and CI using bootstrapping
  get_bootstrap_ci <- function(pred_func, data, n_boot, ci_level) {
    boot_results <- boot(data, pred_func, R = n_boot)
    ci <- boot.ci(boot_results, type = "perc", conf = ci_level)
    return(list(predictions = pred_func(data), ci = ci))
  }
  
  # Initialize output
  output <- NULL
  
  # Switch to select algorithm-specific prediction and bootstrapping methods
  output <- switch(algorithm,
                   "lm" = {
                     pred_func <- function(data) predict(model, newdata = data)
                     get_bootstrap_ci(pred_func, data, n_boot, ci_level)
                   },
                   "lme4" = {
                     # For mixed models, use bootMer for bootstrapping
                     pred_func <- function(fit) predict(fit, newdata = data)
                     boot_results <- lme4::bootMer(model, pred_func, nsim = n_boot)
                     ci <- quantile(boot_results$t, probs = c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2))
                     list(predictions = predict(model, newdata = data), ci = ci)
                   },
                   "glmnet" = {
                     pred_func <- function(data) predict(model, newx = as.matrix(data), s = "lambda.min")
                     get_bootstrap_ci(pred_func, data, n_boot, ci_level)
                   },
                   "rf" = {
                     pred_func <- function(data) predict(model, newdata = data)
                     get_bootstrap_ci(pred_func, data, n_boot, ci_level)
                   },
                   stop(paste("Algorithm", algorithm, "not recognized"))
  )
  
  return(output)
}



# Linear Regression (lm)

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
  srs.plac <- stats::lm(formula, data=placebo_df)
  srs.obj$PlaceboModel <- srs.plac
  
  # Infer Drug Attributable Effects
  active_df$ObsTotalEffect     <- active_df[[ PLACEBO.RESP ]]
  active_df$PredPlaceboEffect  <- stats::predict(srs.plac, newdata=active_df)
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