# predict.lme.R
# 
# Due to subtleties, standard errors calculations for linear mixed effect (lme) 
# models are not available in the lme4 package. We created our own predict.lme 
# function (~predict.lm) that can calculate the confidence and prediction standard 
# errors for new data. One of two methods - analytically or bootstrapping - can be 
# specified when calculating standard errors. References are provided with each
# method function call below. 
#
# Note, these methods were catered towards our specific use case, so this function
# may not be appropriate for replication. Please use your own judgement.   

#------------------------------------------------------------------------------#

library(dplyr)
library(lme4)
library(merTools)
library(parallel)

#------------------------------------------------------------------------------#

predict.lme <- function(object, newdata, interval = c('none','confidence','prediction'), 
           level = 0.95, method = c('analytical','bootstrap'), 
           model.name = c('plac','il12','intg','tnfi'), nsim = 100, 
           parallel = c('no','multicore'), ncpus = 1, seed = 1234) 
{
    if(missing(model.name)) {
      warning("make sure to specify the appropriate model.name = c('plac','il12','intg','tnfi')")
    }
  
    interval = match.arg(interval)
    method = match.arg(method)
    model.name = match.arg(model.name)
    parallel = match.arg(parallel)
    
    if(interval == 'none') {
      newdata[[ paste0(model.name, '.attrib') ]] = predict(object, newdata, re.form = NA)
    }
    
    if(interval == 'confidence' & method == 'bootstrap') {
      newdata = bootstrap.CI(object, newdata, model.name, nsim, parallel, ncpus, seed)
    }
    
    if(interval == 'prediction' & method == 'bootstrap') {
      newdata = bootstrap.PI(object, newdata, model.name, level, nsim, seed)
    }
    
    if(method == 'analytical') {
      newdata = analytical(object, newdata, model.name, interval=interval)
    }
    
    return(newdata)
}

#------------------------------------------------------------------------------#

# REF: [add here]
analytical <- function(object, newdata, model.name, interval = c('confidence','prediction'))
{
  interval = match.arg(interval)
  
  if(model.name == 'plac') {
    x = newdata %>% dplyr::select(Year_Cent:Ileal) # include Year_Cent
  }
  else {
    x = newdata %>% dplyr::select(CDAI_baseline_Cent:Ileal) # remove Year_Cent
  }

  x = as.matrix(cbind(1, x))
  v = vcov(object)
  # SE of CI
  var = diag(x %*% v %*% t(x))
  
  if(interval == 'prediction') {
    var_model = as.data.frame(VarCorr(object))
    var_trial = var_model[1,4] # random effect var
    var_res   = var_model[2,4] # residual var
    var_total = var_trial + var_res
    # SE of PI
    var = var + var_total
  }
  
  # append prediction to data
  newdata[[ paste0(model.name, '.attrib') ]] = predict(object, newdata, re.form = NA)
  newdata[[ paste0(model.name, '.se') ]]     = sqrt(var)
  
  return(newdata)
}

#------------------------------------------------------------------------------#

# REF: [add here]
bootstrap.CI <- function(object, newdata, model.name, nsim=100, parallel='no', ncpus=1, seed=1234)
{
  
  # return predicted values from bootstrap (random effects set to 0)
  myPred <- function(.){ predict(., newdata = newdata, re.form = NA) }
  
  boot.pred <- lme4::bootMer(object, 
                             myPred, 
                             nsim     = nsim,         # number of simulations
                             use.u    = TRUE,         # pre-defined
                             type     = "parametric", # pre-defined
                             seed     = seed, 
                             parallel = parallel,     # no, multicore 
                             ncpus    = ncpus)        # find max # CPUs using parallel::detectCores()
  
  # find median (drug class attrib) and standard deviation (standard error) for bootstrap
  sumBoot <- function(merBoot){
    return(
      data.frame(attrib = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
                 se     = apply(merBoot$t, 2, sd)
      )
    )
  }
  
  # attrib, se
  boot.fit <- sumBoot( boot.pred )
  
  # append prediction to data
  newdata[[ paste0(model.name, '.attrib') ]] = boot.fit$attrib 
  newdata[[ paste0(model.name, '.se') ]]     = boot.fit$se
  
  return(newdata)
}

#------------------------------------------------------------------------------#

# REF: https://github.com/jknowles/merTools/issues/60
bootstrap.PI <- function(object, newdata, model.name, level=0.95, nsim=100, seed=1234) 
{
  
  # set RE to 0 (avg), ensures equal number of columns
  newdata <- newdata %>% mutate(Trial = merTools::averageObs(object)$Trial)
  
  PI <- predictInterval(merMod = object, newdata = newdata,
                        level = level, n.sims = nsim, which = 'fixed',
                        stat = "median", type="linear.prediction",
                        include.resid.var = TRUE, seed = seed)
  
  newdata[[ paste0(model.name, '.attrib') ]] <- PI$fit
  newdata[[ paste0(model.name, '.se') ]]  <- (PI$upr - PI$lwr) / 3.92
  
  newdata <- newdata %>% dplyr::select(-Trial)
  
  return(newdata)
}

#------------------------------------------------------------------------------#

predict.lme.all <- function(object.list, newdata, model.name.list = c('plac','il12','intg','tnfi'), 
                            interval = c('none','confidence','prediction'), 
                            level = 0.95, method = c('analytical','bootstrap'), 
                            nsim = 100, parallel = c('no','multicore'), ncpus = 1, seed = 1234) 
{
  interval = match.arg(interval)
  method = match.arg(method)
  parallel = match.arg(parallel)
  
  for(i in 1:length(object.list)) {
    newdata <- predict.lme(object.list[[i]], newdata=newdata, interval=interval,
                           level=level, method=method, model.name=model.name.list[[i]], 
                           nsim=nsim, parallel=parallel, ncpus=ncpus, seed=seed)
  }
  return(newdata)
}

#------------------------------------------------------------------------------#

