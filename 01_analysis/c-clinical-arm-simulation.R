library(dplyr)
library(lme4)
library(merTools) # predictInterval() for lme4 objects
options(warn=-1)

setwd('UCSF/ipd-ma-cd2-2')

rm_plac <- readRDS('models/rm_plac.rds')
rm_il12 <- readRDS('models/rm_il12.rds')
rm_tnfi <- readRDS('models/rm_tnfi.rds')

## data frame
crohns_data1 <- read.csv("data/crohns_data.csv") %>% 
  # rename
  rename(
    CDAI_baseline_Cent = CDAI_baseline_Norm,
    Age_Cent = Age_Norm, 
    BMI_Cent = BMI_Norm,
    CRP_Cent = CRP_Norm
  ) %>% 
  mutate(Year_Cent = mean(Year_Norm)) %>% 
  dplyr::select(Year_Cent, CDAI_baseline_Cent:Ileal)

crohns_data1 %>% glimpse()

#------------------------------------------------------------------------------#

# SIMULATION FUNCTIONS

sample_cohort <- function(x, size, replace=TRUE) {
  return( x[sample(nrow(x), size, replace=replace), ] )
}

PI.helper <- function(object, newdata, model.name, nsim=100) {
  
  # set RE to 0 (avg), ensures equal number of columns
  # SOURCE: https://github.com/jknowles/merTools/issues/60
  newdata <- newdata %>% mutate(Trial = merTools::averageObs(object)$Trial)
  
  PI <- predictInterval(merMod = object, newdata = newdata,
                        level = 0.95, n.sims = nsim, which = 'fixed',
                        stat = "median", type="linear.prediction",
                        include.resid.var = TRUE)
  
  fit.name <- paste(model.name, 'pi.fit', sep='.')
  sd.name <- paste(model.name, 'pi.sd', sep='.')
  
  newdata[fit.name] <- PI$fit
  newdata[sd.name]  <- (PI$upr - PI$lwr) / 3.92
  
  newdata <- newdata %>% dplyr::select(-Trial)
  
  return(newdata)
}

run_test <- function(ARM1, ARM2, arm.size) {
  
  # Confidence interval = (x1–x2) +/- t*√((sp2/n1)+(sp2/n2))
  
  n = arm.size
  mu1 = mean(ARM1) 
  sd1 = sd(ARM1)
  mu2 = mean(ARM2)
  sd2 = sd(ARM2)
  
  sp = ((n-1)*sd1^2+(n-1)*sd2^2)/(n+n-2) # pooled sd
  margin <- qt(0.975, df=n+n-1)*sqrt(sp/n + sp/n)
  
  lwr = (mu1-mu2) - margin
  upr = (mu1-mu2) + margin
  
  print(c('m1'=mu1,'s1'=sd1,'m2'=mu2,'s2'=sd2,'lwr'=lwr,'upr'=upr))
  
  # confidence interval
  return( ifelse((sign(lwr)==sign(upr)) & (mu1 > mu2), 1, 0) )
}

sim_trial <- function(x, arm.size=100) {
  ARM1 <- sample_cohort(x, size=arm.size)       # random sample w replacement
  ARM1 <- PI.helper(rm_il12, ARM1, 'il12', 100) # drug class fit, PI
  ARM1 <- ARM1 %>% 
    mutate(
      plac.fit = predict(rm_plac, newdata=ARM1, re.form=~0), 
      Response = 1 - pnorm(100, # clinical response
                           mean = plac.fit + il12.pi.fit,
                           sd = il12.pi.sd))
  
  ARM2 <- sample_cohort(x, size=arm.size)
  ARM2 <- PI.helper(rm_tnfi, ARM2, 'tnfi', 100)
  ARM2 <- ARM2 %>% 
    mutate(
      plac.fit = predict(rm_plac, newdata=ARM2, re.form=~0), 
      Response = 1 - pnorm(100, # clinical response
                           mean = plac.fit + tnfi.pi.fit,
                           sd = tnfi.pi.sd))
  
  return(run_test(ARM1$Response, ARM2$Response, arm.size))
}

set.seed(1234)
replicate(5, sim_trial(x=df1, arm.size=100))

#------------------------------------------------------------------------------#

## DATASET1: AGE >= 50
df1 <- crohns_data1 %>% filter((Age_Cent+35) >= 50) 
dim(df1) # 1040 / 5703
quantile(df1$Age_Cent, probs=c(0.25,0.50,0.75))

set.seed(1234)
t1.n100 <- replicate(1000, sim_trial(x=df1, arm.size = 100))
t1.n250 <- replicate(1000, sim_trial(x=df1, arm.size = 250))
t1.n500 <- replicate(1000, sim_trial(x=df1, arm.size = 500))

mean(t1.n100) # 0.585
mean(t1.n250) # 0.873
mean(t1.n500) # 0.974

#------------------------------------------------------------------------------#

## DATASET2: AGE >=50 & SEX_MALE == 0
df3 <- crohns_data1 %>% filter((Age_Cent+35) >= 50 & Sex_Male == 0)
dim(df3) # 577 / 5703
quantile(df3$Age_Cent, probs=c(0.25,0.50,0.75))

set.seed(1234)
t3.n100 <- replicate(1000, sim_trial(x=df3, arm.size = 100))
t3.n250 <- replicate(1000, sim_trial(x=df3, arm.size = 250))
t3.n500 <- replicate(1000, sim_trial(x=df3, arm.size = 500))

mean(t3.n100) # 0.755
mean(t3.n250) # 0.969
mean(t3.n500) # 0.999

#------------------------------------------------------------------------------#

## DATASET3: AGE >=50 | SEX_MALE == 1
df6 <- crohns_data1 %>% filter((Age_Cent+35) >= 50 & Sex_Male == 1)
dim(df6) # 463 / 5703
quantile(df6$Age_Cent, probs=c(0.25,0.50,0.75))

set.seed(1234)
t6.n100 <- replicate(1000, sim_trial(x=df6, arm.size = 100))
t6.n250 <- replicate(1000, sim_trial(x=df6, arm.size = 250))
t6.n500 <- replicate(1000, sim_trial(x=df6, arm.size = 500))

mean(t6.n100) # 0.346
mean(t6.n250) # 0.634
mean(t6.n500) # 0.807

#------------------------------------------------------------------------------#