library(dplyr) 
library(lme4)           # bootstrapping
library(parallel)       # parallel computing
# install.packages('microbenchmark')
library(microbenchmark) # benchmarking

#------------------------------------------------------------------------------#

setwd('UCSF/ipd-ma-cd2-2')

source('helper/predict.lme.R')  # bootstrap.CI.all()

rm_plac <- readRDS('models/rm_plac.rds')
rm_il12 <- readRDS('models/rm_il12.rds')
rm_tnfi <- readRDS('models/rm_tnfi.rds')
rm_intg <- readRDS('models/rm_intg.rds')

## data frame
crohns_data1 <- read.csv('data/crohns_data.csv') %>% 
  # rename
  rename(
    Year_Cent = Year_Norm, 
    CDAI_baseline_Cent = CDAI_baseline_Norm,
    Age_Cent = Age_Norm, 
    BMI_Cent = BMI_Norm,
    CRP_Cent = CRP_Norm
  ) %>% 
  # keep covariate columns
  dplyr::select(Year_Cent:Ileal)

#------------------------------------------------------------------------------#

# isolate one patient
test_pat <- crohns_data1[1,]

# predict on test_pat, re.form = NA (no random effects)
mySumm1 <- function(.) { predict(., newdata=test_pat, re.form=NA) }

# collapses bootstrap into median, 95% CI
sumBoot <- function(merBoot) {
  return(
    data.frame(attrib = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
               se     = apply(merBoot$t, 2, sd)
    )
  )
}

#------------------------------------------------------------------------------#
# test variable cores (1, 5, 10)

mbm.cores = microbenchmark(
  cores_01 = lme4::bootMer(rm_il12, mySumm1, use.u=TRUE, type="parametric", seed = 1234,
                           nsim=1000, 
                           parallel = 'no',
                           ncpus = 1), 
  cores_05 = lme4::bootMer(rm_il12, mySumm1, use.u=TRUE, type="parametric", seed = 1234,
                           nsim=1000, 
                           parallel = 'multicore', 
                           ncpus = 5), 
  cores_10 = lme4::bootMer(rm_il12, mySumm1, use.u=TRUE, type="parametric", seed = 1234,
                           nsim=1000, 
                           parallel = 'multicore', 
                           ncpus = 10),
  times = 10
)

mbm.cores
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# cores_01 6.108185 6.270554 6.349026 6.370023 6.418871 6.626994    50
# cores_05 1.479239 1.515393 1.548999 1.533453 1.553435 1.852336    50
# cores_10 1.024199 1.070884 1.105019 1.095506 1.124753 1.326431    50

#------------------------------------------------------------------------------#
# test fixed cores, different nsim (100, 1000, 10_000)

mbm.nsim = microbenchmark(
  nsim_100__ = lme4::bootMer(rm_il12, mySumm1, use.u=TRUE, type="parametric", seed = 1234,
                             nsim=100,
                             parallel = 'multicore', 
                             ncpus = 8), 
  nsim_1000_ = lme4::bootMer(rm_il12, mySumm1, use.u=TRUE, type="parametric", seed = 1234,
                             nsim=1000,
                             parallel = 'multicore', 
                             ncpus = 8),
  nsim_10000 = lme4::bootMer(rm_il12, mySumm1, use.u=TRUE, type="parametric", seed = 1234,
                             nsim=1000,
                             parallel = 'multicore', 
                             ncpus = 8),
  times = 10
)

mbm.nsim
# Unit: milliseconds
# expr      min         lq       mean     median        uq        max neval
# nsim_100__  146.822   172.8438   203.5946   208.2287   218.867   270.0401    10
# nsim_1000_ 1079.922  1091.5089  1163.1766  1119.2328  1179.994  1491.3509    10
# nsim_10000 9859.251 10200.3076 10305.6172 10239.8251 10477.390 10983.0560    10

#------------------------------------------------------------------------------#
# Fixed core (1 CPU), different nsim (100, 250, 500, 1000)

mbm.nsim = microbenchmark(
  nsim_100_ = bootstrap.CI.all(newdata = test_pat,
                               object_plac=rm_plac,
                               object_il12=rm_il12,
                               object_intg=rm_intg, 
                               object_tnfi=rm_tnfi,
                               nsim = 100, 
                               parallel = 'no', 
                               ncpus = 1), 
  nsim_250_ = bootstrap.CI.all(newdata = test_pat,
                               object_plac=rm_plac,
                               object_il12=rm_il12,
                               object_intg=rm_intg, 
                               object_tnfi=rm_tnfi,
                               nsim = 250, 
                               parallel = 'no', 
                               ncpus = 1), 
  nsim_500_ = bootstrap.CI.all(newdata = test_pat,
                               object_plac=rm_plac,
                               object_il12=rm_il12,
                               object_intg=rm_intg, 
                               object_tnfi=rm_tnfi,
                               nsim = 500, 
                               parallel = 'no', 
                               ncpus = 1), 
  nsim_1000 = bootstrap.CI.all(newdata = test_pat,
                               object_plac=rm_plac,
                               object_il12=rm_il12,
                               object_intg=rm_intg, 
                               object_tnfi=rm_tnfi,
                               nsim = 1000, 
                               parallel = 'no', 
                               ncpus = 1),
  times = 10
)

mbm.nsim
# Unit: seconds
# expr       min        lq      mean    median        uq       max neval
# nsim_100_  3.270677  3.287291  3.306178  3.302081  3.318815  3.366914    10
# nsim_250_  8.154116  8.199958  8.256431  8.229053  8.350920  8.398944    10
# nsim_500_ 16.339106 16.365023 16.414123 16.401052 16.429897 16.539326    10
# nsim_1000 32.625270 32.753363 32.832420 32.785416 32.884738 33.279996    10

