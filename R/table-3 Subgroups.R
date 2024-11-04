#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# table-3 Subgroups.R
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Load functions
source("R/srs.R")
source("R/predict.R")
source("R/subgroups.R")

#------------------------------------------------------------------------------#

# Load data
crohns <- read.csv("data/crohns_biologics_clean.csv")
load("data/formula.RData")

# Train placebo, drug attributable models
srs.obj <- srs(formula         = formula, 
               data            = crohns, 
               treatment.field = Treatment, 
               placebo.label   = "Placebo",
               exclude.cov     = c("Year_Cent"))

# Find Subgroups
srs.cf <- batch.srs.predict(object.list = srs.obj$DrugModels, 
                            newdata     = crohns, 
                            level       = 0.95, 
                            interval    = 'confidence', 
                            method      = 'bootstrap', 
                            nsim        = 10000, 
                            parallel    = 'multicore', 
                            ncpus       = 8, 
                            seed        = 1234)

srs.sbgrps <- srs.summarise( srs.subgroups(srs.cf, srs.obj$DrugModels) )
knitr::kable(srs.sbgrps)

# |Subgroup               | count|
# |:----------------------|-----:|
# |tnfi > integrin = il12 |  1631|
# |tnfi = il12 = integrin |  1511|
# |tnfi = integrin = il12 |   709|
# |il12 = tnfi = integrin |   705|
# |tnfi > il12 = integrin |   389|
# |il12 = tnfi > integrin |   205|
# |il12 = integrin = tnfi |   189|
# |tnfi = il12 > integrin |   149|
# |il12 > tnfi = integrin |    87|
# |il12 > integrin = tnfi |    51|
# |tnfi > integrin > il12 |    43|
# |integrin = tnfi = il12 |    19|
# |integrin = il12 = tnfi |    10|
# |tnfi = integrin > il12 |     4|
# |il12 > tnfi > integrin |     1|

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#