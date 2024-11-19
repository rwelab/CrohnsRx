# Sequential Regression and Simulation (SRS)
Analytical code accompanying the manuscript "[Personalizing treatment selection in Crohn’s disease: a meta-analysis of individual participant data from fifteen randomized controlled trials](https://www.medrxiv.org/content/10.1101/2023.11.10.23291837v1)" by Rudrapatna et al., 2023

## How to Use
```r
# Load data
crohns <- read.csv("data/crohns_biologics_sample.csv")
# crohns %>% glimpse()
load("data/formula.RData")
# print(formula)

# Load functions
source("R/srs.R")
source("R/predict.R")
source("R/subgroups.R")

# Step1: Run Sequential Regression and Simulation (SRS)
srs.obj <- srs(formula         = formula, 
               data            = crohns, 
               treatment.field = Treatment, 
               placebo.label   = "Placebo",
               exclude.cov     = c("Year_Cent"))

# summary(srs.obj$PlaceboModel)
# summary(srs.obj$DrugModels$TNFi)
# summary(srs.obj$DrugModels$IL12)
# summary(srs.obj$DrugModels$Integrin)

# Step 2: Predict treatment counterfactuals (cf) + 95% CI using parametric bootstrapping
srs.cf <- batch.srs.predict(object.list = srs.obj$DrugModels, 
                            newdata     = crohns, 
                            level       = 0.95,
                            interval    = 'confidence', 
                            method      = 'bootstrap', 
                            nsim        = 100,
                            parallel    = "no",
                            ncpus       = 1,
                            seed        = 1234)
# srs.cf %>% glimpse()

# Step 3: Summarize subgroups
srs.sbgrps <- srs.summarise( srs.subgroups(srs.cf, srs.obj$DrugModels) )
knitr::kable(srs.sbgrps)

```
|Subgroup               | count|
|:----------------------|-----:|
|tnfi = integrin = il12 |   358|
|integrin = tnfi = il12 |   208|
|il12 = integrin = tnfi |   198|
|integrin = il12 = tnfi |    83|
|tnfi = il12 = integrin |    73|
|il12 = tnfi = integrin |    64|
|tnfi = integrin > il12 |     7|
|tnfi > integrin = il12 |     4|
|il12 = integrin > tnfi |     3|
|integrin = tnfi > il12 |     1|
|tnfi > integrin > il12 |     1|

## To-do List

* Package code
* Submit to CRAN

## Contact
