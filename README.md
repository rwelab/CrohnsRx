# Sequential Regression and Simulation (SRS)
[Analytical code accompanying the manuscript "Personalizing treatment selection in Crohn’s disease: a meta-analysis of individual participant data from fifteen randomized controlled trials" by Rudrapatna et al., 2023](https://www.medrxiv.org/content/10.1101/2023.11.10.23291837v1)

## How to Use
```r
# Load data
crohns <- read.csv("data/crohns_biologics_sample.csv")
# crohns %>% glimpse()

# Load functions
source("R/srs.R")
source("R/predict.R")
source("R/subgroups.R")

# Step 1 : Train placebo attributable srs model
srs.plac <- srs.placebo(data               = crohns, 
                        rand.var           = Trial,
                        treatment          = Treatment,
                        outcome            = CDAI_reduction,
                        placebo.label      = "Placebo")

# Step 2 : Train drug attributable srs models
srs.drug <- srs.active(data                = crohns, 
                       srs.plac            = srs.plac, 
                       rand.var            = Trial, 
                       treatment           = Treatment,
                       outcome             = CDAI_reduction,
                       placebo.label       = "Placebo",
                       exclude.covariates  = c("Year_Cent"))

# summary(srs.plac)
# summary(srs.drug$tnfi)
# summary(srs.drug$il12)
# summary(srs.drug$integrin)

# Step 3 : Predict treatment counterfactuals (cf) using srs models
#
# Predicting drug attributable effects and 95% CI using parametric bootstrapping
# (nsim = 10k). Leveraging parallel processing (ncpus=8); ~80 sec compute.
srs.cf <- batch.srs.predict(object.list = srs.drug, 
                            newdata     = crohns, 
                            level       = 0.95,
                            interval    = 'confidence', 
                            method      = 'bootstrap', 
                            nsim        = 10000, 
                            parallel    = 'multicore', 
                            ncpus       = 8, 
                            seed        = 1234)

# Step 4 : Summarise subgroups
srs.sbgrps <- srs.summarise( srs.subgroups(srs.cf, srs.drug) )
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
