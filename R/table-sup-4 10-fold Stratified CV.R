#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# table-sup-4 10-fold Stratified CV.R
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(dplyr)
library(tidyr)
library(caret)
source('R/srs.R')
source("R/predict.R")
source("R/subgroups.R")
load("data/formula.RData")

#------------------------------------------------------------------------------#

srs.active(formula, data, srs.plac, treatment.field, exclude.cov) {
  
  srs.obj <- list()
  
  # Constants
  PLACEBO.RESP  <- as.character(attr(terms(formula), "variables")[[2]])
  TRT.LABELS    <- data %>% distinct( {{ treatment.field }} ) %>% pull( {{ treatment.field }} )
  
  # Infer Drug Attributable Effects
  data$ObsTotalEffect     <- data[[ PLACEBO.RESP ]]
  data$PredPlaceboEffect  <- predict(srs.plac, newdata=data, re.form=~0)
  data$InferredDrugEffect <- data[[ PLACEBO.RESP ]] - data$PredPlaceboEffect
  
  # Train Drug Attributable Models
  formula2 <- updateFormula(formula, "InferredDrugEffect", exclude.cov)
  for (trt in TRT.LABELS) {
    trt_df <- data %>% dplyr::filter( {{ treatment.field }} == trt )
    m <- lme4::lmer(formula2, data=trt_df)
    srs.obj$DrugModels[[trt]] <- m
  }
  
  return(srs.obj)
}

subgroup.template <- function() {
  template  <- data.frame(
    drug1  = numeric(0),
    drug2  = numeric(0),
    drug3  = numeric(0),
    p1.ohe = numeric(0),
    p2.ohe = numeric(0)
  )
  template[1,]  <- c("drug1","drug2","drug3","=","=")
  template[2,]  <- c("drug1","drug2","il12", "=",">")
  template[3,]  <- c("drug1","drug2","intg", "=",">")
  template[4,]  <- c("drug1","drug2","tnfi", "=",">")
  template[5,]  <- c("il12", "drug2","drug3",">","=")
  template[6,]  <- c("intg", "drug2","drug3",">","=")
  template[7,]  <- c("tnfi", "drug2","drug3",">","=")
  template[8,]  <- c("il12", "intg", "tnfi", ">",">")
  template[9,]  <- c("il12", "tnfi", "intg", ">",">")
  template[10,] <- c("intg", "il12", "tnfi", ">",">")
  template[11,] <- c("intg", "tnfi", "il12", ">",">")
  template[12,] <- c("tnfi", "intg", "il12", ">",">")
  template[13,] <- c("tnfi", "il12", "intg", ">",">")
  
  #template$p1.ohe <- as.numeric(template$p1.ohe)
  #template$p2.ohe <- as.numeric(template$p2.ohe)
  
  return(template)
}

mask.results <- function(data) {
  masked <- data %>% 
    # case when there is no drug preference (p1.ohe = 0, p2.ohe = 0)
    mutate(drug1 = ifelse(p1.ohe == "=" & p2.ohe == "=", 'drug1', drug1),
           drug2 = ifelse(p1.ohe == "=" & p2.ohe == "=", 'drug2', drug2), 
           drug3 = ifelse(p1.ohe == "=" & p2.ohe == "=", 'drug3', drug3)) %>% 
    
    # case when drug1 > drug2 = drug3 (p1.ohe = 1, p2.ohe = 0)
    mutate(drug2 = ifelse(p1.ohe == ">" & p2.ohe == "=", 'drug2', drug2),
           drug3 = ifelse(p1.ohe == ">" & p2.ohe == "=", 'drug3', drug3)) %>% 
    
    # case when drug1 = drug2 > drug3 (p1.ohe = 0, p2.ohe = 1)
    mutate(drug1 = ifelse(p1.ohe == "=" & p2.ohe == ">", 'drug1', drug1),
           drug2 = ifelse(p1.ohe == "=" & p2.ohe == ">", 'drug2', drug2))
  
  return(masked)
}

subgroup.cohorts <- function(data){
  # overview: find number of cases of 13 possible subgroups
  # 1.  tnfi > il12 > intg
  # 2.  il12 > > 
  # 3.  intg > > 
  # 4.  tnfi > drug2 = drug3
  # ...
  # 13. drug1 = drug2 = drug3
  
  subgroups <- mask.results(data) %>% 
    # collapse to 13 possible subgroups + add count
    group_by(drug1, drug2, drug3, p12_ohe, p23_ohe) %>% 
    summarise(n = n())
  
  return(subgroups)
}

srsStratifiedTestTrainSplit <- function(data, k=10, seed=1234) {
  if(!("strat_var" %in% colnames(data))) { 
    cat("Column strat_var does not exist in the dataframe.\n")
    return(NULL)
  }
  
  # Split into placebo and non-placebo recipients
  placebo_df <- data %>% filter(Treatment == "Placebo")
  active_df  <- data %>% filter(Treatment != "Placebo")
  
  # Train placebo model
  m1 <- lme4::lmer(formula, data=placebo_df)
  
  # Bootstrap / # k-fold stratified split 
  set.seed(seed)
  folds <- caret::createFolds(active_df$strat_var, k = k, list = FALSE)
  
  cohort_all  <- subgroup.template()
  cohort_var  <- subgroup.template()
  cols <- c("drug1","drug2","drug3","p1.ohe","p2.ohe")
  
  for(i in 1:k) {
    cat("Fold", i, "\n")
    crohnsTrain <- active_df[ folds != i, ]
    crohnsTest  <- active_df[ folds == i, ]
    
    m2 <- srs.active(formula, crohnsTrain, m1, Treatment, c("Year_Cent"))
    
    srs.cf <- batch.srs.predict(object.list = m2$DrugModels, newdata = crohnsTest, 
                                level = 0.95, interval = 'confidence', method = 'bootstrap',
                                nsim = 10000, parallel = 'multicore', ncpus = 8, seed = 1234)
    res <- srs.subgroups(srs.cf, m2$DrugModels)
    
    varname1 <- paste(paste("fold",i,sep=""))
    cohort_all <- cohort_all %>% left_join(., subgroup.cohorts(res), by = cols) %>% 
      rename(!!varname1 := n)
    cohort_var <- cohort_var %>% left_join(., subgroup.cohorts(res %>% filter(strat_var==1)), by=cols) %>% 
      rename(!!varname1 := n)
  }
  
  # cohort_all = normal cohort analysis
  # cohort_var = n participants flagged with strat_var==1 who fell into subgroup
  return(list("cohort_all"=cohort_all, "cohort_var"=cohort_var))
  
}

#------------------------------------------------------------------------------#

# IL12-23 Preference Subgroup (approx. female over 50) 
crohns_data1 <- read.csv("data/crohns_biologics_clean.csv") %>% 
  mutate(strat_var = ifelse( Age_Cent >= 15 & Sex_Male == 0, 1, 0))
# crohns_data1 %>% glimpse()
# 
# count(crohns_data1 %>% filter(strat_var == 1)) 
# # 577
# 
# count(crohns_data1 %>% filter(strat_var == 1, Group == 'Placebo')) 
# # 189
# 
# count(crohns_data1 %>% filter(strat_var == 1, Group == 'Active')) 
# # 388

#------------------------------------------------------------------------------#

FO50 <- srsStratifiedTestTrainSplit(crohns_data1, k=10, seed=1234)

num   <- colSums(FO50$cohort_var %>% dplyr::select(-c(drug1:p2.ohe)) %>% slice(5,8,9), na.rm=TRUE)
denom <- colSums(FO50$cohort_var %>% dplyr::select(-c(drug1:p2.ohe)), na.rm=TRUE)
a <- round( num / denom , 1 )
quantile(a)

num   <- colSums(FO50$cohort_var %>% dplyr::select(-c(drug1:p2.ohe)) %>% slice(5,8,9), na.rm=TRUE)
denom <- colSums(FO50$cohort_all %>% dplyr::select(-c(drug1:p2.ohe)) %>% slice(5,8,9), na.rm=TRUE)
b <- round( num / denom , 1 )
quantile(b)

saveRDS(FO50,"data/false_discovery_10fold_FO50.rds")

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
