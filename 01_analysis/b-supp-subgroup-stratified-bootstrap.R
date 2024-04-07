#------------------------------------------------------------------------------#

library(dplyr)
library(tidyr)
library(caret)
options(dplyr.summarise.inform = FALSE) # suppress summarise info

setwd('UCSF-R/ipd-ma-cd2-2')
source('helper/srs.R') # srs.placebo(), srs.active(), subgroup.analysis()

#------------------------------------------------------------------------------#

subgroup.template <- function() {
  template  <- data.frame(
    drug1   = numeric(0),
    drug2   = numeric(0),
    drug3   = numeric(0),
    p12_ohe = numeric(0),
    p23_ohe = numeric(0)
  )
  template[1,]  <- c("drug1","drug2","drug3",0,0)
  template[2,]  <- c("drug1","drug2","il12", 0,1)
  template[3,]  <- c("drug1","drug2","intg", 0,1)
  template[4,]  <- c("drug1","drug2","tnfi", 0,1)
  template[5,]  <- c("il12", "drug2","drug3",1,0)
  template[6,]  <- c("intg", "drug2","drug3",1,0)
  template[7,]  <- c("tnfi", "drug2","drug3",1,0)
  template[8,]  <- c("il12", "intg", "tnfi", 1,1)
  template[9,]  <- c("il12", "tnfi", "intg", 1,1)
  template[10,] <- c("intg", "il12", "tnfi", 1,1)
  template[11,] <- c("intg", "tnfi", "il12", 1,1)
  template[12,] <- c("tnfi", "intg", "il12", 1,1)
  template[13,] <- c("tnfi", "il12", "intg", 1,1)
  
  template$p12_ohe <- as.numeric(template$p12_ohe)
  template$p23_ohe <- as.numeric(template$p23_ohe)
  
  return(template)
}

srsStratifiedTestTrainSplit <- function(data, k=10, seed=1234) {
  if(!("strat_var" %in% colnames(data))) { 
    cat("Column strat_var does not exist in the dataframe.\n")
    return(NULL)
  }
  
  # Split into placebo and non-placebo recipients
  placebo_df <- crohns_data1 %>% filter(Group == "Placebo")
  active_df  <- crohns_data1 %>% filter(Group == "Active")
  
  # Train placebo model
  m1 <- srs.placebo(placebo_df)
  
  # Bootstrap / # k-fold stratified split 
  set.seed(seed)
  folds <- caret::createFolds(active_df$strat_var, k = k, list = FALSE)
  
  cohort_all  <- subgroup.template()
  cohort_var  <- subgroup.template()
  cols <- c("drug1","drug2","drug3","p12_ohe","p23_ohe")
  
  for(i in 1:k) {
    cat("Fold", i, "\n")
    crohnsTrain <- active_df[ folds != i, ]
    crohnsTest  <- active_df[ folds == i, ]
    
    m2 <- srs.active(crohnsTrain, m1$plac_rm)
    res <- subgroup.analysis(crohnsTest, m1$plac_rm, m2$il12_rm, m2$intg_rm, m2$tnfi_rm)
    
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
crohns_data1 <- read.csv("data/crohns_data1.csv") %>% 
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

num   <- colSums(FO50$cohort_var %>% dplyr::select(-c(drug1:p23_ohe)) %>% slice(5,8,9), na.rm=TRUE)
denom <- colSums(FO50$cohort_var %>% dplyr::select(-c(drug1:p23_ohe)), na.rm=TRUE)
a <- round( num / denom , 1 )
quantile(a)

num   <- colSums(FO50$cohort_var %>% dplyr::select(-c(drug1:p23_ohe)) %>% slice(5,8,9), na.rm=TRUE)
denom <- colSums(FO50$cohort_all %>% dplyr::select(-c(drug1:p23_ohe)) %>% slice(5,8,9), na.rm=TRUE)
b <- round( num / denom , 1 )
quantile(b)

saveRDS(FO50,"data/false_discovery_10fold_FO50.rds")

#------------------------------------------------------------------------------#
