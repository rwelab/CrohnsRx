# subgroup.R
# 
# Two tasks are performed:
# 1) The drug class order (from best performing:drug1 to worst performing:drug3)
#    for the three modeled drug classes - interleukin-12/23 (il12), integrin (intg), 
#    and tumor necrosis factor-alpha (tnfi) - are determined for each participant
#    (row). 
# 2) Two-sample t-tests are applied to drug pairs (drug1 vs drug2 and drug2 vs drug3)
#    to determine if a participants significantly prefers one or more drug classes
#    over another.

#------------------------------------------------------------------------------#

library(dplyr)

#------------------------------------------------------------------------------#

patient.drug.preferences <- function(data, df_list) {
  # order drug classes by most to least predicted CDAI reduction (attrib) (drug1, drug2, drug3)
  data_ranked     <- rank.drugs(data)
  # assign drug class superiority (drug1 > drug2 = drug3)
  data_preference <- preferences(data_ranked, df_list) 
}

rank.drugs <- function(data) {
  
  df <- data.frame(
    tnfi = data$tnfi.attrib, 
    il12 = data$il12.attrib, 
    intg = data$intg.attrib
  )
  
  # example df
  ## tnfi  il12  intg
  ##   30    60    45 <- correspond to drug attributable cdai reduction for patient
  ##  ...
  
  ranked <- data.frame(t(apply(-df, 1, rank, ties.method='min')) ) %>% 
    # rank each outcome by row, add patient id (for grouping)
    ## id tnfi  il12  intg
    ##  1    3     1     2
    ##  ...
    mutate(id = row_number()) %>% 
    relocate(id) %>%
    
    # pivot longer (rank based on drug class attributable magnitude)
    ## id drug  rank
    ##  1 tnfi     3
    ##  1 il12     1
    ##  1 intg     2
    ##  ...
    pivot_longer(cols = c(tnfi, il12, intg), names_to = "drug", values_to = "rank") %>% 
    
    # sort drug classes from rank 1 to 3
    ## id drug  rank
    ##  1 il12     1
    ##  1 intg     2
    ##  1 tnfi     3
    ##  ...
    arrange(id, rank) %>% 
    
    # create drug1, drug2, drug3 columns corresponding to rank 1, 2, 3 drug classes
    ## id drug1  drug2  drug3
    ##  1  il12   intg   tnfi
    ##  ...
    pivot_wider(names_from = rank, values_from = drug, names_prefix = "drug") %>%
    ungroup() %>% 
    dplyr::select(-id)
  
  return(data.frame(data, ranked))
}

preferences <- function(data, df_list) {
  
  # calculate p-values between drug1:drug2 and drug2:drug3
  data['p12'] <- apply(X = data, 
                       FUN = two.sample.t.tests, 
                       MARGIN = 1,                # row-wise
                       args = c('drug1','drug2'),
                       df_list = df_list)
  
  data['p23'] <- apply(X = data, 
                       FUN = two.sample.t.tests, 
                       MARGIN = 1,                # row-wise
                       args = c('drug2','drug3'),
                       df_list = df_list)
  
  # one-hot encode (ohe) p12 and p23
  result <- data %>% 
    mutate(
      p12_ohe = ifelse(p12 < 0.05, 1, 0),
      p23_ohe = ifelse(p23 < 0.05, 1, 0)
    )
  
  return(result)
}

two.sample.t.tests <- function(data, args, df_list){
  # data = data.frame of model covariates
  # args = c('drug1','drug2') or c('drug2','drug3')
  # function calculates p-value of a t-score between two drug classes (args)
  
  # initialize 
  X  <- c() # drug class fit
  SE <- c() # drug class standard error
  df <- 0   # sum of drug class degrees of freedom
  
  # for args (drug1, drug2, and/or drug3)
  for(drug in args){
    X  <- c(X,  data[[paste0(data[[drug]], '.attrib')]]) # ex. tnfi.attrib
    SE <- c(SE, data[[paste0(data[[drug]], '.se')]])     # ex. tnfi.se
    df <- df + df_list[[data[[drug]]]]
  }
  
  # convert vectors to numeric
  X <- as.numeric(X)
  SE <- as.numeric(SE)
  
  # Calculate p-value (pt)
  # p = 2 * pt( abs(X1 - X2) / sqrt(SE1^2 + SE2^2) )
  # df = df_X1 + df_X2
  p.value <- 2 * pt(abs(X[1] - X[2]) / sqrt(SE[1]^2 + SE[2]^2),
                    df = df,
                    lower.tail = F)
  
  return(p.value)
}

#------------------------------------------------------------------------------#

# mask drugs where drug order does not matter (p == 0) - makes grouping easier
# ex. tnfi  > drug2 = drug3 
# ex. drug1 = drug2 = drug3
mask.results <- function(data) {
  masked <- data %>% 
    # case when there is no drug preference (p12_ohe = 0, p23_ohe = 0)
    mutate(drug1 = ifelse(p12_ohe == 0 & p23_ohe == 0, 'drug1', drug1),
           drug2 = ifelse(p12_ohe == 0 & p23_ohe == 0, 'drug2', drug2), 
           drug3 = ifelse(p12_ohe == 0 & p23_ohe == 0, 'drug3', drug3)) %>% 
    
    # case when drug1 > drug2 = drug3 (p12_ohe = 1, p23_ohe = 0)
    mutate(drug2 = ifelse(p12_ohe == 1 & p23_ohe == 0, 'drug2', drug2),
           drug3 = ifelse(p12_ohe == 1 & p23_ohe == 0, 'drug3', drug3)) %>% 
    
    # case when drug1 = drug2 > drug3 (p12_ohe = 0, p23_ohe = 1)
    mutate(drug1 = ifelse(p12_ohe == 0 & p23_ohe == 1, 'drug1', drug1),
           drug2 = ifelse(p12_ohe == 0 & p23_ohe == 1, 'drug2', drug2))
  
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

#------------------------------------------------------------------------------#