library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(merTools)
library(parallel)

rm_plac <- readRDS('data/rm_plac.rds')
rm_il12 <- readRDS('data/rm_il12.rds')
rm_intg <- readRDS('data/rm_intg.rds')
rm_tnfi <- readRDS('data/rm_tnfi.rds')

################################################################################
# DISPLAY TEXT
################################################################################

# Disclaimer text
introductory.text.1 <- tags$div(
  "This clinical decision support prototype is based on the work of",
  tags$a(href="https://www.medrxiv.org/content/10.1101/2022.10.07.22280801v2", "Rudrapatna and colleagues"), 
  ", who developed and reported a new model for predicting treatment efficacy in 
  Crohn’s disease based on patient-level characteristics.")

introductory.text.2 <- tags$div("
  These models utilized the raw data from randomized trials in Crohn’s disease, 
  and thus they are only suitable for predictions in real-world patients meeting 
  typical eligibility criteria for registrational trials in Crohn’s. These 
  include active and symptomatic Crohn’s disease involving the ileum and/or colon, 
  and the absence of symptomatic obstructions or untreated infections. See",
                                tags$a(href='https://clinicaltrials.gov/ct2/show/NCT01369342?term=UNITI&cond=Crohn%27s+disease&draw=2&rank=3#eligibility', 
                                       "here"), "for an example of typical inclusion/exclusion criteria.")

introductory.text.3 <- "
  These models were subject to multiple internal validation assessments and 
  sensitivity analyses. However, the model predictions have not yet been 
  prospectively validated. Clinicians desiring to use this tool to inform patient 
  care are advised to exercise an appropriate level of caution. Please note that 
  this model does not yet incorporate recommendations for any FDA-approved 
  treatment as of 2020 or later (e.g. Risankizumab, Upadacitinib)."

introductory.text.4 <- "
  Users should also be aware that this model is only designed to make recommendations 
  pertaining to treatment with an FDA-approved biologic, and should not be used 
  to guide decision making pertaining to steroid or immunomodulator use. For 
  example, if a patient’s predicted probability of clinical response with a given 
  biologic class is higher with a different choice of baseline steroid/immunomodulator 
  use than what the patient is already receiving as part of the standard of care, 
  this does not imply that they will be expected to improve with a switch in these 
  concomitant medications."

# Recommender text
crp.text     <- 'C-reactive protein (CRP) can be reported in different units (mg/dL or mg/L) depending on the clinical laboratory. Please check units before entering or leave blank if a recent CRP is not available.'
tnf.text     <- 'Anti-TNFs include adalimumab, infliximab, or certolizumab pegol.'
steroid.text <- 'Corticosteroids include methylprednisolone, prednisolone, 
                   hydrocortisone, budesonide, or beclomethasone.'
immuno.text  <- 'Immunosuppressants include azathioprine, methotrexate, mercaptopurine, 
                   or an equivalent drug.'

################################################################################
# READ USER INPUTS
################################################################################

# read user inputs
getInput <- function(input) {
  can_afford <- paste(input$q_cost, collapse="|")
  
  raw_data <- data.frame(
    # reading text input: missing values are ""
    CDAI       = ifelse(input$cdai == '', 300, input$cdai), 
    Age        = ifelse(input$age  == '', 35, input$age), 
    BMI        = ifelse(input$bmi  == '', 20, input$bmi), 
    CRP        = ifelse(input$crp  == '', 10, input$crp), 
    Sex        = ifelse(input$sex  == '', 'Female', input$sex), 
    HxOfTNFi   = ifelse(input$tnf  == '', 'No', input$tnf), 
    SteroidUse = ifelse(input$ste  == '', 'No', input$ste), 
    ImmUse     = ifelse(input$imm  == '', 'No', input$imm), 
    Ileal      = ifelse(input$loc  == '', 'No', input$loc),
    # Contraindications
    con_tnfi    = ifelse(input$q_tnfi  == '', 'No', input$q_tnfi),
    con_il12    = ifelse(input$q_il12  == '', 'No', input$q_il12),
    con_intg    = ifelse(input$q_intg  == '', 'No', input$q_intg),
    afford_tnfi = ifelse(grepl("Anti-TNFs", can_afford), 'No', 'Yes'),
    afford_il12 = ifelse(grepl("Anti-IL12/23s", can_afford), 'No', 'Yes'),
    afford_intg = ifelse(grepl("Anti-integrins", can_afford), 'No', 'Yes')
  )
  raw_data
}

# serves as a main function
transformData <- function(raw) {
  data <- centerData(raw)
  pred <- getDrugPreference(data)
  resp <- getDrugResponse(data)
  final_data <- c(data, pred, resp)
  return(final_data)
}

# center data
centerData <- function(raw_data) {
  bin.map <- c("Yes"=1, "No"=0)
  sex.map <- c("Female"=0, "Male"=1)
  
  data <- data.frame(
    Year_Cent          = 6, # average year 2006
    CDAI_baseline_Cent = as.numeric(raw_data['CDAI']) - 300,
    Age_Cent           = as.numeric(raw_data['Age']) - 35,
    BMI_Cent           = as.numeric(raw_data['BMI']) - 20,
    CRP_Cent           = as.numeric(raw_data['CRP']) - 10,
    HxOfTNFi           = bin.map[as.character(raw_data['HxOfTNFi'])],
    Sex_Male           = sex.map[as.character(raw_data['Sex'])],
    SteroidUse         = bin.map[as.character(raw_data['SteroidUse'])],
    ImmUse             = bin.map[as.character(raw_data['ImmUse'])],
    Ileal              = bin.map[as.character(raw_data['Ileal'])],
    # omit drug classes
    omit_tnfi          = as.numeric(ifelse(raw_data['con_tnfi']=='Yes' | raw_data['afford_tnfi']=='No', 1, 0)),
    omit_il12          = as.numeric(ifelse(raw_data['con_il12']=='Yes' | raw_data['afford_il12']=='No', 1, 0)), 
    omit_intg          = as.numeric(ifelse(raw_data['con_intg']=='Yes' | raw_data['afford_intg']=='No', 1, 0))
  )
  row.names(data) <- NULL
  data
}

# predict drug preferences
getDrugPreference <- function(data) {
  # predict drug class attributable (attrib) and standard error of mean
  # response (se) using analytical solution (faster than bootstrapping)
  pred <- predict.lme.all(object.list = list(rm_plac, rm_il12, rm_intg, rm_tnfi),
                          newdata = data,
                          model.name.list = c('plac','il12','intg','tnfi'),
                          interval = 'confidence', method = 'analytical')
  
  # find drug class preferences based on predicted drug class attributable 
  # effects and SEs (CI)
  df_list <- list("il12" = nrow(rm_il12@frame) - 10, #  587 - 10
                  "intg" = nrow(rm_intg@frame) - 10, # 1818 - 10
                  "tnfi" = nrow(rm_tnfi@frame) - 10) # 1677 - 10
  
  pred <- patient.drug.preferences( pred , df_list ) %>%
    dplyr::select(drug1:drug3, p12_ohe, p23_ohe)
}

# get drug responses
getDrugResponse <- function(data) {
  # predict drug class attributable (attrib) and standard error of prediction (se)
  # using analytical solution
  resp <- predict.lme.all(object.list = list(rm_plac, rm_il12, rm_intg, rm_tnfi),
                          newdata = data,
                          model.name.list = c('plac','il12','intg','tnfi'),
                          interval = 'prediction', method = 'analytical')
  
  # probability of reaching clinical response
  resp <- resp %>% mutate(
    il12.response = 1 - pnorm(100, mean = plac.attrib+il12.attrib, sd = il12.se),
    intg.response = 1 - pnorm(100, mean = plac.attrib+intg.attrib, sd = intg.se),
    tnfi.response = 1 - pnorm(100, mean = plac.attrib+tnfi.attrib, sd = tnfi.se)) %>%
    dplyr::select(il12.response:tnfi.response)
}

################################################################################
# FORMAT USER INPUTS
################################################################################

formatDrugPref <- function(data) {
  # Create a pattern: ex. "D1 > D2 = D3" 
  pattern <- paste(data$drug1, ifelse(data$p12_ohe==1,">","="), data$drug2, ifelse(data$p23_ohe==1,">","="), data$drug3)
  
  # Specify which drug classes to omit
  omt.map <- c('il12' = data$omit_il12, 'intg' = data$omit_intg, 'tnfi' = data$omit_tnfi)
  omit <- names(omt.map)[omt.map == 1]
  
  # Split pattern by ">" and "=" symbols: ex.["D1", "D2 = D3"]
  groups <- strsplit(pattern, " > ")[[1]]
  
  title.map <- c('il12' = 'Anti-Interleukin (IL)-12/23', 'intg' = 'Anti-Integrin', 'tnfi' = 'Anti-Tumor Necrosis Factor (TNF)')
  resp.map <- c('il12' = as.integer(100*data$il12.response), 'intg' = as.integer(100*data$intg.response), 'tnfi' = as.integer(100*data$tnfi.response))
  
  res <- list()
  for (i in 1:length(groups)) {
    idx <- setdiff(unlist(strsplit(groups[i], " = ")), omit)    # ex. group[i] - omit
    recommendation <- paste(title.map[idx], collapse = " or ")
    recommendation <- paste(unlist(recommendation), collapse='')
    response <- mean(resp.map[idx])
    if(recommendation != "") { res <- rbind(res, c(recommendation, response)) }
  }
  
  if( length(res)==0 ){
    return(0)
  }
  
  res <- as.data.frame(res)
  colnames(res) <- c("Recommendation","Response")
  res[res == "NaN"] <- NA
  res[res == ""] <- NA
  res$Response <- as.numeric(res$Response)
  res$Recommendation <- as.character(res$Recommendation)
  
  return(res)
  
  # # Populate categories based on the parsed groups
  # First <- Second <- Third <- character(0)
  # if (length(groups) >= 1) { First  <- setdiff(unlist(strsplit(groups[1], " = ")), omit) } # ex. D1 - omit
  # if (length(groups) >= 2) { Second <- setdiff(unlist(strsplit(groups[2], " = ")), omit) } # ex. D2, D3 - omit
  # if (length(groups) >= 3) { Third  <- setdiff(unlist(strsplit(groups[3], " = ")), omit) } # ex. character(0) - omit
  # res <- list(First = First, Second = Second, Third = Third)
}

# data <- data.frame(
#   drug1 = 'tnfi',
#   drug2 = 'intg',
#   drug3 = 'il12',
#   p12_ohe = 1,
#   p23_ohe = 0,
#   tnfi.response = 0.96,
#   intg.response = 0.84,
#   il12.response = 0.76,
#   omit_tnfi = 0,
#   omit_il12 = 0,
#   omit_intg = 0
# )

################################################################################
# GENERATE RESULT SCRIPT
################################################################################

# output drug class recommendation script
GenerateScripts <- function(data) {
  if ( identical(data,0) ) {
    script <- "Based on your responses, none of the following three drug classes (Anti-Tumor Necrosis Factor (TNF), Anti-Integrin, Anti-Interleukin (IL)-12/23) would be suitable for your patient. Consider other options (e.g. anti-JAKs, surgery, immunomodulators, nutritional modification, steroids)."
  } else {
    firstLine = data[1,]
    script <- sprintf("Your patient is predicted as having greatest efficacy with %s, with a %0.2f%% probability of reaching clinical response. Clinical response is defined as CDAI reduction of 100 or more points after 6 weeks of treatment.", 
                      firstLine[[1]], firstLine[[2]]) 
  }
  return(script)
}

################################################################################
# GENERATE RESULT PLOT
################################################################################

GenerateResultPlot <- function(data) {
  
  if ( identical(data,0) ) { return(0) }
  
  responseColor <- "#00AFBB"
  p <- ggplot(data, aes(x=Recommendation)) +
    # effectiveness bar plot
    geom_bar(aes(y=Response), 
             stat='identity', width = 0.3, 
             alpha = 0.8, fill=responseColor) + 
    theme_bw() + 
    theme(panel.grid.major.x = element_blank()) + 
    ggtitle("Drug Class Clinical Response Rates") + 
    ylab('Response Rate') +
    xlab('Drug Class') + 
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100))
  
  return(p)
}

################################################################################
# ABBV. PREDICT FUNCTIONS (ANALYTICAL ONLY)
################################################################################

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

# REF: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6916346/#sim8386-sec-0020
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

################################################################################
# ABBV. SUBGROUP FUNCTIONS
################################################################################

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
    tidyr::pivot_longer(cols = c(tnfi, il12, intg), names_to = "drug", values_to = "rank") %>% 
    
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
    tidyr::pivot_wider(names_from = rank, values_from = drug, names_prefix = "drug") %>%
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
