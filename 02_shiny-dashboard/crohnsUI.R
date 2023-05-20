library(dplyr)
library(ggplot2)

#------------------------------------------------------------------------------#

# reads covariate values from ui
ReadUI <- function(input) {
  data <- data.frame(
    # reading text input: missing values are ""
    CDAI       = input$cdai, 
    Age        = input$age, 
    BMI        = input$bmi, 
    CRP        = input$crp, 
    Sex        = input$sex, 
    HxOfTNFi   = input$tnf, 
    SteroidUse = input$ste, 
    ImmUse     = input$imm, 
    Ileal      = input$loc
  )
  data
}

#------------------------------------------------------------------------------#

Center <- function(raw_data) {
  bin.map <- c("Yes"=1, "No"=0)
  sex.map <- c("Female"=0, "Male"=1)
  
  # account for missing input data
  # reading text input: missing values are ""
  raw_data$CDAI <- ifelse(raw_data$CDAI == '', 300, raw_data$CDAI)
  raw_data$Age  <- ifelse(raw_data$Age  == '', 35, raw_data$Age)
  raw_data$BMI  <- ifelse(raw_data$BMI  == '', 20, raw_data$BMI)
  raw_data$CRP  <- ifelse(raw_data$CRP  == '', 10, raw_data$CRP)
  raw_data$HxOfTFNi   <- ifelse(raw_data$HxOfTNFi   == '', 'No', raw_data$HxOfTNFi)
  raw_data$Sex        <- ifelse(raw_data$Sex        == '', 'Female', raw_data$Sex)
  raw_data$SteroidUse <- ifelse(raw_data$SteroidUse == '', 'No', raw_data$SteroidUse)
  raw_data$ImmUse     <- ifelse(raw_data$ImmUse     == '', 'No', raw_data$ImmUse)
  raw_data$Ileal      <- ifelse(raw_data$Ileal      == '', 'No', raw_data$Ileal)
  
  cent_data <- data.frame(
    Year_Cent          = 6, # average year 2006
    CDAI_baseline_Cent = as.numeric(raw_data['CDAI']) - 300,
    Age_Cent           = as.numeric(raw_data['Age']) - 35,
    BMI_Cent           = as.numeric(raw_data['BMI']) - 20, 
    CRP_Cent           = as.numeric(raw_data['CRP']) - 10, 
    HxOfTNFi           = bin.map[as.character(raw_data['HxOfTNFi'])], 
    Sex_Male           = sex.map[as.character(raw_data['Sex'])], 
    SteroidUse         = bin.map[as.character(raw_data['SteroidUse'])],
    ImmUse             = bin.map[as.character(raw_data['ImmUse'])], 
    Ileal              = bin.map[as.character(raw_data['Ileal'])] 
  )
  
  print(cent_data)
  
  return(cent_data)
}

#------------------------------------------------------------------------------#

# output drug class recommendation script
GenerateScripts <- function(ranking, response) {
  
  titles <- c("il12" = 'Anti-Interleukin (IL)-12/23', 
              "intg" = 'Anti-Integrin', 
              "tnfi" = 'Anti-Tumor Necrosis Factor (TNF)')
  
  probs <- c('il12' = as.integer(100*response$il12.response), 
             'intg' = as.integer(100*response$intg.response), 
             'tnfi' = as.integer(100*response$tnfi.response))
  
  drug_recommendations <- c(titles[ranking$drug1], titles[ranking$drug2], titles[ranking$drug3])
  probability <- c(probs[ranking$drug1], probs[ranking$drug2], probs[ranking$drug3])
  
  if (ranking$p12_ohe == 0 & ranking$p23_ohe == 0) {
    
    prob1 <- as.integer(mean(probability))
    
    recommendation_text <- sprintf("Your patient is predicted as having greatest 
                                   efficacy with %s, %s, or %s treatment, with a 
                                   %d%% probability of reaching 
                                   clinical response. 
                                   Clinical response is defined as CDAI reduction 
                                   of 100 or more points after 6 weeks of treatment.", 
                                   drug_recommendations[1], drug_recommendations[2], drug_recommendations[3], 
                                   prob1)
    
  } else if (ranking$p12_ohe == 0 & ranking$p23_ohe == 1) {
    
    prob1 <- as.integer(mean(probability[1:2]))
    
    recommendation_text <- sprintf("Your patient is predicted as having greatest 
                                   efficacy with %s or %s treatment, with a 
                                   %d%% probability of reaching 
                                   clinical response. 
                                   Clinical response is defined as CDAI reduction 
                                   of 100 or more points after 6 weeks of treatment.", 
                                   drug_recommendations[1], drug_recommendations[2], 
                                   prob1)
  
  } else if (ranking$p12_ohe == 1) {
    
    recommendation_text <- sprintf("Your patient is predicted as having greatest 
                                   efficacy with %s treatment, with a 
                                   %d%% probability of reaching 
                                   clinical response. 
                                   Clinical response is defined as CDAI reduction 
                                   of 100 or more points after 6 weeks of treatment.", 
                                   drug_recommendations[1], 
                                   probability[1])
    
  }
  
  print("Generated scripts")
  
  return( recommendation_text )
}

#------------------------------------------------------------------------------#

GenereateResultPlot <- function(ranking, response) {
  
  titles <- c("il12" = 'Anti-Interleukin (IL)-12/23', 
              "intg" = 'Anti-Integrin', 
              "tnfi" = 'Anti-Tumor Necrosis Factor (TNF)')
  
  probs <- c('il12' = as.integer(100*response$il12.response), 
             'intg' = as.integer(100*response$intg.response), 
             'tnfi' = as.integer(100*response$tnfi.response))
  
  drug_recommendation <- c(titles[ranking$drug1], titles[ranking$drug2], titles[ranking$drug3])
  probability <- c(probs[ranking$drug1], probs[ranking$drug2], probs[ranking$drug3])
  
  if(ranking$p12_ohe == 0 & ranking$p23_ohe == 0) {
    
    dat3 <- data.frame(
      DrugClass = c(paste(drug_recommendation[1], drug_recommendation[2], drug_recommendation[3], sep = ' or ')), 
      Response  = c(mean(probability[1:3]))
    )
    
  } else if(ranking$p12_ohe == 0 & ranking$p23_ohe == 1) {
    
    dat3 <- data.frame(
      DrugClass = c(paste(drug_recommendation[1], drug_recommendation[2], sep = ' or '), 
                    drug_recommendation[3]), 
      Response  = c(mean(probability[1:2]), 
                    probability[3])
    )
    
  } else if(ranking$p12_ohe == 1 & ranking$p23_ohe == 0) {
    
    dat3 <- data.frame(
      DrugClass = c(drug_recommendation[1], 
                    paste(drug_recommendation[2], drug_recommendation[3], sep = ' or ')), 
      Response  = c(probability[1], 
                    mean(probability[2:3]))
    )
    
  } else if(ranking$p12_ohe == 1 & ranking$p23_ohe == 1) {
    dat3 <- data.frame(
      DrugClass = c(drug_recommendation[1], drug_recommendation[2], drug_recommendation[3]),
      Response  = c(probability[1], probability[2], probability[3])
    )
  }
  
  responseColor <- "#00AFBB"
    
  p <- ggplot(dat3, aes(x=DrugClass)) +
    
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

