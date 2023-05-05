# Load the required libraries
library(shiny)
library(shinydashboard)

source('crohnsUI.R')
source('crohnsServer/predict.lme.R') # predict.lme.all()
source('crohnsServer/subgroup.R')    # patient.drug.preferences()

rm_plac <- readRDS('data/rm_plac.rds')
rm_il12 <- readRDS('data/rm_il12.rds')
rm_intg <- readRDS('data/rm_intg.rds')
rm_tnfi <- readRDS('data/rm_tnfi.rds')

# model degrees of freedom (required for two.sample.t.test())
df_list <- list("il12" = nrow(rm_il12@frame) - 10, #  587 - 10
                "intg" = nrow(rm_intg@frame) - 10, # 1818 - 10
                "tnfi" = nrow(rm_tnfi@frame) - 10) # 1677 - 10

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
  treatment as of 2020 or later (e.g. Risankizumab)."

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
# DEFINE UI
################################################################################

ui_disclaimer <- function() {

  column(width = 12,
         h2("Disclaimer"),
         h4(introductory.text.1),
         h4(introductory.text.2),
         h4(introductory.text.3),
         h4(introductory.text.4), 
         # button takes user to recommender tab
         actionButton("understand", "I Understand")
  )

}

ui_recommender <- function() {
  
  column(width = 12, 
         
         tags$hr(style="border-color: grey;"), # grey divider bar
         
         # Age, Sex
         fixedRow(
           column(4, textInput("age", 
                               value = NULL, 
                               label = 'Age', 
                               placeholder = "Age: 18-100")),
           column(4, selectizeInput("sex", 
                                    choices = c('Male / Female' = '', 'Male', 'Female'), 
                                    selected = NULL, 
                                    label = 'Sex'))
         ),
         
         # BMI, CRP
         fixedRow(
           column(4, textInput("bmi", 
                               value = NULL, 
                               label = 'BMI (kg/m2)', 
                               placeholder = "Norm: 18.5-24.9 kg/m2")),
           column(4, textInput("crp", 
                               value = NULL, 
                               label = 'C-Reactive Protein (mg/L)', 
                               placeholder = "Norm: <10 mg/L")),
           helpText(crp.text)
         ),
         
         tags$hr(style="border-color: grey;"),
         
         # History of TFN use, Current steroid use
         fixedRow(
           column(4, selectizeInput("tnf", 
                                    choices = c('No / Yes' = '', "No", "Yes"), 
                                    selected = NULL, 
                                    label = "History of Anti-Tumor Necrosis Factor Use")),
           column(4, selectizeInput("ste", 
                                    choices = c('No / Yes' = '', "No", "Yes"), 
                                    selected = NULL, 
                                    label = "Current Corticosteroid Use")),
           helpText(paste(tnf.text, steroid.text))
         ),
         
         # Current immunosuppressant use, Ileal involvement
         fixedRow(
           column(4, selectizeInput("imm", 
                                    choices = c('No / Yes' = '', "No", "Yes"), 
                                    selected = NULL, 
                                    label = "Current Immunosuppressant Use")),
           column(4, selectizeInput("loc", 
                                    choices = c('No / Yes' = '', "No", "Yes"), 
                                    selected = NULL, 
                                    label = "Ileal Involvement")),
           helpText(immuno.text)
         ),
         
         tags$hr(style="border-color: grey;"),
         
         # CDAI
         textInput('cdai', 
                   value = NULL, 
                   label = "Crohn's Disease Activity Index (CDAI)", 
                   placeholder = "CDAI: 0-600"),
         helpText(tags$div(
           "If not known, leave blank. Alternatively, manually calculate CDAI using this",
           tags$a(href="https://www.mdcalc.com/calc/3318/crohns-disease-activity-index-cdai", 
                  "online calculator.")
           )
         ),
         
         tags$hr(style="border-color: grey;"),
         
         fixedRow(
           column(4, actionButton("submit", "Submit"))
         ),
         
         tags$hr(style="border-color: grey;"),
         
         #--------------------------------------------------------------------------#
         
         # OUTPUTS
         textOutput('rtitle'),
         tags$head(tags$style("#rtitle{font-size: 32px; font-style: bold;}")),
         br(), 
         textOutput('results'),
         tags$head(tags$style("#results{font-size: 20px;}")),
         br(), 
         plotOutput('plot'),
  )
}

ui <- dashboardPage(
  
  # define the dashboard header
  dashboardHeader(title = "Treatment Recommender for Crohn's Disease (Beta)", 
                  titleWidth = 500),
  
  # define the dashboard sidebar
  dashboardSidebar(
    sidebarMenu(
      id = 'tabs',
      menuItem("Disclaimer", tabName = "disclaimer"),
      menuItem("Calculator", tabName = "recommender")
    )
  ),
  
  # define the dashboard body
  dashboardBody(
    # Sets entire app background to the same color
    tags$head(tags$style(HTML('
         .skin-blue .left-side, .skin-blue .wrapper {
                        background-color: #ecf0f5;
                        }
         '))), 
    
    tabItems(
      # Disclaimer tab
      tabItem(tabName = 'disclaimer', ui_disclaimer()),
      
      # Recommender tab
      tabItem(tabName = "recommender", ui_recommender())
    )
  )
)

################################################################################
# DEFINE SERVER
################################################################################

server <- function(input, output, session) {
  
  observeEvent(input$understand, {
    updateTabItems(session, "tabs", selected = "recommender")
  })
  
  # extract and preprocess patient input data
  patient_data <- eventReactive(input$submit, {
    # error handling - ensure all inputs are entered
    validate(
      need(input$age, "Please enter patient's age."),
      need(input$sex, "Please select patient's gender."),
      need(input$bmi, "Please enter patient's BMI."),
      need(input$crp, "Please enter patient's latest c-reactive protein (mg/L) lab result values. If not known, enter `0`."),
      need(input$tnf, "Please select if your patient has taken anti-tumor necrosis factor drugs in the past."),
      need(input$ste, "Please select if your patient is currently on a steroid."),
      need(input$imm, "Please select if your patient is currently on an immunosuppressant"),
      need(input$loc, "Please select if your patient's disease location is in the ileum.")
    )
    
    # crohnsUI::ReadUI, crohnsUI::Center
    Center( ReadUI(input) )
  })
  
  # patient drug class preferences
  preferences <- eventReactive(input$submit, {
    
    # predict drug class attributable (attrib) and standard error of mean
    # response (se) using analytical solution (faster than bootstrapping)
    data <- predict.lme.all(object.list = list(rm_plac, rm_il12, rm_intg, rm_tnfi), 
                            newdata = patient_data(), 
                            model.name.list = c('plac','il12','intg','tnfi'), 
                            interval = 'confidence', method = 'analytical')
    
    # find drug class preferences based on predicted drug class attributable
    # effects and SEs (CI)
    result <- patient.drug.preferences( data , df_list ) %>% 
      dplyr::select(drug1:drug3, p12_ohe, p23_ohe)
    
    result
  })
  
  # probability of drug inducing clinical response (CDAI reduction >100 points)
  response <- eventReactive(input$submit, {
    
    # predict drug class attributable (attrib) and standard error of prediction (se)
    # using analytical solution
    data <- predict.lme.all(object.list = list(rm_plac, rm_il12, rm_intg, rm_tnfi), 
                            newdata = patient_data(), 
                            model.name.list = c('plac','il12','intg','tnfi'), 
                            interval = 'prediction', method = 'analytical')
    
    # probability of reaching clinical response
    result <- data %>% mutate(
      il12.response = 1 - pnorm(100, mean = plac.attrib+il12.attrib, sd = il12.se),
      intg.response = 1 - pnorm(100, mean = plac.attrib+intg.attrib, sd = intg.se),
      tnfi.response = 1 - pnorm(100, mean = plac.attrib+tnfi.attrib, sd = tnfi.se)) %>% 
      dplyr::select(il12.response:tnfi.response)
      
    result
  })

  # generate outputs
  ## more error handling - prevent error messages from being displayed in output section
  valid_input <- eventReactive(input$submit, {
    if (input$age == "" | input$sex == "" | input$bmi == "" | input$crp == "" |
        input$tnf == "" | input$ste == "" | input$imm == "" | input$loc == "") {
      return(FALSE)
    } else {
      return(TRUE)
    }
  })
  
  rtitle <- eventReactive(input$submit, {
    # ensure title isn't printed before submit button is pressed
    "Treatment Recommendation Results:"
  })

  script <- eventReactive(input$submit, {
    # crohnsUI::GenerateScripts
    GenerateScripts( preferences(), response() )
  })
  
  output$rtitle  <- renderText( ifelse(!valid_input(), "", rtitle()) )
  output$results <- renderText( ifelse(!valid_input(), "", paste( script() )) )
  output$plot    <- renderPlot( GenereateResultPlot( response() ) )
}

################################################################################
# LAUNCH APP
################################################################################

shinyApp(ui, server)

################################################################################