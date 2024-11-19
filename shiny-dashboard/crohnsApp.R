# Load the required libraries
library(shiny)
library(shinydashboard)
source('global.R')

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
         
         fluidRow(
           column(8, "Does your patient have any contraindications to anti-TNFs (e.g. increased risk of cancer, serious infections, heart failure, or demyelinating disorders)?"),
           column(4, selectizeInput("q_tnfi", label = NULL, choices = c('No / Yes' = '', "No", "Yes")))
         ),
         
         fluidRow(
           column(8, "Does your patient have any contraindications to anti-IL12/23s (e.g. increased risk of cancer, or of serious infections)?"),
           column(4, selectizeInput("q_il12", label = NULL, choices = c('No / Yes' = '', "No", "Yes")))
         ),
         
         fluidRow(
           column(8, "Does your patient have any contraindications to anti-integrins?"),
           column(4, selectizeInput("q_intg", label = NULL, choices = c('No / Yes' = '', "No", "Yes")))
         ),
         
         fluidRow(
           column(8, "Is your patient unable to access or afford any of these drug classes?"),
           column(4, checkboxGroupInput("q_cost", label = NULL,
                                        choices = list("Anti-TNFs", "Anti-IL12/23s", "Anti-integrins")))
         ),
         
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
         br(),
         plotOutput('plot')
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
  
  # 'I Understand' moves to recommender (Calculator) tab
  observeEvent(input$understand, {
    updateTabItems(session, "tabs", selected = "recommender")
  })
  
  # Preprocess user inputs, predict cohort and drug response
  data <- eventReactive(input$submit, {
    validate(
      need(input$age != '', 'Please enter your age.'),
      need((as.numeric(input$age) >= 18 && as.numeric(input$age) <= 99), "Invalid age. Age must be between 18 and 99."),
      need((as.numeric(input$bmi) >= 12 && as.numeric(input$bmi) <= 45) || (input$bmi == ''), "Invalid BMI. BMI must be between 12 and 45. Leave blank if unknown."),
      need((as.numeric(input$crp) >= 0 && as.numeric(input$crp) <= 200) || (input$crp == ''), "Invalid c-reactive protein (CRP) lab result value. CRP must be between 0 and 200 mg/L. Leave blank if unknown."),
      need((as.numeric(input$cdai) >= 0 && as.numeric(input$cdai) <= 600) || (input$cdai == ''), "Invalid Crohn's Disease Activity Index (CDAI). CDAI must be between 0 and 600. Leave blank if unknown."),
    )
    data_raw    <- getInput( input )
    data_clean  <- transformData( data_raw )
    data_format <- formatDrugPref( data_clean )
    data_format
  })
  
  title <- eventReactive(input$submit, {"Treatment Recommendation Results:"})
  
  # Output
  output$rtitle  <- renderText( title() )
  output$results <- renderText( GenerateScripts( data() ) )
  output$plot    <- renderPlot( GenerateResultPlot( data() ) )
}

################################################################################
# LAUNCH APP
################################################################################

shinyApp(ui, server)

################################################################################