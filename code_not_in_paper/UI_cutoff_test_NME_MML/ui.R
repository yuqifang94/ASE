library(shiny)

# Define UI for slider demo app ----
ui <- fluidPage(
  useShinyjs(),
  # App title ----
  titlePanel(""),
  # Main panel for displaying outputs ----
  
  # Output: Table summarizing the values entered ----
  plotlyOutput("UC_overlap"),
  plotlyOutput("UC_overlap_region"),

  #textOutput("debug"),
  hr(),
  # Sidebar layout with input and output definitions ----
  radioButtons(inputId = "fix_UC", label=NULL,choices  = list("Fixing UC"="UC","Fixing dNME and dMML"="dNME_dMML"), selected  = "dNME_dMML"),
  fluidRow(
    column(12,
     
      # Input: Simple integer interval ----
      selectInput("ts", "Tissue:", 
                  choices=names(cutoff_dt)),
      sliderInput("dNME_cutoff", "dNME cutoff",
                  min = 0, max = 1,
                  value = 0.55,step = 0.05),
      
      # Input: Decimal interval with step value ----
      sliderInput("dMML_cutoff", "dMML cutoff:",
                  min = 0, max = 1,
                  value = 0.5, step = 0.05),
      sliderInput("UC_cutoff", "UC cutoff",
                  min = 0, max = 1,
                  value = 0.1,step = 0.01)
  
    )
    
 

      
  )
)




