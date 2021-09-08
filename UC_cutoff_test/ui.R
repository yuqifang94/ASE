library(shiny)

# Define UI for slider demo app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel(""),
  # Main panel for displaying outputs ----
  
  # Output: Table summarizing the values entered ----
  plotOutput("UC_overlap"),
  plotOutput("UC_overlap_region"),
  textOutput("debug"),
  hr(),
  # Sidebar layout with input and output definitions ----
  fluidRow(
    column(12,
      
      # Input: Simple integer interval ----
      selectInput("ts", "Tissue:", 
                  choices=names(cutoff_dt)),
      sliderInput("dNME_cutoff", "dNME cutoff",
                  min = 0, max = 1,
                  value = 0.5,step = 0.05),
      
      # Input: Decimal interval with step value ----
      sliderInput("dMML_cutoff", "dMML cutoff:",
                  min = 0, max = 1,
                  value = 0.45, step = 0.05)
  
    )
 

      
  )
)




