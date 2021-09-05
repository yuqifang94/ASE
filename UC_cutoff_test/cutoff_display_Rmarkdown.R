library(shiny)

# Define UI for slider demo app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Sliders"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar to demonstrate various slider options ----
    sidebarPanel(
      
      # Input: Simple integer interval ----
      sliderInput("dNME_cutoff", "dNME cutoff",
                  min = 0, max = 1,
                  value = 0.2,step = 0.001),
      
      # Input: Decimal interval with step value ----
      sliderInput("dMML_cutoff", "dMML cutoff:",
                  min = 0, max = 1,
                  value = 0.5, step = 0.001),
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Table summarizing the values entered ----
      plotOutput("UC_overlap")
      
    )
  )
)

# Define server logic for slider examples ----
server <- function(input, output) {
  
  # Reactive expression to create data frame of all input values ----
  sliderValues <- reactive({
    
    data.frame(
      Name = c("dNME_cutoff",
               "dMML_cutoff"
               ),
      Value = as.character(c(input$dNME_cutoff,
                             input$dMML_cutoff)),
      stringsAsFactors = FALSE)
    
  })
  #Processing data
  # Show the values in an HTML table ----
  output$UC_overlap <- renderPlot({
    plot(x = rnorm(input$n), y = rnorm(input$n))
  })
  
}

# Create Shiny app ----
shinyApp(ui, server)