# Define server logic for slider examples ----
server <- function(input, output,session) {
  
  # Reactive expression to create data frame of all input values ----
  dt_in_mt <- reactive({
    
    dt_in=cutoff_dt[[input$ts]][dNME_cutoff == input$dNME_cutoff&dMML_cutoff==input$dMML_cutoff]
    melt.data.table(dt_in[,list(dNME= shared_dNME/(shared_dNME+dNME_specific),
                                         dMML=shared_dMML/(shared_dMML+dMML_specific),
                                         `dMML random` = shared_dMML_rand/(shared_dMML_rand+dMML_specific_rand),
                                         `dNME random`=shared_dNME_rand/(shared_dNME_rand+dNME_specific_rand),
                                         UC=UC)],id.vars='UC',variable.name='data_type')
    
  })
  #Update slide bar
  observe({
    updateSliderInput(session, "dNME_cutoff", 
                    min = min(cutoff_dt[[input$ts]]$dNME_cutoff), max =  max(cutoff_dt[[input$ts]]$dNME_cutoff))
    updateSliderInput(session, "dMML_cutoff", 
                      min = min(cutoff_dt[[input$ts]]$dMML_cutoff), max =  max(cutoff_dt[[input$ts]]$dMML_cutoff))
  })
  #Processing data

  # Show the values in an HTML table ----
  output$UC_overlap <- renderPlot({
    ggplot(dt_in_mt(),aes(x=UC,y=value))+
      geom_point(aes(color=data_type),size=2)+
      xlab("UC cutoff")+ylab("Porportion of overlapped regions")+ylim(c(0,max(dt_in_mt()$value)*1.2))+
      guides(color=guide_legend(title="",override.aes = list(size=10)))+
      scale_x_reverse()+geom_vline(xintercept = 0.1)+
      theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                            axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                            axis.title.y=element_text(hjust=0.5,size=18,face="bold"), 
                            axis.text.x=element_text(size=16),
                            axis.text.y=element_text(size=16),
                            legend.text=element_text(size=16))
  })
}
  