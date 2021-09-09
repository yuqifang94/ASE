# Define server logic for slider examples ----
server <- function(input, output,session) {
  
  # Reactive expression to create data frame of all input values ----
  dt_in_mt <- reactive({
    #fixing dNME and dMML cutoff
    dt_in=cutoff_dt[[input$ts]][dNME_cutoff == as.character(input$dNME_cutoff)&dMML_cutoff==as.character(input$dMML_cutoff)]
    melt.data.table(dt_in[,list(dNME= shared_dNME/(shared_dNME+dNME_specific),
                                         dMML=shared_dMML/(shared_dMML+dMML_specific),
                                         `dMML random` = shared_dMML_rand/(shared_dMML_rand+dMML_specific_rand),
                                         `dNME random`=shared_dNME_rand/(shared_dNME_rand+dNME_specific_rand),
                                         UC=UC,UC_region=UC_specific_dNME+shared_dNME,dNME_region=shared_dNME+dNME_specific,dMML_region=shared_dMML+dMML_specific)],
                    id.vars=c('UC','UC_region',
                              'dNME_region','dMML_region'),variable.name='data_type')
    #fixing UC cutoff
    # dt_in=cutoff_dt[[input$ts]][UC_cutoff == as.character(input$UC_cutoff)]
    # dNME_dMML_mt=melt.data.table(dt_in[,list(dNME= shared_dNME/(shared_dNME+dNME_specific),
    #                             dMML=shared_dMML/(shared_dMML+dMML_specific),
    #                             `dMML random` = shared_dMML_rand/(shared_dMML_rand+dMML_specific_rand),
    #                             `dNME random`=shared_dNME_rand/(shared_dNME_rand+dNME_specific_rand),
    #                             dNME_region=shared_dNME+dNME_specific,dMML_region=shared_dMML+dMML_specific)],
    #                 id.vars=c('dNME_region','dMML_region'),variable.name='data_type')
    #melt.data.table(dNME_dMML_mt,id.vars=c('data_type','value'),value.name='total_region')

    
  })

  #Processing data

  # Show the values in an HTML table ----
  output$UC_overlap <- renderPlotly({
    ggplotly(ggplot(dt_in_mt(),aes(x=UC,y=value))+
      geom_point(aes(color=data_type),size=2)+
      xlab("UC cutoff")+ylab("Porportion of overlapped regions")+
      guides(color=guide_legend(title="",override.aes = list(size=10)))+
      scale_x_reverse()+geom_vline(xintercept = 0.1)+
      theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                            axis.title.x=element_text(hjust=0.5,size=14,face="bold"),
                            axis.title.y=element_text(hjust=0.5,size=14,face="bold"),
                            axis.text.x=element_text(size=12),
                            axis.text.y=element_text(size=12),
                            legend.text=element_text(size=12)))
  })
  output$UC_overlap_region <- renderPlotly({
    ggplotly(ggplot(dt_in_mt(),aes(x=UC_region,y=value))+
      geom_point(aes(color=data_type),size=2)+
      xlab("Number of regions")+ylab("Porportion of overlapped regions")+
      guides(color=guide_legend(title="",override.aes = list(size=10)))+
      geom_vline(xintercept = dt_in_mt()[UC==0.1]$UC_region)+
      theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                            axis.title.x=element_text(hjust=0.5,size=14,face="bold"),
                            axis.title.y=element_text(hjust=0.5,size=14,face="bold"),
                            axis.text.x=element_text(size=12),
                            axis.text.y=element_text(size=12),
                            legend.text=element_text(size=12)))
  })
  # output$UC_overlap_region_fix_UC <- renderPlotly({
  #   ggplotly(ggplot(dt_in_mt(),aes(x=UC_region,y=value))+
  #              geom_point(aes(color=data_type),size=2)+
  #              xlab("Number of regions")+ylab("Porportion of overlapped regions")+
  #              guides(color=guide_legend(title="",override.aes = list(size=10)))+
  #              geom_vline(xintercept = dt_in_mt()[UC==0.1]$UC_region)+
  #              theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
  #                                    axis.title.x=element_text(hjust=0.5,size=14,face="bold"),
  #                                    axis.title.y=element_text(hjust=0.5,size=14,face="bold"),
  #                                    axis.text.x=element_text(size=12),
  #                                    axis.text.y=element_text(size=12),
  #                                    legend.text=element_text(size=12)))
  # })
  # output$debug <-renderText({dt_in_mt()$value})
  #Update slide bar
  observe({
    updateSliderInput(session, "dNME_cutoff", 
                      min = min(cutoff_dt[[input$ts]]$dNME_cutoff), max =  max(cutoff_dt[[input$ts]]$dNME_cutoff))
    updateSliderInput(session, "dMML_cutoff", 
                      min = min(cutoff_dt[[input$ts]]$dMML_cutoff), max =  max(cutoff_dt[[input$ts]]$dMML_cutoff))
  })
}
  