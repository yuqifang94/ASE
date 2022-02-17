# Define server logic for slider examples ----
server <- function(input, output,session) {
  # Reactive expression to create data frame of all input values ----
  dt_in_mt <- reactive({
    if(input$fix_UC=="dNME_dMML"){
      #fixing dNME and dMML cutoff
      dt_in=cutoff_dt[[input$ts]][dNME_cutoff == as.character(input$dNME_cutoff)&dMML_cutoff==as.character(input$dMML_cutoff)]
      melt.data.table(dt_in[,list(dNME= shared_dNME/(shared_dNME+dNME_specific),
                                           dMML=shared_dMML/(shared_dMML+dMML_specific),
                                           `dMML random` = dMML_rand_shared_dMML /(dMML_rand_shared_dMML +dMML_rand_dMML_specific ),
                                           `dNME random`=dNME_rand_shared_dNME /(dNME_rand_shared_dNME +dNME_rand_dNME_specific ),
                                           UC=UC,UC_region=UC_specific_dNME+shared_dNME,dNME_region=shared_dNME+dNME_specific,dMML_region=shared_dMML+dMML_specific)],
                      id.vars=c('UC','UC_region',
                                'dNME_region','dMML_region'),variable.name='data_type')
    }
    #fixing UC cutoff
    else if(input$fix_UC=="UC"){
      dt_in=cutoff_dt[[input$ts]][UC_cutoff == as.character(input$UC_cutoff)]
    
      dNME_dMML_mt=melt.data.table(dt_in[,list(dNME= shared_dNME/(shared_dNME+dNME_specific),
                                  dMML=shared_dMML/(shared_dMML+dMML_specific),
                                  `dMML random` = UC_rand_shared_dMML/(UC_rand_shared_dMML +UC_rand_dMML_specific),
                                  `dNME random`=UC_rand_shared_dNME/(UC_rand_shared_dNME +UC_rand_dNME_specific ),
                                  dMML_cutoff=dMML_cutoff,dNME_cutoff=dNME_cutoff,
                                  dNMEQuant=dNMEQuant,dMMLQuant=dNMEQuant,
                                  dNME_region=shared_dNME+dNME_specific,dMML_region=shared_dMML+dMML_specific)],
                      id.vars=c('dNME_region','dMML_region','dMMLQuant','dNMEQuant',
                                "dMML_cutoff","dNME_cutoff"),variable.name='data_type')
    
  
      #There're some randomness for the same cutoff for random
      dNME_mt=dNME_dMML_mt[data_type%in%c("dNME","dNME random"),list(value=mean(value),total_region=mean(dNME_region)),
                                      by=list(dNME_cutoff,data_type)]
      colnames(dNME_mt)[1]="cutoff"
      dMML_mt=dNME_dMML_mt[data_type%in%c("dMML","dMML random"),list(value=mean(value),total_region=mean(dMML_region)),
                                      by=list(dMML_cutoff,data_type)]
      colnames(dMML_mt)[1]="cutoff"
      rbind(dNME_mt[!is.na(value)],dMML_mt[!is.na(value)])
    }
                       
    
  })

  #Processing data

  # Show the values in an HTML table ----
 
    output$UC_overlap <- renderPlotly({
      if(input$fix_UC=="dNME_dMML"){
        print(dt_in_mt())
      ggplotly(ggplot(dt_in_mt(),aes(x=as.numeric(as.character(UC)),y=value))+
        geom_point(aes(color=data_type),size=1)+
        xlab("UC cutoff")+ylab("Porportion of overlapped regions")+
        guides(color=guide_legend(title="",override.aes = list(size=10)))+
        scale_x_reverse()+geom_vline(xintercept = 0.1)+
        theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                              axis.title.x=element_text(hjust=0.5,size=14,face="bold"),
                              axis.title.y=element_text(hjust=0.5,size=14,face="bold"),
                              axis.text.x=element_text(size=12),
                              axis.text.y=element_text(size=12),
                              legend.text=element_text(size=12)))
      }else if(input$fix_UC=="UC"){
        ggplotly(ggplot(dt_in_mt(),aes(x=as.numeric(as.character(cutoff)),y=value))+
                   geom_point(aes(color=data_type),size=1)+
                   xlab("dMML or dNME cutoff")+ylab("Porportion of overlapped regions")+
                   guides(color=guide_legend(title="",override.aes = list(size=10)))+
                   scale_x_reverse()+
                   theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                                         axis.title.x=element_text(hjust=0.5,size=14,face="bold"),
                                         axis.title.y=element_text(hjust=0.5,size=14,face="bold"),
                                         axis.text.x=element_text(size=12),
                                         axis.text.y=element_text(size=12),
                                         legend.text=element_text(size=12)))
        
      }
    })
    output$UC_overlap_region <- renderPlotly({
      if(input$fix_UC=="dNME_dMML"){
        ggplotly(ggplot(dt_in_mt(),aes(x=UC_region,y=value))+
          geom_point(aes(color=data_type),size=1)+
          xlab("Number of selected regions")+ylab("Porportion of overlapped regions")+
          guides(color=guide_legend(title="",override.aes = list(size=10)))+
          geom_vline(xintercept = dt_in_mt()[UC==0.1]$UC_region)+
          theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                                axis.title.x=element_text(hjust=0.5,size=14,face="bold"),
                                axis.title.y=element_text(hjust=0.5,size=14,face="bold"),
                                axis.text.x=element_text(size=12),
                                axis.text.y=element_text(size=12),
                                legend.text=element_text(size=12)))
      }else if(input$fix_UC=="UC"){
        ggplotly(ggplot(dt_in_mt(),aes(x=total_region,y=value))+
                   geom_point(aes(color=data_type),size=1)+
                   xlab("Number of selected regions")+ylab("Porportion of overlapped regions")+
                   guides(color=guide_legend(title="",override.aes = list(size=10)))+
                   scale_x_continuous(trans=reverselog_trans(10))+
                   theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                                         axis.title.x=element_text(hjust=0.5,size=14,face="bold"),
                                         axis.title.y=element_text(hjust=0.5,size=14,face="bold"),
                                         axis.text.x=element_text(size=12),
                                         axis.text.y=element_text(size=12),
                                         legend.text=element_text(size=12)))
      }
    })
    observe({
      if(input$fix_UC=="dNME_dMML"){
        updateSliderInput(session, "dNME_cutoff",
                          min = min(as.numeric(as.character(cutoff_dt[[input$ts]]$dNME_cutoff))), 
                          max =  max(as.numeric(as.character(cutoff_dt[[input$ts]]$dNME_cutoff))))
        updateSliderInput(session, "dMML_cutoff",
                          min = min(as.numeric(as.character(cutoff_dt[[input$ts]]$dMML_cutoff))), 
                          max =  max(as.numeric(as.character(cutoff_dt[[input$ts]]$dMML_cutoff))))
      }else if(input$fix_UC=="UC"){
        updateSliderInput(session, "UC_cutoff",
                          min = min(as.numeric(as.character(cutoff_dt[[input$ts]]$UC_cutoff))), 
                          max =  max(as.numeric(as.character(cutoff_dt[[input$ts]]$UC_cutoff))))
      }
        
    })
  #Hide and show based on toggle value
    observeEvent(eventExpr = input$fix_UC, handlerExpr = {
      
      
      if(input$fix_UC=="dNME_dMML"){ 
          shinyjs::hide(id = "UC_cutoff") 
          shinyjs::show(id="UC_overlap")
          shinyjs::show(id = "dNME_cutoff") 
          shinyjs::show(id = "dMML_cutoff") 
          shinyjs::show(id="UC_overlap_region")
    
      }else if(input$fix_UC=="UC"){
        shinyjs::hide(id = "dNME_cutoff") 
        shinyjs::hide(id = "dMML_cutoff") 
        shinyjs::show(id = "UC_cutoff") 
        shinyjs::show(id="UC_overlap_region")
        shinyjs::show(id="UC_overlap")
        
      }
      
      
      })

 
  
  output$debug <-renderText({dt_in_mt()$value})
  #Update slide bar


}
  