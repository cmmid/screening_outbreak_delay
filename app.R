
suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(tidyverse)
  library(dplyr)
  library(furrr)
  #install.packages("emojifont")
  library(emojifont)
  #library(cowplot)
  library(gridExtra)
  library(knitr)
  library(kableExtra)
  library(shinyWidgets)
  library(ggrepel)
  library(scales)
  
})

source("utils.R")


ui <- shinyUI(fluidPage(
  titlePanel("Effectiveness of airport-based interventions at detecting travellers and delaying an outbreak of COVID-19 (formerly 2019-nCoV)"),
  shiny::includeMarkdown("date_stamp.md"),
  sidebarLayout(
    sidebarPanel(
      #tags$label(class="h3",)
      sliderInput('dur.flight', 
                  label = "Travel duration (hours)",
                  min = 1, max = 20, step = 1, value = 12),
      sliderInput(inputId = 'sens.exit',
                  label = 'Sensitivity of exit screening',
                  value = 86, min = 0, max = 100, step = 1, post  = " %"),
      sliderInput(inputId = 'sens.entry',
                  label = 'Sensitivity of entry screening',
                  value = 86, min = 0, max = 100, step = 1, post  = " %"),
      hr(),
      
      selectInput('pathogen',label='Pathogen',
                  choices = unique(pathogen$name),
                  selected = pathogen$name[1]),
      
      conditionalPanel(
        condition = "input.conditionedPanels=='Delaying an outbreak'",
        sliderTextInput('n_infected', 
                        label = "Number of infected travellers per week",
                        choices = c(1,10,100), selected=1, grid=T), 
        sliderInput('sensitisation',
                    label = "Effectiveness of traveller sensitisation",
                    min = 0, max = 70, step = 10, value = 50, post= " %"), 
        sliderInput(inputId = "r0",
                    label = "R0 (lower and upper limits of 95% interval)",
                    min=0.1,max=5,step=0.1,value=c(1.4,3.9)),               
        sliderInput("k",
                    label="Dispersion (k)",
                    min=0.1,max=3.0,step=0.01,value = 0.54),
        checkboxInput(inputId = "logScale", 
                      label = "Logarithmic x scale",
                      value = TRUE)),
      conditionalPanel(
        condition = "input.conditionedPanels %in% c('Delaying an outbreak')",
        numericInput("mu_inc",
                     'Days from infection to symptom onset (mean)',
                     value = 5.2, min = 0.1, max = 20, step = 0.1),
        numericInput("sigma_inc",
                     "Days from infection to symptom onset (variance)",
                     value = 4.1, min = 0.1, max = 20, step = 0.1),
        numericInput("mu_inf",
                     'Days from symptom onset to severe symptoms e.g hospitalisation (mean)',
                     value = 9.1, min = 0.1, max = 20, step = 0.1),
        numericInput("sigma_inf",
                     "Days from symptom onset to severe symptoms e.g hospitalisation (variance)",
                     value = 14.7, min = 0.1, max = 20, step = 0.1),
        sliderInput(inputId = 'prop.asy',
                    label = 'Proportion of cases that are asymptomatic',
                    value = 17, min = 0, max = 100, step = 1)
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(title = "Delaying an outbreak",
                 shiny::includeMarkdown("delay_description.md"),
                 fluidRow(
                   uiOutput("delay_plot"))),
        tabPanel(title = "Model",
                 fluidRow(shiny::includeMarkdown("assumptions.md"))),
        tabPanel(title = "References",
                 shiny::includeMarkdown("references.md")),
        id="conditionedPanels"
      )
    )
  )
)
)



server <- function(input, output, session){
  
  
  
  observe({
    default=input$pathogen
    
    updateNumericInput(session,"prop.asy",
                       value = pathogen %>% 
                         filter(name==default) %>% pull(prop.asy))
    updateNumericInput(session,"mu_inc",
                       value = pathogen %>% 
                         filter(name==default) %>% pull(mu_inc))
    updateNumericInput(session,"sigma_inc",
                       value = pathogen %>% 
                         filter(name==default) %>% pull(sigma_inc))
    updateNumericInput(session,"mu_inf",
                       value = pathogen %>% 
                         filter(name==default) %>% pull(mu_inf))
    updateNumericInput(session,"sigma_inf",
                       value = pathogen %>% 
                         filter(name==default) %>% pull(sigma_inf))
    
    
    pathv <- pathogen %>% filter(name == default)
    
    
    updateSliderInput(session, inputId = "r0",
                      value = as.numeric(round(pathv[,c("r0l","r0u")],1)))
    
    updateSliderInput(session,inputId = "k", value = pathv$k)
    
  }
  
  )
  
  
  
  nat_hist_periods <- reactive({
    
    
    #convert to gamma parameters
    periods <- data.frame(
      inc_period = time_to_event(1e4, input$mu_inc, input$sigma_inc),
      inf_period = time_to_event(1e4, input$mu_inf, input$sigma_inf))
    return(periods)
  })
  
  delays <- reactive({
    
    #browser()
    
    pathv <- filter(pathogen, name == input$pathogen)
    
    histories <- generate_histories(input) %>%
      filter(inf + incu >= flight.departure)
    
    outcomes  <- calc_outcomes(histories, input)
    
    theta <- list(`No screening` = 1,
                  `Exit only` = 0.46,
                  `Exit and entry` = 0.42) %>% 
      map_df(~data.frame(frac_undetected = .x),.id = "Screening") %>%
      mutate(Screening = fct_inorder(Screening))
    
    n_R0 <- 100
    
    
    
    if (input$r0[1] == input$r0[2]){
      R0 <- rep(input$r0[1], n_R0)
    } else {
      R0_parms <- gamma.parms.from.quantiles(q = input$r0,p = c(0.025, 0.975)) 
      R0 <- rgamma(n     = n_R0,
                   shape = R0_parms$shape,
                   rate  = R0_parms$rate)
    }
    
    
    
    delay_parameter_grid_new <- 
      
      ## PARAMETERS
      data.frame(R0 = R0,
                 rowid = 1:n_R0) %>%
      #crossing(k_labels) %>%
      arrange(rowid) %>%
      mutate(T0 = T0(R0 = R0, k = input$k, c = 0.5),
             T0 = ifelse(T0 <= 0, Inf, T0)) %>%
      
      ## INTERVENTIONS AND ARRIVAL RATES
      crossing(lambda           = input$n_infected,
               theta,
               contact_tracing  = input$sensitisation/100,
               k = input$k,
      ) %>%
      
      ## DELAYS
      mutate(Rq      = (1 - contact_tracing)*R0,
             q0      = runif(n()),
             t0      = qgamma(q0, shape=T0, rate=lambda/7),
             Tq      = T0(Rq, k, 0.5),
             Tq      = ifelse(Tq <= 0, Inf, Tq),
             tdash   = qgamma(q0, shape=Tq, rate=frac_undetected*lambda/7),
             delta_t = tdash - t0,
             delta_t = ifelse(is.infinite(t0), Inf, delta_t)
      ) 
    
    
    
    
    delays <- delay_parameter_grid_new %>%
      filter(!(contact_tracing == 0 & Screening == "No screening"))
    
    return(delays)
  })
  
  
  
  height = function() {
    session$clientData$output_delayplot_width
  }
  
  
  output$delayplot <- renderPlot(expr={
    
    delays <- delays()
    
    #add .54 estimate for COVID-19
    # use factor to order with COVID-19 first
    
    delay_summaries <- delays %>%  
      filter(contact_tracing==(input$sensitisation/100)) %>% 
      filter(lambda==(input$n_infected)) %>% 
      filter(k==input$k) %>% 
      #filter(Screening=="Exit and entry") %>% 
      group_by(Screening, contact_tracing, lambda, k) %>%
      summarise(median = median(delta_t),
                lower_025  = quantile(delta_t, 0.025,na.rm=T),
                lower_250  = quantile(delta_t, 0.25,na.rm=T),
                median     = median(delta_t, na.rm = T),
                upper_750  = quantile(delta_t, 0.75,na.rm=T),
                upper_975  = quantile(delta_t, 0.975,na.rm=T),
                averted    = mean(is.infinite(delta_t))) 
    
    
    delays_breaks <- c(0.5, 10^seq(0,3, length.out = 100))
    
    
    nesting_vars <- c("lambda",
                      "Screening",
                      "contact_tracing")
    
    delta_t_out   <- 10^seq(-1, 3, length.out = 101)
    
    
    delays_to_plot <- delays %>% 
      nest(data =  -all_of(nesting_vars))
    
    
    delays_to_plot_processed <- 
      delays_to_plot %>%
      # filter out when there's no intervention
      filter(!(contact_tracing == 0 & Screening == "No screening")) %>%
      # split and unnest so we don't lose grouping variables
      mutate(data = map(data, ~mutate(.x, 
                                      delta_t_ecdf = ecdf(delta_t)(delta_t)))) %>%
      mutate(data_interp = map(
        data,
        make_interp,
        delta_t_out = delta_t_out)) %>%
      unnest(data_interp) %>%
      # allocate to quantile ranges in case we want to fill under line
      mutate(prob = cut(delta_t_ecdf,
                        breaks = c(0,probs,1),
                        include.lowest = T))
    
    
    palette_interventions <- c("No screening"   = "grey60",
                               "Exit only"      = "#41b6c4",
                               "Exit and entry" = "#225ea8")
    
    #    browser()
    
    
    delays_plot <- 
      delays_to_plot_processed %>% 
      ggplot() +
      geom_line(aes(x=delta_t,#ymin=0,
                    y=1-delta_t_ecdf,
                    color = Screening
                    #fill=prob
      ),
      size = 1,
      show.legend = FALSE)+
      
      scale_y_continuous(labels=percent,
                         limits = c(0,1))+
      geom_segment(data = delays_to_plot_processed %>%
                     filter(delta_t == min(delta_t)),
                   aes(xend  = delta_t,
                       color = Screening,
                       y     = 1 - delta_t_ecdf,
                       yend  = 1-  delta_t_ecdf),
                   x    = -Inf,
                   size = 1)  +
      geom_segment(data = delays_to_plot_processed %>%
                     group_by(Screening) %>%
                     filter(is.finite(delta_t),
                            !is.na(delta_t_ecdf),
                            delta_t < 1e3) %>%
                     filter(delta_t == max(delta_t, na.rm=T)),
                   aes(x  = delta_t,
                       color = Screening,
                       y     = 1 - delta_t_ecdf,
                       yend  = 1 - delta_t_ecdf),
                   xend    = 999,
                   size = 1)  +
      scale_color_manual(values = palette_interventions) +
      scale_fill_manual(values = palette_interventions) +
      
      geom_label_repel(
        data = delay_summaries,
        aes(
          x = median,
          y = 0.5,
          fill = Screening,
          
          label = ifelse(
            is.infinite(median),
            "50%: no outbreak",
            paste0("50%: ≥ ", round(median), " days delay")
          )
        ),
        nudge_x = 0.5,
        nudge_y = 0.1,
        segment.colour = "black",
        color = "white"
      ) +
      geom_label_repel(
        data = delay_summaries,
        aes(
          x = lower_250,
          y = 0.75,
          fill = Screening,
          
          label = ifelse(
            is.infinite(lower_250),
            "75%: no outbreak",
            paste0("75%: ≥ ", round(lower_250), " days delay")
          )
        ),
        nudge_x = 0.5,
        nudge_y = 0.1,
        hjust = 0,
        segment.colour = "black",
        color = "white"
      ) +
      geom_label_repel(
        data = delay_summaries,
        aes(
          x = upper_750,
          y = 0.25,
          fill = Screening,
          
          label = ifelse(
            is.infinite(upper_750),
            "25%: no outbreak",
            paste0("25%: ≥ ", round(upper_750), " days delay")
          )
        ),
        nudge_x = 0.5,
        nudge_y = 0.1,
        hjust = 0,
        segment.colour = "black",
        color = "white"
      ) +
      geom_label_repel(
        data = delay_summaries,
        aes(
          x = lower_025,
          y = 0.975,
          fill = Screening,
          
          label = ifelse(
            is.infinite(lower_025),
            "97.5%: no outbreak",
            paste0("97.5%: ≥ ", round(lower_025), " days delay")
          )
        ),
        nudge_x = -0.5,
        nudge_y = -0.2,
        hjust = 1,
        segment.colour = "black",
        color = "white"
      )+
      geom_label_repel(
        data = delay_summaries,
        aes(
          x = upper_975,
          y = 0.025,
          fill = Screening,
          
          label = ifelse(
            is.infinite(upper_975),
            "2.5%: no outbreak",
            paste0("2.5%: ≥ ", round(upper_975), " days delay")
          )
        ),
        nudge_x = 0.5,
        nudge_y = 0.1,
        hjust = 1,
        segment.colour = "black",
        color = "white"
      )+
      
      xlab(
        "Delay to outbreak (days)")+
      ylab("Percentage of simulations\n with delays at least this long")+
      theme_bw()+
      theme(panel.grid.minor.y=element_blank(),
            legend.position = "none",
            text=element_text(size=18)) +
      guides(fill = FALSE) +
      facet_grid(Screening ~ .)
    
    if (input$logScale){
      delays_plot <- 
        delays_plot + 
        scale_x_log10(breaks=log_breaks(10),
                      #labels = scales::label_number(),
                      expand = c(0,0), 
                      labels=label_number_si(),
                      limits=c(1,1000))
    } else{
      delays_plot <- 
        delays_plot + 
        scale_x_continuous(labels=label_number_si(),
                           expand = c(0,0), 
                           limits = c(0,1000)) 
    }
    
    delays_plot
    
  }
  
  )
  
  output$delay_plot <- renderUI({
    plotOutput("delayplot")
  })
  
}

shinyApp(ui=ui,server=server)
