
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#    http://shiny.rstudio.com/

library(tidyverse)
library(ggplot2)
library(usdata)
library(ggthemes)
library(maps)
library(ggtext)
library(ggnewscale)
library(ggpubr)
library(plotly)
library(gridExtra)

library(haven)
library(lmw)
library(data.table) ## For some minor data wrangling
library(fixest)     ## NB: Requires version >=0.9.0
library(cdlTools)
library(dplyr)
library(grDevices)
library(fastDummies)
library(lfe)
library(sbw)
library(parallel)

library(shiny)
library(sf)
library(leaflet)
library(shinyTree)
source("util.R")

#######################################################################
##    Global 
#######################################################################

obs_colors <- c(Treat = "#9dc4e4", 
                Control = "#c6e0b4", 
                Invalid = "white")

# ess_obs_colors <- c(IdealTreat = "#3293e3", 
#                     IdealControl = "#69b334", 
#                     InvTreat = "#c5dbed", 
#                     InvControl = "#ddedd1", 
#                     LimitAnticipate = "#ffc001", 
#                     DelayOnset = "#D2691E", 
#                     Dissipate = "#f56c5b",
#                     InvalidTreat = "#e1e1e1", 
#                     InvalidControl = "#fafafa")

ess_obs_colors <- c(IdealExperiment= "#3293e3", 
                    TimeInvariance = "#69b334", 
                    LimitedAnticipation = "#ffc001", 
                    DelayedOnset = "#D2691E", 
                    EffectDissipation = "#f56c5b",
                    Invalid = "#e1e1e1")

panelview_plot_tyle <- 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="white", size=0.5, linetype="solid"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title=element_text(size=10),
        axis.title.x = element_text(size = 12, margin = margin(t = 8, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 8, b = 0, l = 0)),
        axis.text = element_text(color="black", face="bold", size=12),
        axis.text.x = element_text(face="bold", size = 10, angle = 0, hjust=0.5, vjust=0),
        axis.text.y = element_text(face="bold", size = 10),
        plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.position = "bottom",
        legend.margin = margin(c(0, 5, 5, 0)),
        legend.text = element_text(margin = margin(r = 10, unit = "pt"), size = 10),
        plot.title = element_text(size=15, hjust = 0.5, face="bold", margin = margin(8, 0, 8, 0)))
guides(fill = guide_legend(order = 1), 
       shape = guide_legend(order = 2))

implied_plot_tyle <- 
  theme(plot.title = element_text(family = "Times", size = rel(2), hjust = 0.5),
        text = element_text(family = "Times"),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = "white"),
        axis.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.margin = unit(0, "cm"),
        plot.margin=unit(c(1,0,0,0),unit = "mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold", family = "Times"))

#######################################################################
##    UI
#######################################################################

ui <- fluidPage(
  
  theme = shinythemes::shinytheme("sandstone"), 
  
  # Application title
  titlePanel("Visualize Event Study Dataset"),
  
  navbarPage("Visualizations",
             
             tabPanel("Data",
                      
                      sidebarLayout(
                        
                        sidebarPanel(
                          
                          h4("Data"), 
                          fileInput(inputId = "file", 
                                    label = "Upload data from your event study:",
                                    accept = ".csv"), 
                          tags$hr(), 
                          
                          h4("Variable"), 
                          textInput(inputId = "outcome_name", 
                                    label = "What's the variable name of your outcome:", 
                                    value = "asmrs"), 
                          textInput(inputId = "treat_name", 
                                    label = "What's the variable name of your binary treatment:", 
                                    value = "treat_status"),
                          textInput(inputId = "time_ind_name", 
                                    label = "What's the variable name of your time indicator:", 
                                    value = "year"),
                          textInput(inputId = "unit_ind_name", 
                                    label = "What's the variable name of your unit indicator:", 
                                    value = "state"), 
                          textInput(inputId = "covar_name", 
                                    label = "What's the name(s) of your control covariates (seperated by comma, e.g. 'var1,var2,var3'): ", 
                                    value = "pcinc,asmrh,cases")
                        ), 
                        
                        mainPanel(
                          navbarPage(
                            title = "", 
                            tabPanel("Visualize the data", 
                                     plotlyOutput("panel_fig", height = "600px", width = "900px")), 
                            tabPanel("Proportion of units under treatment each year", 
                                     plotOutput("treat_prop_fig", height = "600px", width = "900px"))
                            
                          )
                        )
                      )
             ), 
             
             tabPanel("Estimand and Assumptions",
                      
                      sidebarLayout(
                        
                        sidebarPanel(
                          
                          h4("Estimand"), 
                          radioButtons("estimand1", "Choose the estimand:",
                                       c("ATE\\(^{P}_{t_y} (t_1, t_p)\\)" = "new",
                                         "ATE\\(^{P}_{t_y} (t_1, \\infty)\\)" = "std")),
                          withMathJax(),
                          helpText('Estimand can be ATE\\(^{P}_{t_y} (t_1, t_p)\\), which is the average causal 
                      effect the target population observed at time \\(t_y\\) from adopting the treatment 
                      for the first time at time \\(t_1\\) to adopting the treatment sometime after time \\(t_1\\).'), 
                          helpText('A special case of the estimand is ATE\\(^{P}_{t_y} (t_1, \\infty)\\), which is is the average causal 
                      effect in the target population observed at time \\(t_y\\) from adopting the treatment for the first 
                     time at time \\(t_1\\) to never adopting the treatment.'), 
                          
                          numericInput("t1_1", "What is the \\(t_1\\) in the target estimand? ", value = 1975), 
                          withMathJax(),
                          helpText('\\(t_1\\) denotes the target treatment initiation time.'),
                          
                          numericInput("ty_1", "What is the \\(t_y\\) in the target estimand? ", value = 1980), 
                          withMathJax(),
                          helpText('\\(t_y\\) denotes the outcome measurement time.'),
                          
                          # Horizontal line 
                          tags$hr(), 
                          h4("Assumptions"), 
                          helpText('Consider scenarios when the following assumptions hold or not hold.'),
                          
                          checkboxInput("invariance", 'Time Shift Invariance. ', FALSE), 
                          checkboxInput("anticipation", 'Limited Treatment Anticipation.', FALSE), 
                          conditionalPanel(
                            condition = "input.estimand1 == 'std'",
                            checkboxInput("delay", 'Delayed Onset of Treatment Effect.', FALSE)
                          ), 
                          checkboxInput("washout", 'Treatment Effect Dissipation.', FALSE), 
                          checkboxInput("highlight", 'Show observations with different colors.', FALSE)
                        ),
                        mainPanel(plotOutput("panelview", height = "600px", width = "900px"))
                      )
             ), 
             
             
             tabPanel("Two-way Fixed Effect Regression",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Estimand"), 
                          radioButtons("estimand6", "Choose the estimand:",
                                       c("ATE\\(^{P}_{t_y} (t_1, \\infty)\\)" = "std")),
                          withMathJax(),
                          helpText('Estimand is ATE\\(^{P}_{t_y} (t_1, \\infty)\\), which is is the average causal 
                      effect in the target population observed at time \\(t_y\\) from adopting the treatment for the first 
                     time at time \\(t_1\\) to never adopting the treatment.'), 
                          
                          numericInput("t1_6", "What is the \\(t_1\\) in the target estimand? ", value = 1975), 
                          withMathJax(),
                          helpText('\\(t_1\\) denotes the target treatment initiation time.'),
                          
                          numericInput("ty_6", "What is the \\(t_y\\) in the target estimand? ", value = 1980), 
                          withMathJax(),
                          helpText('\\(t_y\\) denotes the outcome measurement time.'), 
                          
                          # Horizontal line 
                          tags$hr(), 
                          h4("Regression"), 
                          h5("If to fit the dynamic TWFE model with the following form, please specify the minimum and maximum relative treatment periods \\(a\\) and \\(b\\)."), 
                          h5("\\(Y_{it} = \\alpha_i + \\beta_t + \\sum_{l = a}^{b} \\tau_l \\mathbf{1}_{l = t_i - G_i}\\)."), 
                          h5("\\(a\\) and \\(b\\) in the dynamic specification define the treatment horizon."), 
                          h5("Specifically, we allow anticipation effects due to the treatment \\(-a\\) periods 
                          after the current time point and carryover effects due to the treatment \\(b\\) periods 
                          before the current time point. To specify the minimum relative treatment period \\(a\\), 
                          just input a positive nummber."), 
                          numericInput("l_min", "What is the treatment horizon: minimum lead?", value = 15), 
                          numericInput("l_max", "What is the treatment horizon: maximum lag?", value = 20)
                        ), 
                        mainPanel(
                          plotOutput("twfe_event_study_plot", height = "600px", width = "900px"), 
                          plotOutput("twfe_panelview", height = "600px", width = "900px")
                        )
                      )
             ), 
             
             
             tabPanel("Implied Population",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Estimand"), 
                          radioButtons("estimand2", "Choose the estimand:",
                                       c("ATE\\(^{P}_{t_y} (t_1, \\infty)\\)" = "std")),
                          withMathJax(),
                          helpText('Estimand is ATE\\(^{P}_{t_y} (t_1, \\infty)\\), which is is the average causal 
                      effect in the target population observed at time \\(t_y\\) from adopting the treatment for the first 
                     time at time \\(t_1\\) to never adopting the treatment.'), 
                          
                          numericInput("t1_2", "What is the \\(t_1\\) in the target estimand? ", value = 1975), 
                          withMathJax(),
                          helpText('\\(t_1\\) denotes the target treatment initiation time.'),
                          
                          numericInput("ty_2", "What is the \\(t_y\\) in the target estimand? ", value = 1980), 
                          withMathJax(),
                          helpText('\\(t_y\\) denotes the outcome measurement time.'),
                          
                          # Horizontal line 
                          tags$hr(), 
                          h4("Implied Population of Units"), 
                          checkboxInput("implied_unit", 'Show implied Units.', FALSE), 
                          
                          conditionalPanel(
                            condition = "input.implied_unit",
                            checkboxInput("implied_map", "If the Units are US States, show a US map.", FALSE), 
                            
                            checkboxInput("implied_unit_by_time", 'Show implied Units by Time.', FALSE), 
                            
                            conditionalPanel(
                              condition = "input.implied_map",
                              checkboxInput("shade", 'Show differential implied with shade.', FALSE)
                            ), 
                            
                            conditionalPanel(
                              condition = "input.implied_unit_by_time",
                              textInput("time_value", "Show implied Units by the following Time point:", value = "1980")
                            )
                            
                          ),
                          
                          # Horizontal line 
                          tags$hr(), 
                          h4("Implied Population of Time Periods"), 
                          checkboxInput("implied_time", 'Show implied Time periods.', FALSE), 
                          conditionalPanel(
                            condition = "input.implied_time",
                            checkboxInput("implied_time_by_unit", 'Show implied Time periods by Unit.', FALSE), 
                            conditionalPanel(
                              condition = "input.implied_time_by_unit",
                              textInput("unit_value", "Show implied Time periods by the following Unit:", value = "MA")
                            )
                          ),  
                          
                          # Horizontal line 
                          tags$hr(), 
                          h4("Download Dataset with Implied Weights"), 
                          downloadButton("download_lmw", "Download")
                          
                        ), 
                        mainPanel(
                          conditionalPanel(
                            condition =  "input.implied_unit",
                            plotlyOutput("implied_unit_plot", height = "600px", width = "1000px")
                          ), 
                          conditionalPanel(
                            condition = "input.implied_map",
                            plotOutput("implied_map_plot", height = "600px", width = "1000px")
                          ), 
                          conditionalPanel(
                            condition = "input.implied_time",
                            plotlyOutput("implied_time_plot", height = "600px", width = "1000px")
                          )
                        )
                        
                      )
             ),
             
             tabPanel("Effective Sample Size",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Estimand"), 
                          radioButtons("estimand3", "Choose the estimand:",
                                       c("ATE\\(^{P}_{t_y} (t_1, \\infty)\\)" = "std")),
                          withMathJax(),
                          helpText('Estimand is ATE\\(^{P}_{t_y} (t_1, \\infty)\\), which is is the average causal 
                      effect in the target population observed at time \\(t_y\\) from adopting the treatment for the first 
                     time at time \\(t_1\\) to never adopting the treatment.'), 
                          
                          numericInput("t1_4", "What is the \\(t_1\\) in the target estimand? ", value = 1975), 
                          withMathJax(),
                          helpText('\\(t_1\\) denotes the target treatment initiation time.'),
                          
                          numericInput("ty_4", "What is the \\(t_y\\) in the target estimand? ", value = 1980), 
                          withMathJax(),
                          helpText('\\(t_y\\) denotes the treatment initiation time.'),
                          
                          # Horizontal line 
                          tags$hr(), 
                          checkboxInput("twfe", 'Compute ESS for TWFE approach.', TRUE), 
                          checkboxInput("balance_experiment", 'Compute ESS for balancing Approach under Hypothetical Experiment.', FALSE), 
                          checkboxInput("balance_invariance", 'Compute ESS for balancing Approach if further invoke Invariance to Time Shifts assumption.', FALSE), 
                          checkboxInput("balance_anticipate", 'Compute ESS for balancing Approach if further invoke Limited Treatment Anticipation assumption.', FALSE), 
                          checkboxInput("balance_delay", 'Compute ESS for balancing Approach if further invoke Delayed Treatment Onset assumption.', FALSE), 
                          checkboxInput("balance_dissipate", 'Compute ESS for balancing Approach using all observations.', FALSE)
                          
                        ),
                        mainPanel(
                          plotOutput("ess_plot", height = "600px", width = "1000px"),
                          conditionalPanel(
                            condition = "input.twfe",
                            tableOutput("twfe_ess")
                          ),
                          conditionalPanel(
                            condition = "input.balance_experiment",
                            tableOutput("balance_experiment_ess")
                          ),
                          conditionalPanel(
                            condition = "input.balance_invariance",
                            tableOutput("balance_invariance_ess")
                          ),
                          conditionalPanel(
                            condition = "input.balance_anticipate",
                            tableOutput("balance_anticipate_ess")
                          ), 
                          conditionalPanel(
                            condition = "input.balance_delay",
                            tableOutput("balance_delay_ess")
                          ),  
                          conditionalPanel(
                            condition = "input.balance_dissipate",
                            tableOutput("balance_dissipate_ess")
                          ))
                      )
             ), 
             
             tabPanel("Sample Influence",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Estimand"),
                          radioButtons("estimand5", "Choose the estimand:",
                                       c("ATE\\(^{P}_{t_y} (t_1, \\infty)\\)" = "std")),
                          withMathJax(),
                          helpText('Estimand is ATE\\(^{P}_{t_y} (t_1, \\infty)\\), which is is the average causal 
                      effect in the target population observed at time \\(t_y\\) from adopting the treatment for the first 
                     time at time \\(t_1\\) to never adopting the treatment.'), 
                          
                          numericInput("t1_5", "What is the \\(t_1\\) in the target estimand? ", value = 1975),
                          withMathJax(),
                          helpText('\\(t_1\\) denotes the target treatment initiation time.'),

                          numericInput("ty_5", "What is the \\(t_y\\) in the target estimand? ", value = 1980),
                          withMathJax(),
                          helpText('\\(t_y\\) denotes the outcome measurement time.'),

                          # Horizontal line
                          tags$hr(),
                          selectInput("influence_metric", "What influence metric to show:",
                                      c("Change of point estimate due to each observation." = "est_change",
                                        "Exact SIC of each observation." = "sic",
                                        "Scaled SIC of each observation." = "sic_scaled")),

                          # Horizontal line
                          tags$hr(),
                          h4("Download Table with Influence of Each Observation"),
                          downloadButton("download_inf", "Download")
                        ),
                        mainPanel(
                          plotOutput("infuence_plot_panel", height = "600px", width = "1000px"),
                          plotlyOutput("infuence_plot", height = "600px", width = "1000px")
                        )
                      )
             )
  )
)


#######################################################################
##    Backend Control 
#######################################################################

server = function(input, output){
  
  # reactive data frame of the study data, without all dummy columns; will be used for panel plotting
  input_data <- reactive({
    infile <- input$file
    if (is.null(infile)) {
      return(NULL)
    }
    data = read.csv(infile$datapath, header = TRUE)
    
    outcome = data[, input$outcome_name]
    treat = data[, input$treat_name]
    time = data[, input$time_ind_name]
    unit = data[, input$unit_ind_name]
    idx = 1:nrow(data)
    covar_names = str_split(input$covar_name, ",")[[1]]
    covars_df = data.frame(data[, covar_names])
    col_names = paste("X", 1:ncol(covars_df), sep = "")
    colnames(covars_df) = col_names
    
    df = data.frame(cbind(idx, unit, time, outcome, treat, covars_df))
    colnames(df) = c("X", "Unit", "Time", "Outcome", "Treatment", col_names)
    
    df = df %>% 
      mutate(Time = as.numeric(as.character(Time)), 
             Treatment = factor(Treatment, level = c(1, 0))) %>% 
      group_by(Unit) %>% 
      # no non-missing arguments to min; returning Inf
      mutate(TreatStartTime = min(Time[Treatment == 1])) %>% ungroup() %>% 
      mutate(
        time_til = ifelse(is.finite(TreatStartTime), Time - TreatStartTime, -999), 
        # construct lead and lag columns (in the range of our defined horizon)
        lead = ifelse(time_til < 0, -1*time_til, "lag"), 
        lag = ifelse(time_til >= 0, time_til, "lead")) 
    df
  })
  
  # reactive data frame of the study data, with all dummy columns; will be used for fitting models and color decomposition 
  input_data_augment <- reactive({
    df_augment = input_data() %>% 
      fastDummies::dummy_cols(select_columns = c("Unit", "Time", "lead", "lag")) %>% 
      dplyr::select(-c("lead_lag", "lag_lead", "lead", "lag")) %>% 
      dplyr::rename(lag_inf = lead_999)
    df_augment
  })
  
  # reactive data frame of the study data with balanced format; will be used for panel plotting 
  panel_data <- reactive({
    panel_data = data.frame(
      cbind(
        rep(unique(input_data()$Unit), each = length(unique(input_data()$Time))), 
        rep(unique(input_data()$Time), length(unique(input_data()$Unit)))
      )
    )
    colnames(panel_data) = c("Unit", "Time")
    
    panel_data = panel_data %>% 
      mutate(Time = as.numeric(Time), 
             Unit = factor(Unit)) %>% 
      left_join(., input_data(), by = c("Unit", "Time")) %>% 
      dplyr::select(X, Unit, Time, Outcome, Treatment, TreatStartTime, time_til) %>% 
      mutate(Unit = factor(Unit, level = unique(Unit[order(TreatStartTime, decreasing = T)]) )) 
    panel_data
  })
  
  # reactive data frame of the study data with calculated sample influences; will be used for SIC plotting 
  inf_res = reactive({ 
    data = GetImpliedWeights(input_data_augment(), t0 = input$t1_5, t1 = input$ty_5, l_min = input$l_min, l_max = input$l_max)
    data = GetObsGroup(data = data, t0 = input$t1_5, t1 = input$ty_5, estimand = "std")
    GetInfluence(data = data, t0 = input$t1_5, t1 = input$ty_5, l_min = input$l_min, l_max = input$l_max) 
  })
  
  # Panel plotting 
  output$panel_fig = renderPlotly({
    panel_fig = ggplot(panel_data(), 
                       aes(x = Time, y = Unit, label=Outcome, label2=TreatStartTime), position="identity")  +
      geom_tile(color="gray80", fill="white", size=0.1, stat="identity") + 
      geom_point(aes(x = Time, y = Unit, shape = Treatment)) + 
      scale_shape_manual(values = c(19, 1)) + 
      labs(x = "Time", y = "Unit", 
           title = "", 
           shape = "Treatment Status"
      ) + panelview_plot_tyle
    ggplotly(panel_fig, height = 600, width = 880) %>% 
      plotly::layout(legend=list(xanchor='center', yanchor='bottom', y = -0.35, x = 0.5, orientation='h'))
  })
  
  # Panel plotting 
  output$treat_prop_fig = renderPlot({
    
    treat_prop_fig = panel_data() %>% 
      dplyr::group_by(Time) %>% 
      summarise(Treat = mean(Treatment == "1"), 
                Untreat = 1-Treat) %>% ungroup() %>% 
      gather(.,
             key = "TreatmentStatus",
             value = "Proportion",
             Treat, Untreat) %>% 
      ggplot() + 
      geom_bar(aes(x=Time, y=Proportion, fill=TreatmentStatus), position="fill", stat="identity") + 
      scale_fill_manual(values = c("#9dc4e4", "#c6e0b4")) +
      labs(x = "Time", 
           y = "Proportion", 
           title = "", 
           fill = "Treatment Status"
      ) 
    treat_prop_fig
  }, height = 600, width = 880)
  
  
  output$panelview = renderPlot({
    
    if(input$estimand1 == "new") {
      # define the valid control lead lag range
      time_til_ubound <- input$ty_1-input$t1_1
      time_til_lbound <- input$ty_1-max(panel_data()$Time)
      
      if (input$invariance) {
        panelview_data = panel_data() %>% 
          mutate(group_ind = ifelse(time_til == input$ty_1-input$t1_1, "Treat", 
                                    ifelse(between(time_til, time_til_lbound, time_til_ubound) | time_til==-999, "Control", "Invalid")))
        if (input$anticipation) {
          if (input$washout) {
            panelview_data = panel_data() %>% 
              mutate(group_ind = ifelse(time_til == input$ty_1-input$t1_1, "Treat", 
                                        ifelse(time_til > max(time_til_lbound:time_til_ubound), "Control", 
                                               ifelse(between(time_til, time_til_lbound, time_til_ubound) | time_til==-999, "Control", "Control"))))
          } else {
            panelview_data = panel_data() %>% 
              mutate(group_ind = ifelse(time_til == input$ty_1-input$t1_1, "Treat", 
                                        ifelse(time_til > max(time_til_lbound:time_til_ubound), "Invalid", 
                                               ifelse(between(time_til, time_til_lbound, time_til_ubound) | time_til==-999, "Control", "Control"))))
          }
        } else {
          if (input$washout) {
            panelview_data = panel_data() %>% 
              mutate(group_ind = ifelse(time_til == input$ty_1-input$t1_1, "Treat", 
                                        ifelse(time_til > max(time_til_lbound:time_til_ubound), "Control", 
                                               ifelse(between(time_til, time_til_lbound, time_til_ubound) | time_til==-999, "Control", "Invalid"))))
          } else {
            panelview_data = panel_data() %>% 
              mutate(group_ind = ifelse(time_til == input$ty_1-input$t1_1, "Treat", 
                                        ifelse(time_til > max(time_til_lbound:time_til_ubound), "Invalid", 
                                               ifelse(between(time_til, time_til_lbound, time_til_ubound) | time_til==-999, "Control", "Invalid"))))
          }
        }
      } else {
        panelview_data = panel_data() %>% 
          mutate(group_ind = ifelse(TreatStartTime == input$t1_1 & Time == input$ty_1, "Treat", 
                                    ifelse(Time == input$ty_1 & (between(time_til, time_til_lbound, time_til_ubound) | time_til==-999), "Control", "Invalid")))
      }
      
      if (!input$highlight) {
        panelview_data = panelview_data %>% 
          mutate(Time = factor(Time), 
                 group_ind = factor(group_ind, levels = c("Treat", "Control", "Invalid")))
        
        panelview = ggplot(panelview_data, aes(x = Time, y = Unit), position = "identity")  +
          geom_tile(aes(fill = group_ind), colour="gray90", size=0.1, stat="identity") + 
          scale_fill_manual(values=obs_colors) + 
          geom_point(aes(x = Time, y = Unit, shape = Treatment)) + 
          scale_shape_manual(values = c(19, 1)) + 
          labs(x = "Time", y = "Unit", title = "", 
               fill = "Observation Group", 
               shape = "Treatment Status") + 
          panelview_plot_tyle
        
      } else {
        panelview_data = GetObsGroup(panel_data(), t0 = input$t1_1, t1 = input$ty_1, estimand = "new") 
        panelview_data = panelview_data %>% 
          mutate(Time = factor(Time))
        
        panelview = ggplot(panelview_data, aes(x = Time, y = Unit), position = "identity")  +
          geom_tile(aes(fill = group_ind), colour="gray90", size=0.1, stat="identity") + 
          scale_fill_manual(values=ess_obs_colors) + 
          geom_point(aes(x = Time, y = Unit, shape = Treatment)) + 
          scale_shape_manual(values = c(19, 1)) + 
          labs(x = "Time", y = "Unit", title = "", 
               fill = "Observation Group", 
               shape = "Treatment Status") + 
          panelview_plot_tyle + 
          guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2))
      } 
    } else {
      panelview_data = panel_data() %>% 
        mutate(group_ind = ifelse(TreatStartTime == input$t1_1 & Time == input$ty_1, "Treat", 
                                  ifelse(TreatStartTime == Inf & Time == input$ty_1, "Control", "Invalid")))
      if (input$invariance) {
        panelview_data = panelview_data %>% 
          mutate(group_ind = ifelse(time_til == input$ty_1-input$t1_1, "Treat", 
                                    ifelse(time_til == -999, "Control", "Invalid")))
        if (input$anticipation) {
          panelview_data = panelview_data %>% 
            mutate(group_ind = ifelse(time_til == input$ty_1-input$t1_1, "Treat", 
                                      ifelse(Treatment == 0, "Control", "Invalid")))
          if (input$washout) {
            panelview_data = panelview_data %>% 
              mutate(group_ind = ifelse(time_til == input$ty_1-input$t1_1, "Treat", 
                                        ifelse((Treatment == 0 | time_til > input$ty_1-input$t1_1), "Control", "Invalid")))
            if (input$delay) {
              panelview_data = panelview_data %>% 
                mutate(group_ind = ifelse(time_til == input$ty_1-input$t1_1, "Treat", "Control"))
            }
          } else {
            if (input$delay) {
              panelview_data = panelview_data %>% 
                mutate(group_ind = ifelse(time_til == input$ty_1-input$t1_1, "Treat", 
                                          ifelse(time_til < input$ty_1-input$t1_1, "Control", "Invalid")))
            }
          }
        } else {
          # if limited anticipation does not hold 
          if (input$washout) {
            panelview_data = panelview_data %>% 
              mutate(group_ind = ifelse(time_til == input$ty_1-input$t1_1, "Treat", 
                                        ifelse(time_til > input$ty_1-input$t1_1 | time_til == -999, "Control", "Invalid")))
            if (input$delay) {
              panelview_data = panelview_data %>% 
                mutate(group_ind = ifelse(time_til == input$ty_1-input$t1_1, "Treat", 
                                          ifelse(Treatment == 1 | time_til == -999, "Control", "Invalid")))
            }
          } else {
            if (input$delay) {
              panelview_data = panelview_data %>% 
                mutate(group_ind = ifelse(time_til == input$ty_1-input$t1_1, "Treat", 
                                          ifelse((Treatment == 1 & time_til < input$ty_1-input$t1_1) | time_til == -999, "Control", "Invalid")))
            }
          }
        }
      }
      
      if (input$highlight) {
        panelview_data = GetObsGroup(panel_data(), t0 = input$t1_1, t1 = input$ty_1, estimand = "std")
        panelview_data = panelview_data %>% 
          mutate(Time = factor(Time))
        panelview = ggplot(panelview_data,
                           aes(x = Time, y = Unit), position = "identity")  +
          geom_tile(aes(fill = group_ind), colour="gray90", size=0.1, stat="identity") +
          scale_fill_manual(values=ess_obs_colors) +
          geom_point(data = panelview_data, aes(x = Time, y = Unit, shape = Treatment)) +
          scale_shape_manual(values = c(19, 1)) +
          labs(x = "Time", y = "Unit", title = "",
               fill = "Observation Group",
               shape = "Treatment Status") +
          panelview_plot_tyle
      } else {
        panelview_data = panelview_data %>% 
          mutate(Time = factor(Time), 
                 group_ind = factor(group_ind, levels = c("Treat", "Control", "Invalid")))
        panelview = ggplot(panelview_data, aes(x = Time, y = Unit), position = "identity")  +
          geom_tile(aes(fill = group_ind), colour="gray90", size=0.1, stat="identity") +
          scale_fill_manual(values=obs_colors) +
          geom_point(data = panelview_data, aes(x = Time, y = Unit, shape = Treatment)) +
          scale_shape_manual(values = c(19, 1)) +
          labs(x = "Time", y = "Unit", title = "",
               fill = "Observation Group",
               shape = "Treatment Status") +
          panelview_plot_tyle
      }
    }
    panelview}, 
    height = 600, width = 880)
  
  
  output$twfe_event_study_plot = renderPlot({
    event_study_formula <- as.formula(
      paste("Outcome ~",
            paste(
              paste(paste("lead_", 1:input$l_min, sep = ""), collapse = " + "),
              paste(paste("lag_", 0:input$l_max, sep = ""), collapse = " + "), sep = " + "),
            "| Time + Unit | 0 | Unit"
      ),
    )
    event_study_reg <- felm(event_study_formula, data = input_data_augment())
    # order of the coefficients for the plot
    plot_order <- c(paste("lead_", input$l_min:1, sep = ""), paste("lag_", 0:input$l_max, sep = ""))
    coef_lead1 <- coef(event_study_reg)["lead_1"]
    
    # grab the clustered standard errors and average coefficient estimates
    # from the regression, label them accordingly add a zero'th lag for plotting purposes
    leadslags_plot <- tibble(
      sd = c(event_study_reg$cse[plot_order]),
      mean = c(coef(event_study_reg)[plot_order]),
      label = c(-input$l_min:-1, 0:input$l_max)
    )
    # make the sd of lead -1 (reference group) as NA 
    leadslags_plot$sd[leadslags_plot$label == -1] = NA
    # An interval that traces the confidence intervals of coefficients
    p = leadslags_plot %>%
      ggplot(aes(x = label, y = mean-coef_lead1,
                 ymin = mean-coef_lead1-1.96*sd, 
                 ymax = mean-coef_lead1+1.96*sd)) +
      geom_line() + 
      geom_errorbar(col = "gray70", width = 0.5, lty = 1) +
      geom_point() +
      theme_bw() +
      labs(x = "Relative Time Periods", y = "Estimated Effects", title = "Event Study Plot by TWFE Model") +
      geom_hline(yintercept = 0, lty = 2) +
      geom_vline(xintercept = -1, lty = 2) + 
      theme(axis.title.x = element_text(size = 12, margin = margin(t = 8, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 8, b = 0, l = 0)),
            axis.text.x = element_text(face="bold", size = 10, angle = 0, hjust=0.5, vjust=0),
            axis.text.y = element_text(face="bold", size = 10),
            plot.title = element_text(size=15, hjust = 0.5, face="bold", margin = margin(8, 0, 8, 0)))
    p
  }, height = 600, width = 880)
  
  output$twfe_panelview = renderPlot({
    panelview_data = panel_data() %>% 
      mutate(group_ind = ifelse(time_til == input$ty_6-input$t1_6, "Treat", "Control"), 
             Time = factor(Time),
             group_ind = factor(group_ind, levels = c("Treat", "Control", "Invalid")))
    
    panelview = ggplot(panelview_data, aes(x = Time, y = Unit), position = "identity")  +
      geom_tile(aes(fill = group_ind), colour="gray90", size=0.1, stat="identity") +
      scale_fill_manual(values=obs_colors) +
      geom_point(data = panelview_data, aes(x = Time, y = Unit, shape = Treatment)) +
      scale_shape_manual(values = c(19, 1)) +
      labs(x = "Time", y = "Unit", title = "Treatment and Control Component Implied by TWFE Model",
           fill = "",
           shape = "Treatment Status") +
      panelview_plot_tyle
    panelview
  }, height = 600, width = 880)
  
  output$ess_plot = renderPlot({
    panelview_data = GetObsGroup(panel_data(), t0 = input$t1_4, t1 = input$ty_4, estimand = "std") 
    panelview_data = panelview_data %>% 
      mutate(Time = factor(Time))
    
    panelview = ggplot(panelview_data,
                       aes(x = Time, y = Unit), position = "identity")  +
      geom_tile(aes(fill = group_ind), colour="gray90", size=0.1, stat="identity") +
      scale_fill_manual(values=ess_obs_colors) +
      geom_point(data = panelview_data, aes(x = Time, y = Unit, shape = Treatment)) +
      scale_shape_manual(values = c(19, 1)) +
      labs(x = "Time", y = "Unit", title = "", fill = "Observation Group", shape = "Treatment Status") +
      panelview_plot_tyle
    panelview
  }, height = 600, width = 880)
  
  output$twfe_ess = renderTable({
    data_augment_weight = GetImpliedWeights(data = input_data_augment(), t0 = input$t1_4, t1 = input$ty_4, l_min = input$l_min, l_max = input$l_max)
    data = GetObsGroup(data = data_augment_weight, t0 = input$t1_4, t1 = input$ty_4, estimand = "std")
    covar_names = str_split(input$covar_name, ",")[[1]]
    GetESS(data = data, t0 = input$t1_4, t1 = input$ty_4, bal_covariates = covar_names, method = "twfe")
  }, 
  striped = TRUE,
  spacing = "l",
  digits = 3,
  rownames = TRUE, 
  width = "100%",
  caption = ""
  )
  
  output$balance_experiment_ess = renderTable({
    data_augment_weight = GetImpliedWeights(data = input_data_augment(), t0 = input$t1_4, t1 = input$ty_4, l_min = input$l_min, l_max = input$l_max)
    data = GetObsGroup(data = data_augment_weight, t0 = input$t1_4, t1 = input$ty_4, estimand = "std")
    covar_names = str_split(input$covar_name, ",")[[1]]
    GetESS(data = data, t0 = input$t1_4, t1 = input$ty_4, bal_covariates = covar_names, method = "experiment")
  }, 
  striped = TRUE,
  spacing = "l",
  digits = 3,
  rownames = TRUE, 
  width = "100%",
  caption = ""
  )
  
  output$balance_invariance_ess = renderTable({
    data_augment_weight = GetImpliedWeights(data = input_data_augment(), t0 = input$t1_4, t1 = input$ty_4, l_min = input$l_min, l_max = input$l_max)
    data = GetObsGroup(data = data_augment_weight, t0 = input$t1_4, t1 = input$ty_4, estimand = "std")
    covar_names = str_split(input$covar_name, ",")[[1]]
    GetESS(data = data, t0 = input$t1_4, t1 = input$ty_4, bal_covariates = covar_names, method = "invariance")
  }, 
  striped = TRUE,
  spacing = "l",
  digits = 3,
  rownames = TRUE, 
  width = "100%",
  caption = ""
  )
  
  output$balance_anticipate_ess = renderTable({
    data_augment_weight = GetImpliedWeights(data = input_data_augment(), t0 = input$t1_4, t1 = input$ty_4, l_min = input$l_min, l_max = input$l_max)
    data = GetObsGroup(data = data_augment_weight, t0 = input$t1_4, t1 = input$ty_4, estimand = "std")
    covar_names = str_split(input$covar_name, ",")[[1]]
    GetESS(data = data, t0 = input$t1_4, t1 = input$ty_4, bal_covariates = covar_names, method = "anticipate")
  }, 
  striped = TRUE,
  spacing = "l",
  digits = 3,
  rownames = TRUE, 
  width = "100%",
  caption = ""
  )
  
  output$balance_delay_ess = renderTable({
    data_augment_weight = GetImpliedWeights(data = input_data_augment(), t0 = input$t1_4, t1 = input$ty_4, l_min = input$l_min, l_max = input$l_max)
    data = GetObsGroup(data = data_augment_weight, t0 = input$t1_4, t1 = input$ty_4, estimand = "std")
    covar_names = str_split(input$covar_name, ",")[[1]]
    GetESS(data = data, t0 = input$t1_4, t1 = input$ty_4, bal_covariates = covar_names, method = "delay")
  }, 
  striped = TRUE,
  spacing = "l",
  digits = 3,
  rownames = TRUE, 
  width = "100%",
  caption = ""
  )
  
  output$balance_dissipate_ess = renderTable({
    data_augment_weight = GetImpliedWeights(data = input_data_augment(), t0 = input$t1_4, t1 = input$ty_4, l_min = input$l_min, l_max = input$l_max)
    data = GetObsGroup(data = data_augment_weight, t0 = input$t1_4, t1 = input$ty_4, estimand = "std")
    covar_names = str_split(input$covar_name, ",")[[1]]
    GetESS(data = data, t0 = input$t1_4, t1 = input$ty_4, bal_covariates = covar_names, method = "dissipate")
  }, 
  striped = TRUE,
  spacing = "l",
  digits = 3,
  rownames = TRUE, 
  width = "100%",
  caption = ""
  )
  
  output$infuence_plot_panel = renderPlot({
    panelview_data = GetObsGroup(panel_data(), t0 = input$t1_5, t1 = input$ty_5, estimand = "std") 
    panelview_data = panelview_data %>% 
      mutate(Time = factor(Time))
    
    panelview = ggplot(panelview_data,
                       aes(x = Time, y = Unit), position = "identity")  +
      geom_tile(aes(fill = group_ind), colour="gray90", size=0.1, stat="identity") +
      scale_fill_manual(values=ess_obs_colors) +
      geom_point(data = panelview_data, aes(x = Time, y = Unit, shape = Treatment)) +
      scale_shape_manual(values = c(19, 1)) +
      labs(x = "Time", y = "Unit", title = "", fill = "Observation Group", shape = "Treatment Status") +
      panelview_plot_tyle
    panelview
  }, height = 600, width = 880)
  
  output$infuence_plot = renderPlotly({
    
    data = GetImpliedWeights(input_data_augment(), t0 = input$t1_5, t1 = input$ty_5, l_min = input$l_min, l_max = input$l_max)
    data = GetObsGroup(data, t0 = input$t1_5, t1 = input$ty_5, estimand = "std")
    
    if (input$influence_metric == "sic") {
      change_plt = PlotInfluence(data = data, inf_all = inf_res(), input$t1_5, input$ty_5, metric = "sic")
      ggplotly(change_plt, tooltip = c("group_ind", "Unit", "Time", "TreatStartTime","Outcome","SIC"))
    } else if (input$influence_metric == "sic_scaled") {
      change_plt = PlotInfluence(data = data, inf_all = inf_res(), input$t1_5, input$ty_5, metric = "sic_scaled")
      ggplotly(change_plt, tooltip = c("group_ind", "Unit", "Time", "TreatStartTime","Outcome","SIC.Scaled"))
    } else {
      change_plt = PlotInfluence(data = data, inf_all = inf_res(), input$t1_5, input$ty_5, metric = "est_change")
      ggplotly(change_plt, tooltip = c("group_ind", "Unit", "Time", "TreatStartTime","Outcome","Est.Change"))
    }
  })
  
  output$implied_unit_plot = renderPlotly({
    data_augment_weight = GetImpliedWeights(data = input_data_augment(), t0 = input$t1_2, t1 = input$ty_2, l_min = input$l_min, l_max = input$l_max)
    
    if (input$implied_unit_by_time) {
      implied_unit_plot = data_augment_weight %>% 
        dplyr::select(X, Unit, Time, Treatment, lmw_weight, time_til, treatment_component) %>%
        dplyr::filter(Time == as.numeric(input$time_value)) %>% 
        group_by(Unit, treatment_component) %>%
        summarise(Proportion = round(mean(lmw_weight), 4)) %>% ungroup() %>%
        ggplot(aes(x = Unit, xend = Unit, y = 0, yend = Proportion)) +
        geom_segment(aes(col = treatment_component), linetype = "solid", linewidth = 0.3) + 
        scale_color_manual(values=c("1" = "#3293e3", "0" = "#69b334")) +
        labs(x = "Unit", y = "Implied Proportion", title = paste0("Implied Unit Proportion of Time ", input$time_value)) +
        # geom_text(aes(x = Unit, y = Proportion, label = Proportion), size = 5, vjust = 0) +
        scale_x_discrete(breaks=unique(data_augment_weight$Unit)) + 
        theme_bw() + 
        theme(plot.title = element_text(size = 15, face = "bold"), 
              axis.title = element_text(size = 12), 
              axis.text = element_text(size = 10, face = "bold"))
    } else {
      implied_unit_plot = data_augment_weight %>% 
        dplyr::select(X, Unit, Time, Treatment, lmw_weight, time_til, treatment_component) %>%
        group_by(Unit) %>%
        summarise(Proportion = round(mean(lmw_weight), 4)) %>% ungroup() %>%
        ggplot(aes(x = Unit, xend = Unit, y = 0, yend = Proportion)) +
        geom_segment(linetype = "solid", linewidth = 0.3) +  
        labs(x = "Unit", y = "Implied Proportion", title = "Implied Unit Proportion") +
        theme_bw() + 
        theme(plot.title = element_text(size = 15, face = "bold"), 
              axis.title = element_text(size = 12), 
              axis.text = element_text(size = 10, face = "bold"))
    }
    ggplotly(implied_unit_plot, tooltip = c("Unit", "Proportion"), height = 600, width = 900) 
  })
  
  
  output$implied_time_plot = renderPlotly({
    data_augment_weight = GetImpliedWeights(data = input_data_augment(), t0 = input$t1_2, t1 = input$ty_2, l_min = input$l_min, l_max = input$l_max)
    if (input$implied_time_by_unit) {
      implied_time_plot = data_augment_weight %>% 
        dplyr::select(X, Unit, Time, Treatment, lmw_weight, time_til, treatment_component) %>%
        dplyr::filter(Unit == input$unit_value) %>% 
        group_by(Time, treatment_component) %>%
        summarise(Proportion = round(mean(lmw_weight), 4)) %>% ungroup() %>%
        ggplot(aes(x = Time, xend = Time, y = 0, yend = Proportion)) +
        geom_segment(aes(col = treatment_component), linetype = "solid", linewidth = 0.3) + 
        scale_color_manual(values=c("1" = "#3293e3", "0" = "#69b334")) +
        labs(x = "Time", y = "Implied Proportion", title = paste0("Implied Time Proportion of Unit ", input$unit_value)) +
        theme_bw() + 
        theme(plot.title = element_text(size = 15, face = "bold"), 
              axis.title = element_text(size = 12), 
              axis.text = element_text(size = 10, face = "bold"))
    } else {
      implied_time_plot = data_augment_weight %>% 
        dplyr::select(X, Unit, Time, Treatment, lmw_weight, time_til, treatment_component) %>%
        group_by(Time) %>%
        summarise(Proportion = round(mean(lmw_weight), 3)) %>% ungroup() %>%
        ggplot(aes(x = Time, xend = Time, y = 0, yend = Proportion)) +
        geom_segment(linetype = "solid", linewidth = 0.3) +  
        labs(x = "Time", y = "Implied Proportion", title = "Implied Time Proportion") +
        theme_bw() + 
        theme(plot.title = element_text(size = 15, face = "bold"), 
              axis.title = element_text(size = 12), 
              axis.text = element_text(size = 10, face = "bold"))
    }
    ggplotly(implied_time_plot, tooltip = c("Time", "Proportion"), height = 600, width = 900) 
  })
  
  output$implied_map_plot = renderPlot({
    ################ Map Data ################
    # Map state names and state abbreviations
    state_names = str_to_lower(state.name)
    state_abbr = state.abb
    state_info = data.frame(cbind(state = state_abbr, state_names = state_names)) %>%
      filter(!state %in% c("AK", "HI"))
    
    # Longitude and Latitude of state center (for plotting)
    state_center = tibble(
      state = datasets::state.abb,
      c_long = datasets::state.center$x,
      c_lat = datasets::state.center$y
    ) %>%
      filter(!state %in% c("AK", "HI")) %>%
      left_join(state_info, "state")
    
    # Complete map data with state abbv. & Longitude and Latitude of state center
    state_data = map_data("state") %>%
      dplyr::select(long, lat, group, region) %>%
      dplyr::rename(state_names = region) %>%
      left_join(state_center, by = "state_names") %>%
      filter(!state %in% c("AK", "HI")) %>%
      distinct() %>%
      dplyr::select(state, long, lat, c_long, c_lat, group)
    
    ################ Implied Weighted States ################
    # Join the two datasets: map data and the value (lmw weight) for each observation
    
    data_augment_weight = GetImpliedWeights(data = input_data_augment(), t0 = input$t1_2, t1 = input$ty_2, l_min = input$l_min, l_max = input$l_max)
    map_plot_data = state_data %>% 
      left_join(data_augment_weight, by = c("state" = "Unit")) %>%
      dplyr::rename(Unit = state)
    
    MapPlot = function(by_time, time_value, shade) {
      if (by_time) {
        map_center_plot_data = map_plot_data %>% 
          filter(Time == time_value) %>% 
          dplyr::select(Time, Unit, group, c_long, c_lat, Treatment) %>% 
          mutate(Treatment = factor(Treatment, level = c(0, 1))) %>% 
          distinct()
        
        map_year_plot_data = map_plot_data %>%
          filter(Time == time_value) %>%
          dplyr::select(X, Time, Unit, long, lat, group, c_long, c_lat,
                        lmw_weight, time_til, Treatment, treatment_component) %>%
          mutate(Treatment = factor(Treatment, level = c(0, 1)))
        
        # Plot
        if (shade) {
          map_plot =
            ggplot() +
            geom_polygon(data=map_year_plot_data[map_year_plot_data$treatment_component == 1, ],
                         aes(x=long, y=lat, group=group, fill=lmw_weight), color="black", size = 0.2) +
            scale_fill_continuous(name= paste0("Treated observation implied weight"),
                                  low = "#9dc4e4", high = "#3787c8",
                                  limits = c(
                                    min(map_plot_data[map_plot_data$treatment_component == 1, "lmw_weight"]),
                                    max(map_plot_data[map_plot_data$treatment_component == 1, "lmw_weight"])),
                                  na.value = "grey50") +
            new_scale_fill() +
            geom_polygon(data=map_year_plot_data[map_year_plot_data$treatment_component == 0,],
                         aes(x=long, y=lat, group=group, fill=lmw_weight), color="black", size = 0.2) +
            scale_fill_gradient2(name= paste0("Control observation implied weight"),
                                 low = "darkgreen", mid = "white", high = "darkgreen",
                                 limits = c(
                                   min(map_plot_data[map_plot_data$treatment_component == 0, "lmw_weight"]),
                                   max(map_plot_data[map_plot_data$treatment_component == 0, "lmw_weight"])),
                                 na.value = "grey50") +
            geom_point(data = map_center_plot_data[map_center_plot_data$Treatment == 1, ],
                       aes(x = c_long, y = c_lat), size = 10, shape = 19) +
            geom_text(data = map_center_plot_data[map_center_plot_data$Treatment == 1, ],
                      aes(c_long, c_lat, label = Unit), color = "white", size = 5) +
            geom_text(data = map_center_plot_data[map_center_plot_data$Treatment == 0, ],
                      aes(c_long, c_lat, label = Unit), color = "black", size = 5) +
            labs(title = paste0("Time ", time_value)) +
            ggthemes::theme_map() + implied_plot_tyle
        } else {
          map_plot =
            ggplot() +
            geom_polygon(data=map_year_plot_data[map_year_plot_data$treatment_component == 1, ],
                         aes(x=long, y=lat, group=group), fill = "#9dc4e4", color="black", size = 0.2) +
            new_scale_fill() +
            geom_polygon(data=map_year_plot_data[map_year_plot_data$treatment_component == 0,],
                         aes(x=long, y=lat, group=group), fill="#c6e0b4", color="black", size = 0.2) +
            geom_point(data = map_center_plot_data[map_center_plot_data$Treatment == 1, ],
                       aes(x = c_long, y = c_lat), size = 10, shape = 19) +
            geom_text(data = map_center_plot_data[map_center_plot_data$Treatment == 1, ],
                      aes(c_long, c_lat, label = Unit), color = "white", size = 5) +
            geom_text(data = map_center_plot_data[map_center_plot_data$Treatment == 0, ],
                      aes(c_long, c_lat, label = Unit), color = "black", size = 5) +
            labs(title = paste0("Time ", time_value)) +
            ggthemes::theme_map() + implied_plot_tyle
        }
      } else {
        
        map_center_plot_data = map_plot_data %>% 
          dplyr::select(Time, Unit, group, c_long, c_lat, Treatment) %>% 
          mutate(Treatment = factor(Treatment, level = c(0, 1))) %>% 
          distinct()
        
        map_year_plot_data = map_plot_data %>%
          dplyr::select(X, Time, Unit, long, lat, group, c_long, c_lat,
                        lmw_weight, time_til, Treatment, treatment_component) %>%
          mutate(Treatment = factor(Treatment, level = c(0, 1)))
        
        if (shade) {
          map_plot = map_plot_data %>%
            na.omit() %>%
            dplyr::select(X, Time, Unit, long, lat, group, c_long, c_lat,
                          lmw_weight, time_til, Treatment, treatment_component) %>%
            filter(treatment_component == 1) %>%
            group_by(Unit) %>%
            summarise(lmw_weight = mean(lmw_weight)) %>% ungroup() %>%
            right_join(state_data, by = c("Unit" = "state")) %>%
            ggplot() +
            geom_polygon(aes(x=long, y=lat, group=group, fill=lmw_weight), color="black", size = 0.2) +
            scale_fill_continuous(name= paste0("Observation implied weight"),
                                  low = "#9dc4e4", high = "#3787c8", na.value = "white") +
            geom_text(data = map_plot_data, aes(c_long, c_lat, label = Unit),
                      color = "black", size = 5) +
            ggthemes::theme_map() + implied_plot_tyle
        } else {
          map_plot =
            state_data %>%
            na.omit() %>%
            ggplot() +
            geom_polygon(aes(x=long, y=lat, group=group), fill="white", color="black", size = 0.2) +
            geom_text(data = state_data, aes(c_long, c_lat, label = state),
                      color = "black", size = 5) +
            ggthemes::theme_map() + implied_plot_tyle  
        }
      }
      return(map_plot)
    }
    MapPlot(by_time = input$implied_unit_by_time, time_value = as.numeric(input$time_value), shade=input$shade)
  })
  
  lmw_output_download <- reactive({
    data_augment_weight = GetImpliedWeights(data = input_data_augment(), t0 = input$t1_2, t1 = input$ty_2, l_min = input$l_min, l_max = input$l_max)
    data = GetObsGroup(data = data_augment_weight, t0 = input$t1_2, t1 = input$ty_2, estimand = input$estimand2)
    
    covar_names = str_split(input$covar_name, ",")[[1]]
    col_names = paste("X", 1:length(covar_names), sep = "")
    data_output = data %>% 
      dplyr::select(Unit, Time, Outcome, Treatment, any_of(col_names), lmw_weight, group_ind) 
    colnames(data_output) = c(input$unit_ind_name, input$time_ind_name, 
                              input$outcome_name, input$treat_name, covar_names, 
                              "implied.weights", "obs.group")
    data_output
  })
  
  inf_output_download <- reactive({
    inf_output = inf_res()
    user_df = input_data()
    
    covar_names = str_split(input$covar_name, ",")[[1]]
    col_names = paste("X", 1:length(covar_names), sep = "")
    inf_output = user_df %>% 
      left_join(inf_output, by = "X") %>% 
      dplyr::select(Unit, Time, Outcome, Treatment, any_of(col_names), inf, inf_scaled, est_change) 
    colnames(inf_output) = c(input$unit_ind_name, input$time_ind_name, 
                             input$outcome_name, input$treat_name, covar_names, 
                             "SIC", "SIC.Scaled", "Est.Change") 
    inf_output
  })
  
  output$download_lmw <- downloadHandler(
    filename = function() { paste("lmw-", Sys.Date(), ".csv", sep="") },
    content = function(file) { write.csv( lmw_output_download(), file) }
  )
  
  output$download_inf <- downloadHandler(
    filename = function() { paste("inf-", Sys.Date(), ".csv", sep="") },
    content = function(file) { write.csv(inf_output_download(), file) }
  )
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
