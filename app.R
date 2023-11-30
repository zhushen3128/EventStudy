#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(bacondecomp)
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

library(shiny)
library(sf)
library(leaflet)
library(shinyTree)

#######################################################################
##    Global 
#######################################################################

obs_colors <- c(treat = "#9dc4e4", control = "#c6e0b4", 
                early ="#f7a59b", late = "#ffc001", invalid = "white")


panelview_plot_tyle <- 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA,color="white", size=0.5, linetype="solid"),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title=element_text(size=12),
          axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
          axis.text = element_text(color="black", size=8),
          axis.text.x = element_text(size = 8, angle = 0, hjust=0.5, vjust=0),
          axis.text.y = element_text(size = 8),
          plot.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white"),
          legend.position = "bottom",
          legend.margin = margin(c(0, 5, 5, 0)),
          legend.text = element_text(margin = margin(r = 10, unit = "pt"), size = 12),
          plot.title = element_text(size=15, hjust = 0.5, face="bold",
                                    margin = margin(8, 0, 8, 0)))

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


# define the entire follow-up period 
lead_min = 20 
lag_max = 20

# read and preprocess data 
dat <- read.csv("bacon_example.csv")
data_study <- dat %>% 
    # for each state, contruct the baseline covariates 
    group_by(stfips) %>%
    summarise(
        pcinc_baseline = pcinc[which.min(year)], 
        asmrh_baseline = asmrh[which.min(year)], 
        cases_baseline = cases[which.min(year)], 
        # construct column for treatment start year 
        treat_start_year = ifelse(sum(treat) == 0, 9999, year[time_to_treat == 0])
    ) %>% ungroup() %>% 
    right_join(dat, by = "stfips") %>% 
    rename(time_til = time_to_treat) %>%  
    mutate(# replace with the state names 
        state = cdlTools::fips(stfips, to = "Abbreviation"), 
        # if time_til = -999 for never treated units 
        time_til = replace(time_til, (treat == 0 & time_til==0), -999), 
        # construct column for treatment status 
        treat_status = ifelse(time_til >= 0 & treat == 1, 1, 0), 
        # construct lead and lag columns (in the range of our defined horizon)
        lead = ifelse((time_til < 0) & (time_til >= -1*lead_min), -1*time_til, "lag"), 
        lag = ifelse((time_til >= 0) & (time_til <= lag_max), time_til, "lead"), 
        # one column correspond to never treated 
        lag_inf = as.factor(ifelse(time_til == -999, 1, 0)), 
        X = row_number()) 

data_study_augment = data_study %>% 
    fastDummies::dummy_cols(select_columns = c("state", "year", "lead", "lag")) %>% 
    dplyr::select(-c("lead_lag", "lag_lead"))  

# define names of variables we would use 
baseline_covariates <- c(
    "pcinc_baseline", "asmrh_baseline", "cases_baseline"
)

state_covariates <- c(
    paste("state", unique(data_study$state), sep = "_")
)

year_covariates <- c(
    paste("year", unique(data_study$year), sep = "_")
)

lead_lag_covariates <- c(
    paste("lead", seq(lead_min,1), sep = "_"), 
    paste("lag", seq(0,lag_max), sep = "_"), 
    "lag_inf"
)

#######################################################################
##    UI
#######################################################################

ui <- fluidPage(
    
    # Application title
    titlePanel("Visualize Event Study Dataset"),
    
    navbarPage("Visualizations",
               
               tabPanel("Estimand and Assumptions",
                        
                        sidebarLayout(
                            sidebarPanel(
                                h4("Estimand"), 
                                selectInput("estimand", "Choose the estimand:",
                                            c("ATE_{t0, >t0, t1}" = "new1",
                                              "ATE_{t0, inf, t1}" = "std1")), 
                                withMathJax(),
                                helpText('Estimand can be ATE\\(_{t_0, >t_0, t_1}\\), which is the average causal 
                      effect the target population observed at time \\(t_1\\) from adopting the treatment 
                      for the first time at time \\(t_0\\) to adopting the treatment sometime after time \\(t_0\\); '), 
                                helpText('A special case of the estimand is ATE\\(_{t_0, \\infty, t_1}\\), which is is the average causal 
                      effect in the target population observed at time \\(t_1\\) from adopting the treatment for the first 
                     time at time \\(t_0\\) to never adopting the treatment.'), 
                                
                                sliderInput("t0_1",
                                            "What is the t0 in the target estimand? ", 
                                            value = 1980, 
                                            min = 1964,
                                            max = 1996), 
                                withMathJax(),
                                helpText('\\(t_0\\) denotes the treatment initiation time.'),
                                
                                sliderInput("t1_1",
                                            "What is the t1 in the target estimand? ", 
                                            value = 1985, 
                                            min = 1964,
                                            max = 1996), 
                                withMathJax(),
                                helpText('\\(t_1\\) denotes the time of observation.'),
                                
                                # Horizontal line 
                                tags$hr(), 
                                h4("Assumptions"), 
                                helpText('Consider scenarios when the following assumptions hold or not hold.'),
                                
                                checkboxInput("homogeneity", 'Treatment Effect Homogeneity. ', FALSE), 
                                #withMathJax(),
                                #helpText('\\( E \\{ Y_{it_1}(t_0) - Y_{it_1}(\\infty) \\mid X_i \\}\\)'),
                                #helpText('\\( = E \\{ Y_{it_1+l}(t_0+l) - Y_{it_1+l}(\\infty) \\mid X_i \\} \\)'),
                                #helpText('for all \\( t_0, t_1. \\)'), 
                                
                                
                                # Only show this panel if the assumption of treatment effect homogeneity holds 
                                conditionalPanel(
                                    condition = "input.homogeneity",
                                    checkboxInput("anticipation", 'Limited Treatment Anticipation.', FALSE), 
                                    
                                    checkboxInput("washout", 'Washout of Treatment Effect.', FALSE)
                                ), 
                                
                                conditionalPanel(
                                    condition = "input.homogeneity",
                                    helpText('Represent groups of observation corresponding to the assumptions above.'),
                                    checkboxInput("highlight", 'Show observations in the four groups with different colors.', FALSE)
                                ), 
                                
                                # Horizontal line 
                                tags$hr(), 
                                
                                h5("Strength of assumptions"), 
                                conditionalPanel(
                                    condition = "input.anticipation",
                                    checkboxInput("anticipation_strength", 'Show strength of limited treatment anticipation assumption.', FALSE) 
                                ), 
                                conditionalPanel(
                                    condition = "input.washout",
                                    checkboxInput("washout_strength", 'Show strength of washout of treatment effect assumption.', FALSE)
                                )
                                #selectInput("strength", "Represent the strength of the anticipation and washout assumptions:",
                                #            c("default" = "default", 
                                #              "Strength of anticipation assumption" = "anticipation",
                                #              "Strength of washout assumption" = "washout", 
                                #              "Strength of both" = "both"))
                            ),
                            mainPanel(plotOutput("panelview"))
                        )
               ), 
               
               tabPanel("Implied Weights Population",
                        
                        sidebarLayout(
                            sidebarPanel(
                                h4("Estimand"), 
                                selectInput("estimand_2", "Choose the estimand:",
                                            #c("ATE_{t0, >t0, t1}" = "new",
                                            #  "ATE_{t0, inf, t1}" = "std")), 
                                            c("ATE_{t0, inf, t1}" = "std2")), 
                                withMathJax(),
                                helpText('Estimand can be ATE\\(_{t_0, >t_0, t_1}\\), which is the average causal 
                      effect the target population observed at time \\(t_1\\) from adopting the treatment 
                      for the first time at time \\(t_0\\) to adopting the treatment sometime after time \\(t_0\\); '), 
                                helpText('A special case of the estimand is ATE\\(_{t_0, \\infty, t_1}\\), which is is the average causal 
                      effect in the target population observed at time \\(t_1\\) from adopting the treatment for the first 
                     time at time \\(t_0\\) to never adopting the treatment.'), 
                                
                                sliderInput("t0_2",
                                            "What is the t0 in the target estimand? ", 
                                            value = 1980, 
                                            min = 1964,
                                            max = 1996), 
                                withMathJax(),
                                helpText('\\(t_0\\) denotes the treatment initiation time.'),
                                
                                sliderInput("t1_2",
                                            "What is the t1 in the target estimand? ", 
                                            value = 1985, 
                                            min = 1964,
                                            max = 1996), 
                                withMathJax(),
                                helpText('\\(t_1\\) denotes the time of observation.'),
                                
                                # Horizontal line 
                                tags$hr(), 
                                
                                
                                h4("Implied Weights Population"), 
                                
                                h5("Implied Weighted States"), 
                                
                                checkboxInput("implied_state", 
                                              'Show implied weighted states', FALSE), 
                                
                                conditionalPanel(
                                    condition = "input.implied_state",
                                    checkboxInput("implied_state_by_year", 
                                                  'Show implied weighted states by year', FALSE), 
                                    conditionalPanel(
                                        condition = "input.implied_state_by_year",
                                        sliderInput("year_value",
                                                    "Show implied weighted states by the following year:", 
                                                    value = 1980, 
                                                    min = 1964,
                                                    max = 1996)
                                    ), 
                                    checkboxInput("shade", 'Show differential implied weights (colored with shade)', FALSE)
                                ), 
                                
                                tags$hr(), 
                                h5("Implied Weighted Years"), 
                                checkboxInput("implied_year", 'Show implied weighted years', FALSE)
                            ), 
                            mainPanel(plotOutput("implied"))
                        )
               )
    )
)


#######################################################################
##    Backend Control 
#######################################################################

server = function(input, output){
    
    output$panelview = renderPlot({
        # define the valid control lead lag range
        time_til_ubound <- input$t1_1-input$t0_1
        time_til_lbound <- input$t1_1-max(data_study$year)
        
        # generate visualization of the panel data 
        panelview_data = data.frame(cbind(rep(unique(data_study$state), 
                                              each = length(unique(data_study$year))), 
                                          rep(unique(data_study$year), 
                                              length(unique(data_study$state)))))
        colnames(panelview_data) = c("state", "year")
        
        if(input$estimand == "new1") {
            
            panelview_data = panelview_data %>% 
                mutate(year = as.numeric(year)) %>% 
                left_join(., data_study, by = c("state", "year")) %>% 
                dplyr::select(state, year, treat_status, time_til, treat_start_year) %>% 
                mutate(state = factor(state, level = unique(state[order(treat_start_year, decreasing = T)]) ), 
                       year = factor(year, level = min(year):max(year)), 
                       treat_status = factor(treat_status, level = c(1, 0)))
            
            if (input$homogeneity) {
                panelview_data = panelview_data %>% 
                    mutate(four_group_ind = 
                               ifelse(time_til == input$t1_1-input$t0_1, "treat", 
                                      ifelse(between(time_til, time_til_lbound, time_til_ubound) 
                                             | time_til==-999, "control", "invalid")))
                if (input$anticipation) {
                    if (input$washout) {
                        panelview_data = panelview_data %>% 
                            mutate(four_group_ind = ifelse(time_til == input$t1_1-input$t0_1, "treat", 
                                                           ifelse(time_til > max(time_til_lbound:time_til_ubound), "control", 
                                                                  ifelse(between(time_til, time_til_lbound, time_til_ubound) 
                                                                         | time_til==-999, "control", "control"))))
                    } else {
                        panelview_data = panelview_data %>% 
                            mutate(four_group_ind = ifelse(time_til == input$t1_1-input$t0_1, "treat", 
                                                           ifelse(time_til > max(time_til_lbound:time_til_ubound), "invalid", 
                                                                  ifelse(between(time_til, time_til_lbound, time_til_ubound) 
                                                                         | time_til==-999, "control", "control"))))
                    }
                } else {
                    if (input$washout) {
                        panelview_data = panelview_data %>% 
                            mutate(four_group_ind = ifelse(time_til == input$t1_1-input$t0_1, "treat", 
                                                           ifelse(time_til > max(time_til_lbound:time_til_ubound), "control", 
                                                                  ifelse(between(time_til, time_til_lbound, time_til_ubound) 
                                                                         | time_til==-999, "control", "invalid"))))
                    } else {
                        panelview_data = panelview_data %>% 
                            mutate(four_group_ind = ifelse(time_til == input$t1_1-input$t0_1, "treat", 
                                                           ifelse(time_til > max(time_til_lbound:time_til_ubound), "invalid", 
                                                                  ifelse(between(time_til, time_til_lbound, time_til_ubound) 
                                                                         | time_til==-999, "control", "invalid"))))
                    }
                }
            } else {
                panelview_data = panelview_data %>% 
                    mutate(
                        four_group_ind = 
                            ifelse(treat_start_year == input$t0_1 & year == input$t1_1, "treat", 
                                   ifelse(year == input$t1_1 & (between(time_til, time_til_lbound, time_til_ubound) 
                                                              | time_til==-999), "control", "invalid")))
            }
            
            if (!input$highlight) {
                panelview = ggplot(panelview_data, 
                                   aes(x = year, y = state, fill = four_group_ind), position = "identity")  +
                    geom_tile(colour="gray90", size=0.1, stat="identity") + 
                    scale_fill_manual(values=obs_colors) + 
                    geom_point(data = panelview_data, aes(x = year, y = state, shape = treat_status)) + 
                    scale_shape_manual(values = c(19, 1)) + 
                    labs(x = "Year", y = "State", title = "", 
                         fill = "Observation Group", 
                         shape = "Treatment Status"
                    ) + panelview_plot_tyle
            } else {
                panelview_data = panelview_data %>%
                    mutate(
                        four_group_ind = ifelse(time_til == input$t1_1-input$t0_1, "treat",
                                                ifelse(time_til > max(time_til_lbound:time_til_ubound), "early",
                                                       ifelse(between(time_til, time_til_lbound, time_til_ubound)
                                                              | time_til==-999, "control", "late"))))
                
                panelview = ggplot(panelview_data, 
                                   aes(x = year, y = state, fill = four_group_ind), position = "identity")  +
                    geom_tile(colour="gray90", size=0.1, stat="identity") + 
                    scale_fill_manual(values=obs_colors) + 
                    geom_point(data = panelview_data, aes(x = year, y = state, shape = treat_status)) + 
                    scale_shape_manual(values = c(19, 1)) + 
                    labs(x = "Year", y = "State", title = "", 
                         fill = "Observation Group", 
                         shape = "Treatment Status"
                    ) + panelview_plot_tyle
            } 
            
            
            # Strength of the assumptions 
            panelview_data = panelview_data %>%
                mutate(
                    four_group_ind = ifelse(time_til == input$t1_1-input$t0_1, "treat",
                                            ifelse(time_til > max(time_til_lbound:time_til_ubound), "early",
                                                   ifelse(between(time_til, time_til_lbound, time_til_ubound)
                                                          | time_til==-999, "control", "late")))) 
            
            if (input$anticipation_strength) {
                if (input$washout_strength) {
                    panelview = ggplot() +
                        geom_tile(
                            data = panelview_data[panelview_data$four_group_ind == "late", ], 
                            aes(x = year, y = state, fill = time_til), position = "identity", 
                            colour="gray90", size=0.1, stat="identity") + 
                        scale_fill_continuous(name= paste0("Limited Treatment Anticipation"), 
                                              low = "white", high = "#ffc001",
                                              breaks=c(
                                                  min(panelview_data[panelview_data$four_group_ind == "late", "time_til"]),
                                                  max(panelview_data[panelview_data$four_group_ind == "late", "time_til"])),
                                              labels=c("Weaker","Stronger"), 
                                              limits = c(
                                                  min(panelview_data[panelview_data$four_group_ind == "late", "time_til"]),
                                                  max(panelview_data[panelview_data$four_group_ind == "late", "time_til"])),
                                              na.value = "white") + 
                        new_scale_fill() +
                        geom_tile(
                            data = panelview_data[panelview_data$four_group_ind == "early", ], 
                            aes(x = year, y = state, fill = time_til), position = "identity", 
                            colour="gray90", size=0.1, stat="identity") + 
                        scale_fill_continuous(name= paste0("Washout Treatment Effect"), 
                                              low = "#f7a59b", high = "white",
                                              breaks=c(
                                                  min(panelview_data[panelview_data$four_group_ind == "early", "time_til"]),
                                                  max(panelview_data[panelview_data$four_group_ind == "early", "time_til"])),
                                              labels=c("Stronger", "Weaker"), 
                                              limits = c(
                                                  min(panelview_data[panelview_data$four_group_ind == "early", "time_til"]),
                                                  max(panelview_data[panelview_data$four_group_ind == "early", "time_til"])),
                                              na.value = "white") +
                        geom_point(data = panelview_data, aes(x = year, y = state, shape = treat_status)) +
                        scale_shape_manual(values = c(19, 1)) + 
                        labs(x = "Year", y = "State", title = "", 
                             fill = "Observation Group", 
                             shape = "Treatment Status") + 
                        ylim(levels(panelview_data$state)) + 
                        panelview_plot_tyle
                    
                } else {
                    panelview = ggplot() +
                        geom_tile(
                            data = panelview_data[panelview_data$four_group_ind == "late", ], 
                            aes(x = year, y = state, fill = time_til), position = "identity", 
                            colour="gray90", size=0.1, stat="identity") + 
                        scale_fill_continuous(name= paste0("Limited Treatment Anticipation"), 
                                              low = "white", high = "#ffc001",
                                              breaks=c(
                                                  min(panelview_data[panelview_data$four_group_ind == "late", "time_til"]),
                                                  max(panelview_data[panelview_data$four_group_ind == "late", "time_til"])),
                                              labels=c("Weaker","Stronger"), 
                                              limits = c(
                                                  min(panelview_data[panelview_data$four_group_ind == "late", "time_til"]),
                                                  max(panelview_data[panelview_data$four_group_ind == "late", "time_til"])),
                                              na.value = "white") + 
                        geom_point(data = panelview_data, aes(x = year, y = state, shape = treat_status)) +
                        scale_shape_manual(values = c(19, 1)) + 
                        labs(x = "Year", y = "State", title = "", 
                             fill = "Observation Group", 
                             shape = "Treatment Status") + 
                        ylim(levels(panelview_data$state)) + 
                        panelview_plot_tyle
                    
                }
                
            } else {
                
                if (input$washout_strength) {
                    panelview = ggplot() +
                        geom_tile(
                            data = panelview_data[panelview_data$four_group_ind == "late", ], 
                            aes(x = year, y = state, fill = time_til), position = "identity", 
                            colour="gray90", size=0.1, stat="identity") + 
                        geom_tile(
                            data = panelview_data[panelview_data$four_group_ind == "early", ], 
                            aes(x = year, y = state, fill = time_til), position = "identity", 
                            colour="gray90", size=0.1, stat="identity") + 
                        scale_fill_continuous(name= paste0("Washout Treatment Effect"), 
                                              low = "#f7a59b", high = "white",
                                              breaks=c(
                                                  min(panelview_data[panelview_data$four_group_ind == "early", "time_til"]),
                                                  max(panelview_data[panelview_data$four_group_ind == "early", "time_til"])),
                                              labels=c("Stronger", "Weaker"), 
                                              limits = c(
                                                  min(panelview_data[panelview_data$four_group_ind == "early", "time_til"]),
                                                  max(panelview_data[panelview_data$four_group_ind == "early", "time_til"])),
                                              na.value = "white") +
                        geom_point(data = panelview_data, aes(x = year, y = state, shape = treat_status)) +
                        scale_shape_manual(values = c(19, 1)) + 
                        labs(x = "Year", y = "State", title = "", 
                             fill = "Observation Group", 
                             shape = "Treatment Status") + 
                        ylim(levels(panelview_data$state)) + 
                        panelview_plot_tyle
                }
                
            }
            
        } else {
            
            panelview_data = panelview_data %>% 
                mutate(year = as.numeric(year)) %>% 
                left_join(., data_study, by = c("state", "year")) %>% 
                dplyr::select(state, year, treat_status, time_til, treat_start_year) %>% 
                mutate(state = factor(state, level = unique(state[order(treat_start_year, decreasing = T)]) ), 
                       treat_status = factor(treat_status, level = c(1, 0)))
            
            ctrl_year_ubound <- max(as.numeric(data_study$year[data_study$time_til==input$t1_1-input$t0_1]))
            ctrl_year_lbound <- min(as.numeric(data_study$year[data_study$time_til==input$t1_1-input$t0_1]))
            
            
            if (input$homogeneity) {
                panelview_data = panelview_data %>% 
                    mutate(four_group_ind = 
                               ifelse(time_til == input$t1_1-input$t0_1, "treat", 
                                      ifelse(between(year, ctrl_year_lbound, ctrl_year_ubound) & time_til == -999, "control", "invalid")), 
                           year = factor(year, level = min(year):max(year))
                    ) 
                
                
            } else {
                panelview_data = panelview_data %>% 
                    mutate(four_group_ind = 
                               ifelse(treat_start_year == input$t0_1 & year == input$t1_1, "treat", 
                                      ifelse(year == input$t1_1 & time_til==-999, "control", "invalid")), 
                           year = factor(year, level = min(year):max(year)))
            }
            
            panelview = ggplot(panelview_data,
                               aes(x = year, y = state, fill = four_group_ind), position = "identity")  +
                geom_tile(colour="gray90", size=0.1, stat="identity") +
                scale_fill_manual(values=obs_colors) +
                geom_point(data = panelview_data, aes(x = year, y = state, shape = treat_status)) +
                scale_shape_manual(values = c(19, 1)) +
                labs(x = "Year", y = "State", title = "",
                     fill = "Observation Group",
                     shape = "Treatment Status"
                ) +
                panelview_plot_tyle
        }
        panelview}, 
        height = 500, width = 730)
    
    
    output$implied = renderPlot({
        
        ################ Implied Weights ################ 
        # define the implied weights formula
        event_study_formula <- as.formula(
            paste("asmrs ~ ",
                  paste(state_covariates,collapse = " + "), " + ",
                  paste(year_covariates,collapse = " + "), " + ",
                  paste(lead_lag_covariates[-length(lead_lag_covariates)],collapse = " + ")
            )
        )
        
        all_covariates_df <- data_study_augment %>%
            dplyr::select(c(year_covariates,state_covariates,lead_lag_covariates))
        
        if (input$t1_2 >= input$t0_2) {
            drop <- c(paste0("lag",input$t1_2-input$t0_2))
            treat = paste0("lag_",input$t1_2-input$t0_2)
        } else {
            drop <- c(paste0("lead",input$t0_2-input$t1_2))
            treat = paste0("lead_",input$t0_2-input$t1_2)
        }
        all_covariates_df = all_covariates_df[,!(names(all_covariates_df) %in% drop)]
        all_covariates_df <- apply(all_covariates_df, 2, as.numeric)
        
        lmw_att_out_original = lmw::lmw(event_study_formula,
                                        data = data_study_augment,
                                        estimand = "ATT", method = "URI",
                                        treat = treat)
        weight_original = lmw_att_out_original$weights
        
        act_index = data_study_augment$X[data_study_augment$time_til == input$t1_2-input$t0_2]
        non_index = data_study_augment$X[data_study_augment$time_til != input$t1_2-input$t0_2]
        
        data_study_augment = data_study_augment %>%
            mutate(lmw_weight = weight_original, treat_group = ifelse(X %in% act_index, 1, 0))
        
        ################ Map Data ################
        # Map state names and state abbreviations
        state_names = str_to_lower(state.name)
        state_abbr = state.abb
        state_info = data.frame(cbind(state = state_abbr,
                                      state_names = state_names)) %>%
            # we do not include AK and HI
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
            rename(state_names = region) %>%
            left_join(state_center, by = "state_names") %>%
            filter(!state %in% c("AK", "HI")) %>%
            distinct() %>%
            dplyr::select(state, long, lat, c_long, c_lat, group)
        
        ################ Implied Weighted States ################
        # Join the two datasets: map data and the value (lmw weight) for each observation
        map_plot_data = state_data %>%
            left_join(data_study_augment, by="state")
        
        MapPlot = function(year_value, shade=TRUE) {
            
            if (!is.null(year_value)) {
                
                map_center_plot_data = map_plot_data %>% 
                    filter(year == year_value) %>% 
                    dplyr::select(year, state, group, c_long, c_lat, treat_status) %>% 
                    mutate(treat_status = factor(treat_status, level = c(0, 1))) %>% 
                    distinct()
                
                map_year_plot_data = map_plot_data %>%
                    filter(year == year_value) %>%
                    dplyr::select(X, year, state, long, lat, group, c_long, c_lat,
                                  lmw_weight, time_til, treat_status, treat_group) %>%
                    mutate(treat_status = factor(treat_status, level = c(0, 1)))
                
                # Plot
                if (shade) {
                    map_plot =
                        ggplot() +
                        geom_polygon(data=map_year_plot_data[map_year_plot_data$treat_group == 1, ],
                                     aes(x=long, y=lat, group=group, fill=lmw_weight), color="black", size = 0.2) +
                        scale_fill_continuous(name= paste0("Treated observation implied weight"),
                                              low = "#9dc4e4", high = "#3787c8",
                                              limits = c(
                                                  min(map_plot_data[map_plot_data$treat_group == 1, "lmw_weight"]),
                                                  max(map_plot_data[map_plot_data$treat_group == 1, "lmw_weight"])),
                                              na.value = "grey50") +
                        new_scale_fill() +
                        geom_polygon(data=map_year_plot_data[map_year_plot_data$treat_group == 0,],
                                     aes(x=long, y=lat, group=group, fill=lmw_weight), color="black", size = 0.2) +
                        scale_fill_gradient2(name= paste0("Control observation implied weight"),
                                             low = "darkgreen", mid = "white", high = "darkgreen",
                                             limits = c(
                                                 min(map_plot_data[map_plot_data$treat_group == 0, "lmw_weight"]),
                                                 max(map_plot_data[map_plot_data$treat_group == 0, "lmw_weight"])),
                                             na.value = "grey50") +
                        geom_point(data = map_center_plot_data[map_center_plot_data$treat_status == 1, ],
                                   aes(x = c_long, y = c_lat), size = 10, shape = 19) +
                        geom_text(data = map_center_plot_data[map_center_plot_data$treat_status == 1, ],
                                  aes(c_long, c_lat, label = state), color = "white", size = 5) +
                        geom_text(data = map_center_plot_data[map_center_plot_data$treat_status == 0, ],
                                  aes(c_long, c_lat, label = state), color = "black", size = 5) +
                        labs(title = paste0("Year ", year_value)) +
                        ggthemes::theme_map() + implied_plot_tyle
                } else {
                    map_plot =
                        ggplot() +
                        geom_polygon(data=map_year_plot_data[map_year_plot_data$treat_group == 1, ],
                                     aes(x=long, y=lat, group=group), fill = "#9dc4e4", color="black", size = 0.2) +
                        new_scale_fill() +
                        geom_polygon(data=map_year_plot_data[map_year_plot_data$treat_group == 0,],
                                     aes(x=long, y=lat, group=group), fill="#c6e0b4", color="black", size = 0.2) +
                        geom_point(data = map_center_plot_data[map_center_plot_data$treat_status == 1, ],
                                   aes(x = c_long, y = c_lat), size = 10, shape = 19) +
                        geom_text(data = map_center_plot_data[map_center_plot_data$treat_status == 1, ],
                                  aes(c_long, c_lat, label = state), color = "white", size = 5) +
                        geom_text(data = map_center_plot_data[map_center_plot_data$treat_status == 0, ],
                                  aes(c_long, c_lat, label = state), color = "black", size = 5) +
                        labs(title = paste0("Year ", year_value)) +
                        ggthemes::theme_map() + implied_plot_tyle
                    
                }
                
            } else {
                
                if (shade) {
                    map_plot = map_plot_data %>%
                        na.omit() %>%
                        dplyr::select(X, year, state, long, lat, group, c_long, c_lat,
                                      lmw_weight, time_til, treat_status, treat_group) %>%
                        filter(treat_group == 1) %>%
                        group_by(state) %>%
                        summarise(lmw_weight = mean(lmw_weight)) %>% ungroup() %>%
                        right_join(state_data, by = "state") %>%
                        ggplot() +
                        geom_polygon(aes(x=long, y=lat, group=group, fill=lmw_weight), color="black", size = 0.2) +
                        scale_fill_continuous(name= paste0("Observation implied weight"),
                                              low = "#9dc4e4", high = "#3787c8", na.value = "white") +
                        geom_text(data = map_center_plot_data, aes(c_long, c_lat, label = state),
                                  color = "black", size = 5) +
                        ggthemes::theme_map() + implied_plot_tyle
                } else {
                    map_plot =
                        state_data %>%
                        na.omit() %>%
                        ggplot() +
                        geom_polygon(aes(x=long, y=lat, group=group), fill="white", color="black", size = 0.2) +
                        geom_text(data = map_center_plot_data, aes(c_long, c_lat, label = state),
                                  color = "black", size = 5) +
                        ggthemes::theme_map() + implied_plot_tyle  
                }
                
            }
            
            return(map_plot)
        }
        
        ################ Implied Weighted Years ################
        
        
        ################ Render Plot ################
        if (input$implied_state) {
            if (input$implied_state_by_year) {
                MapPlot(year_value = input$year_value, shade=input$shade)
            } else {
                if (input$shade) {
                    MapPlot(year_value = NULL, shade=TRUE)
                } else {
                    MapPlot(year_value = NULL, shade=FALSE)
                } 
            }
        } 
        
    })
    
    output$balance = renderPlot({
        
        
    })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
