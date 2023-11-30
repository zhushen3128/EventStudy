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

obs_colors <- c(treat = "#9dc4e4", control = "#c6e0b4", 
                early ="#f7a59b", late = "#ffc001", invalid = "white")


ui <- fluidPage(

    # Application title
    titlePanel("Visualize Assumptions for Event Study Dataset"),
    
    
    
    # Sidebar with a slider input for number of bins 
    # The following code will create a slider with one draggable point, set by default to 0 
    sidebarLayout(
        sidebarPanel(
            h4("Estimand"), 
            selectInput("estimand", "Choose the estimand:",
                        c("ATE_{t0, >t0, t1}" = "new",
                          "ATE_{t0, inf, t1}" = "std")), 
            withMathJax(),
            helpText('Estimand can be either ATE\\(_{t_0, >t_0, t_1}\\), which is a contrast between 
                      treated versus not-yet-treated or ATE\\(_{t_0, \\infty, t_1}\\), which is a contrast 
                      between treated versus never-treated.'),
            
            sliderInput("t0",
                        "What is the t0 in the target estimand? ", 
                        value = 1980, 
                        min = 1964,
                        max = 1996), 
            withMathJax(),
            helpText('\\(t_0\\) denotes the treatment initiation time.'),
            
            sliderInput("t1",
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
                #withMathJax(),
                #helpText('\\( E \\{ Y_{it_1}(t_0) - Y_{it_1}(\\infty) \\mid X_i \\}\\)'),
                #helpText('\\( = E \\{ Y_{it_1+l}(t_0+l) - Y_{it_1+l}(\\infty) \\mid X_i \\} \\)'),
                #helpText('for all \\( t_0, t_1. \\)'), 
                
                # Horizontal line 
                # tags$hr(), 
                
                checkboxInput("washout", 'Washout of Treatment Effect.', FALSE)
                #withMathJax(),
                #helpText('\\( E \\{ Y_{it_1}(t_0) - Y_{it_1}(\\infty) \\mid X_i \\}\\)'),
                #helpText('\\( = E \\{ Y_{it_1+l}(t_0+l) - Y_{it_1+l}(\\infty) \\mid X_i \\} \\)'),
                #helpText('for all \\( t_0, t_1. \\)')
                
            ), 
            
            # Horizontal line 
            # tags$hr(), 
            
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
                checkboxInput("anticipation_strength", 'Show strength of limited treatment anticipation assumption.', FALSE), 
            ), 
            
            conditionalPanel(
                condition = "input.washout",
                checkboxInput("washout_strength", 'Show strength of washout of treatment effect assumption.', FALSE), 
            ), 
            
            
            #selectInput("strength", "Represent the strength of the anticipation and washout assumptions:",
            #            c("default" = "default", 
            #              "Strength of anticipation assumption" = "anticipation",
            #              "Strength of washout assumption" = "washout", 
            #              "Strength of both" = "both"))
        ),

        # Show a map plot of the values 
        mainPanel(plotOutput("plot2"))
    )
)


server = function(input, output){
    
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
    
    output$plot2 = renderPlot({
        # define the valid control lead lag range
        time_til_ubound <- input$t1-input$t0
        time_til_lbound <- input$t1-max(data_study$year)
        
        # generate visualization of the panel data 
        panelview_data = data.frame(cbind(rep(unique(data_study$state), 
                                              each = length(unique(data_study$year))), 
                                          rep(unique(data_study$year), 
                                              length(unique(data_study$state)))))
        colnames(panelview_data) = c("state", "year")
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
                           ifelse(time_til == input$t1-input$t0, "treat", 
                                  ifelse(between(time_til, time_til_lbound, time_til_ubound) 
                                         | time_til==-999, "control", "invalid")))
            if (input$anticipation) {
                if (input$washout) {
                    panelview_data = panelview_data %>% 
                        mutate(four_group_ind = ifelse(time_til == input$t1-input$t0, "treat", 
                                                       ifelse(time_til > max(time_til_lbound:time_til_ubound), "control", 
                                                              ifelse(between(time_til, time_til_lbound, time_til_ubound) 
                                                                     | time_til==-999, "control", "control"))))
                } else {
                    panelview_data = panelview_data %>% 
                        mutate(four_group_ind = ifelse(time_til == input$t1-input$t0, "treat", 
                                                       ifelse(time_til > max(time_til_lbound:time_til_ubound), "invalid", 
                                                              ifelse(between(time_til, time_til_lbound, time_til_ubound) 
                                                                     | time_til==-999, "control", "control"))))
                }
            } else {
                if (input$washout) {
                    panelview_data = panelview_data %>% 
                        mutate(four_group_ind = ifelse(time_til == input$t1-input$t0, "treat", 
                                                       ifelse(time_til > max(time_til_lbound:time_til_ubound), "control", 
                                                              ifelse(between(time_til, time_til_lbound, time_til_ubound) 
                                                                     | time_til==-999, "control", "invalid"))))
                } else {
                    panelview_data = panelview_data %>% 
                        mutate(four_group_ind = ifelse(time_til == input$t1-input$t0, "treat", 
                                                       ifelse(time_til > max(time_til_lbound:time_til_ubound), "invalid", 
                                                              ifelse(between(time_til, time_til_lbound, time_til_ubound) 
                                                                     | time_til==-999, "control", "invalid"))))
                }
            }
            
            
        } else {
            panelview_data = panelview_data %>% 
                mutate(
                    four_group_ind = 
                        ifelse(treat_start_year == input$t0 & year == input$t1, "treat", 
                               ifelse(year == input$t1 & (between(time_til, time_til_lbound, time_til_ubound) 
                                                          | time_til==-999), "control", "invalid")))
        }
        
        
        
        # panelview_data = panelview_data %>% 
        #     mutate(year = as.numeric(year)) %>% 
        #     left_join(., data_study, by = c("state", "year")) %>% 
        #     dplyr::select(state, year, treat_status, time_til, treat_start_year) %>% 
        #     mutate(
        #         four_group_ind = ifelse(time_til == input$t1-input$t0, "treat", 
        #                                 ifelse(time_til > max(time_til_lbound:time_til_ubound), "early", 
        #                                        ifelse(between(time_til, time_til_lbound, time_til_ubound) 
        #                                               | time_til==-999, "control", "late"))), 
        #         state = factor(state, level = unique(state[order(treat_start_year, decreasing = T)]) ), 
        #         year = factor(year, level = min(year):max(year)), 
        #         treat_status = factor(treat_status, level = c(1, 0)))
        
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
                ) + theme_bw() + 
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
        } else {
            
            panelview_data = panelview_data %>%
                mutate(
                    four_group_ind = ifelse(time_til == input$t1-input$t0, "treat",
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
                ) + theme_bw() + 
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
        } 
        
        
        # Strength of the assumptions 
        panelview_data = panelview_data %>%
            mutate(
                four_group_ind = ifelse(time_til == input$t1-input$t0, "treat",
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
                          legend.text = element_text(margin = margin(r = 10, unit = "pt"), size = 8),
                          plot.title = element_text(size=15, hjust = 0.5, face="bold",
                                                    margin = margin(8, 0, 8, 0)))
            
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
                          legend.text = element_text(margin = margin(r = 10, unit = "pt"), size = 8),
                          plot.title = element_text(size=15, hjust = 0.5, face="bold",
                                                    margin = margin(8, 0, 8, 0)))
                
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
                          legend.text = element_text(margin = margin(r = 10, unit = "pt"), size = 8),
                          plot.title = element_text(size=15, hjust = 0.5, face="bold",
                                                    margin = margin(8, 0, 8, 0)))
            }
            
        }
    panelview}, height = 500, width = 730)
    
}

# Run the application 
shinyApp(ui = ui, server = server)
