#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#    http://shiny.rstudio.com/

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
library(lfe)
library(sbw)
library(parallel)

library(shiny)
library(sf)
library(leaflet)
library(shinyTree)

#######################################################################
##    Global 
#######################################################################

obs_colors <- c(treat = "#9dc4e4", control = "#c6e0b4", 
                early ="#f7a59b", late = "#ffc001", invalid = "white")

ess_obs_colors <- c(ideal_treat = "#3293e3", ideal_control = "#69b334", 
                    invariance_treat = "#c5dbed", invariance_control = "#ddedd1", 
                    prior_treat_of_treat = "#ffc001", prior_control_of_treat = "#fae9b6", 
                    post_treat_of_treat = "#f56c5b", post_control_of_treat = "#fadbd7", 
                    treat_of_invalid = "#e1e1e1", control_of_invalid = "#fafafa")

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
dat <- read.csv("../bacon_example.csv")
data_study <- dat %>% 
    # for each state, contruct the baseline covariates 
    group_by(stfips) %>%
    dplyr::summarise(
        pcinc_baseline = pcinc[which.min(year)], 
        asmrh_baseline = asmrh[which.min(year)], 
        cases_baseline = cases[which.min(year)], 
        # construct column for treatment start year 
        treat_start_year = ifelse(sum(treat) == 0, 9999, year[time_to_treat == 0])
    ) %>% ungroup() %>% 
    right_join(dat, by = "stfips") %>% 
    dplyr::rename(time_til = time_to_treat) %>%  
    dplyr::mutate(# replace with the state names 
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
##    Helper functions 
#######################################################################

ESS <- function(w) {
    num <- sum(abs(w))^2
    den <- sum(w^2)
    ess <- num/den
    return(ess)
}

PlotESS = function(estimand, t0, t1) {
    target_name = ifelse(t1 >= t0, paste0("lag_", t1-t0, sep=""), paste0("lead_", t0-t1, sep=""))
    
    if (estimand == "std") {
        # Define different groups of observations 
        treat_set1 = data_study_augment %>% filter(year == t1 & treat_start_year == t0)
        control_set1 = data_study_augment %>% filter(year == t1 & treat_start_year == 9999) 
        invalid_set1 = data_study_augment %>% filter(year == t1 & !treat_start_year %in% c(t0, 9999))
        
        treat_set2 = data_study_augment %>% filter(time_til == t1-t0) 
        control_set2 = data_study_augment %>% filter(treat_start_year == 9999) 
        invalid_set2 = data_study_augment %>% filter(time_til != t1-t0 & treat_start_year != 9999) 
        
        treat_ideal_idx = treat_set1$X
        control_ideal_idx = control_set1$X
        treat_inv_idx = setdiff(treat_set2$X, treat_set1$X)
        control_inv_idx = setdiff(control_set2$X, control_set1$X)
        
        prior_treat_treat_set3 = data_study_augment %>% 
            filter(state %in% treat_set2$state) %>% 
            filter(treat_status == 1 & time_til < t1-t0)
        
        prior_control_treat_set3 = data_study_augment %>% 
            filter(state %in% treat_set2$state) %>% 
            filter(treat_status == 0 & time_til < t1-t0)
        
        post_treat_treat_set3 = data_study_augment %>% 
            filter(state %in% treat_set2$state) %>% 
            filter(treat_status == 1 & time_til > t1-t0)
        
        post_control_treat_set3 = data_study_augment %>% 
            filter(state %in% treat_set2$state) %>% 
            filter(treat_status == 0 & time_til > t1-t0)
        
        treat_invalid_set3 = data_study_augment %>% 
            filter(!state %in% c(treat_set2$state, control_set2$state)) %>% 
            filter(treat_status == 1)
        
        control_invalid_set3 = data_study_augment %>% 
            filter(!state %in% c(treat_set2$state, control_set2$state)) %>% 
            filter(treat_status == 0)
        
        # Create one column indicating the group for each observation 
        data_study_augment = data_study_augment %>% 
            mutate(group_ind = case_when(X %in% treat_ideal_idx ~ "ideal_treat",
                                         X %in% control_ideal_idx ~ "ideal_control",
                                         X %in% treat_inv_idx ~ "invariance_treat", 
                                         X %in% control_inv_idx ~ "invariance_control", 
                                         X %in% prior_treat_treat_set3$X ~ "prior_treat_of_treat",
                                         X %in% prior_control_treat_set3$X ~ "prior_control_of_treat",
                                         X %in% post_treat_treat_set3$X ~ "post_treat_of_treat", 
                                         X %in% post_control_treat_set3$X ~ "post_control_of_treat", 
                                         X %in% treat_invalid_set3$X ~ "treat_of_invalid", 
                                         X %in% control_invalid_set3$X ~ "control_of_invalid")) %>% 
            mutate(group_ind = factor(group_ind, levels = c("ideal_treat", "ideal_control", 
                                                            "invariance_treat", "invariance_control",
                                                            "prior_treat_of_treat", "prior_control_of_treat", 
                                                            "post_treat_of_treat", "post_control_of_treat", 
                                                            "treat_of_invalid", "control_of_invalid")))
        
        panelview_data = data.frame(cbind(rep(unique(data_study$state), 
                                              each = length(unique(data_study$year))), 
                                          rep(unique(data_study$year), 
                                              length(unique(data_study$state)))))
        colnames(panelview_data) = c("state", "year")
        
        panelview_data = panelview_data %>% 
            mutate(year = as.numeric(year)) %>% 
            left_join(., data_study_augment, by = c("state", "year")) %>% 
            dplyr::select(state, year, treat_status, time_til, treat_start_year, group_ind) %>% 
            mutate(
                state = factor(state, level = unique(state[order(treat_start_year, decreasing = T)]) ), 
                year = factor(year, level = min(year):max(year)), 
                treat_status = factor(treat_status, level = c(1, 0)))
        
        p <- ggplot(panelview_data, aes(x = year, y = state), position = "identity")  +
            geom_tile(aes(fill = group_ind), colour="gray90", size=0.1, stat="identity") + 
            geom_point(aes(x = year, y = state, shape = treat_status)) + 
            scale_fill_manual(values = ess_obs_colors) + 
            scale_shape_manual(values = c(19, 1)) + 
            labs(x = "Year", y = "State", title = "", 
                 fill = "Observation Group", 
                 shape = "Treatment Status") + panelview_plot_tyle
        p
    }
}

GetESS = function(estimand, t0, t1, method = c("twfe", "balance1", "balance2", "balance3")) {
    
    target_name = ifelse(t1 >= t0, paste0("lag_", t1-t0, sep=""), paste0("lead_", t0-t1, sep=""))
    
    if (estimand == "std") {
        # Define different groups of observations 
        treat_set1 = data_study_augment %>% filter(year == t1 & treat_start_year == t0)
        control_set1 = data_study_augment %>% filter(year == t1 & treat_start_year == 9999) 
        invalid_set1 = data_study_augment %>% filter(year == t1 & !treat_start_year %in% c(t0, 9999))
        
        treat_set2 = data_study_augment %>% filter(time_til == t1-t0) 
        control_set2 = data_study_augment %>% filter(treat_start_year == 9999) 
        invalid_set2 = data_study_augment %>% filter(time_til != t1-t0 & treat_start_year != 9999) 
        
        treat_ideal_idx = treat_set1$X
        control_ideal_idx = control_set1$X
        treat_inv_idx = setdiff(treat_set2$X, treat_set1$X)
        control_inv_idx = setdiff(control_set2$X, control_set1$X)
        
        prior_treat_treat_set3 = data_study_augment %>% 
            filter(state %in% treat_set2$state) %>% 
            filter(treat_status == 1 & time_til < t1-t0)
        
        prior_control_treat_set3 = data_study_augment %>% 
            filter(state %in% treat_set2$state) %>% 
            filter(treat_status == 0 & time_til < t1-t0)
        
        post_treat_treat_set3 = data_study_augment %>% 
            filter(state %in% treat_set2$state) %>% 
            filter(treat_status == 1 & time_til > t1-t0)
        
        post_control_treat_set3 = data_study_augment %>% 
            filter(state %in% treat_set2$state) %>% 
            filter(treat_status == 0 & time_til > t1-t0)
        
        treat_invalid_set3 = data_study_augment %>% 
            filter(!state %in% c(treat_set2$state, control_set2$state)) %>% 
            filter(treat_status == 1)
        
        control_invalid_set3 = data_study_augment %>% 
            filter(!state %in% c(treat_set2$state, control_set2$state)) %>% 
            filter(treat_status == 0)
        
        # Create one column indicating the group for each observation 
        data_study_augment = data_study_augment %>% 
            mutate(group_ind = case_when(X %in% treat_ideal_idx ~ "ideal_treat",
                                         X %in% control_ideal_idx ~ "ideal_control",
                                         X %in% treat_inv_idx ~ "invariance_treat", 
                                         X %in% control_inv_idx ~ "invariance_control", 
                                         X %in% prior_treat_treat_set3$X ~ "prior_treat_of_treat",
                                         X %in% prior_control_treat_set3$X ~ "prior_control_of_treat",
                                         X %in% post_treat_treat_set3$X ~ "post_treat_of_treat", 
                                         X %in% post_control_treat_set3$X ~ "post_control_of_treat", 
                                         X %in% treat_invalid_set3$X ~ "treat_of_invalid", 
                                         X %in% control_invalid_set3$X ~ "control_of_invalid")) %>% 
            mutate(group_ind = factor(group_ind, levels = c("ideal_treat", "ideal_control", 
                                                            "invariance_treat", "invariance_control",
                                                            "prior_treat_of_treat", "prior_control_of_treat", 
                                                            "post_treat_of_treat", "post_control_of_treat", 
                                                            "treat_of_invalid", "control_of_invalid")))
        
        if(method == "twfe") {
            # Now compute the TWFE weights 
            ## Define the formula 
            event_study_formula2 <- as.formula(
                paste("asmrs ~ ",
                      paste(state_covariates,collapse = " + "), " + ",
                      paste(year_covariates,collapse = " + "), " + ",
                      paste(lead_lag_covariates[-length(lead_lag_covariates)],collapse = " + ")
                )
            )
            
            ## TWFE regression uses all observations 
            lmw_att_out2_original = lmw::lmw(event_study_formula2, 
                                             data = data_study_augment, 
                                             estimand = "ATT", method = "URI", 
                                             treat = target_name)
            weight_original = lmw_att_out2_original$weights
            data_study_augment = data_study_augment %>% mutate(twfe_lmw_weight = weight_original)
            
        } else if (method == "balance1") {
            # Now compute the Balancing weights 
            ### Define moment covariates
            bal = list()
            bal$bal_cov = baseline_covariates
            ### Set tolerances
            bal$bal_tol = 0.02
            bal$bal_std = "group"
            bal$bal_alg = FALSE
            
            ## Setting 1: ideal experiment
            data_study_augment_subset = data_study_augment %>% 
                filter(year == t1) %>% 
                filter(treat_start_year %in% c(t0, 9999)) %>% 
                mutate(treat_ind = ifelse(treat_start_year == t0, 1, 0)) %>% 
                dplyr::select(X, treat_ind, any_of(baseline_covariates), asmrs)
            df = as.data.frame(cbind(data_study_augment_subset$X, 
                                     data_study_augment_subset$treat_ind, 
                                     data_study_augment_subset[, baseline_covariates], 
                                     data_study_augment_subset$asmrs))
            names(df) = c("X", "t", baseline_covariates, "Y")
            # Solve for the Average Treatment Effect on the Treated, ATT (default)
            sbwatt_object_1 = tryCatch(sbw(dat = df, ind = "t", out = "Y", bal = bal), error=function(e) {NULL})
            if (!is.null(sbwatt_object_1)) {
                sbwatt_weight_set1 = sbwatt_object_1$dat_weights %>% dplyr::select(X, sbw_weights) 
                data_study_augment = data_study_augment %>% left_join(sbwatt_weight_set1, by = "X")
                data_study_augment$sbw_weights = replace(data_study_augment$sbw_weights, 
                                                         is.na(data_study_augment$sbw_weights), 0)
            } 
        } else if (method == "balance2") {
            # Now compute the Balancing weights 
            ### Define moment covariates
            bal = list()
            bal$bal_cov = baseline_covariates
            ### Set tolerances
            bal$bal_tol = 0.02
            bal$bal_std = "group"
            bal$bal_alg = FALSE
            
            ## Setting 2: ideal experiment + invariance to time shifts 
            data_study_augment_subset = data_study_augment %>% 
                filter(time_til == t1-t0 | treat_start_year == 9999) %>% 
                mutate(treat_ind = ifelse(time_til == t1-t0, 1, 0)) %>% 
                dplyr::select(X, treat_ind, any_of(baseline_covariates), asmrs)
            df = as.data.frame(cbind(data_study_augment_subset$X, 
                                     data_study_augment_subset$treat_ind, 
                                     data_study_augment_subset[, baseline_covariates], 
                                     data_study_augment_subset$asmrs))
            names(df) = c("X", "t", baseline_covariates, "Y")
            # Solve for the Average Treatment Effect on the Treated, ATT (default)
            sbwatt_object_2 = tryCatch(sbw(dat = df, ind = "t", out = "Y", bal = bal), error=function(e) {NULL})
            if (!is.null(sbwatt_object_2)) {
                sbwatt_weight_set2 = sbwatt_object_2$dat_weights %>% dplyr::select(X, sbw_weights)  
                data_study_augment = data_study_augment %>% left_join(sbwatt_weight_set2, by = "X")
                data_study_augment$sbw_weights = replace(data_study_augment$sbw_weights, 
                                                         is.na(data_study_augment$sbw_weights), 0)
            } 
        } else if (method == "balance3") {
            # Now compute the Balancing weights 
            ### Define moment covariates
            bal = list()
            bal$bal_cov = baseline_covariates
            ### Set tolerances
            bal$bal_tol = 0.02
            bal$bal_std = "group"
            bal$bal_alg = FALSE
            
            ## Setting 3: all observations 
            data_study_augment_subset = data_study_augment %>% 
                mutate(treat_ind = ifelse(time_til == t1-t0, 1, 0)) %>% 
                dplyr::select(X, treat_ind, any_of(baseline_covariates), asmrs)
            df = as.data.frame(cbind(data_study_augment_subset$X, 
                                     data_study_augment_subset$treat_ind, 
                                     data_study_augment_subset[, baseline_covariates], 
                                     data_study_augment_subset$asmrs))
            names(df) = c("X", "t", baseline_covariates, "Y")
            sbwatt_object_3 = tryCatch(sbw(dat = df, ind = "t", out = "Y", bal = bal), error=function(e) {NULL})
            if (!is.null(sbwatt_object_3)) {
                sbwatt_weight_set3 = sbwatt_object_3$dat_weights %>% 
                    dplyr::select(X, sbw_weights)
                data_study_augment = data_study_augment %>% 
                    left_join(sbwatt_weight_set3, by = "X")
                data_study_augment$sbw_weights = replace(data_study_augment$sbw_weights, 
                                                         is.na(data_study_augment$sbw_weights), 0)
            } 
        }
        
        if (method == "twfe") {
            data_study_augment_ess = data_study_augment %>% 
                group_by(group_ind) %>% 
                summarise(ess = ESS(twfe_lmw_weight), 
                          ss = n(), 
                          ratio_ess_to_ss = ess / ss) %>% ungroup() %>% 
                mutate(percent_of_each = ratio_ess_to_ss / sum(ratio_ess_to_ss, na.rm = T))
            data_study_augment_ess =  rbind(data_study_augment_ess, 
                                            data.frame(group_ind='total', t(colSums(data_study_augment_ess[, -1]))))
            data_study_augment_ess[nrow(data_study_augment_ess), "ess"] = NA
        } else {
            data_study_augment_ess = data_study_augment %>% 
                group_by(group_ind) %>% 
                summarise(ess = ESS(sbw_weights), 
                          ss = n(), 
                          ratio_ess_to_ss = ess / ss) %>% ungroup() %>% 
                mutate(percent_of_each = ratio_ess_to_ss / sum(ratio_ess_to_ss, na.rm = T))
            data_study_augment_ess =  rbind(data_study_augment_ess, 
                                            data.frame(group_ind='total', t(colSums(data_study_augment_ess[, -1], na.rm = T))))
            data_study_augment_ess[nrow(data_study_augment_ess), "ess"] = NA
        }
        rownames(data_study_augment_ess) = NULL 
        colnames(data_study_augment_ess) = c("Group of Observation", 
                                             "Effective Sample Size (ESS)", 
                                             "Actual Sample Size (N)", 
                                             "ESS/N", 
                                             "Proportion of Information Borrowing")
        data_study_augment_ess
    }
}

PlotInfluence = function(inf_all, estimand, t0, t1, metric) {
    
    target_name = ifelse(t1 >= t0, paste0("lag_", t1-t0, sep=""), paste0("lead_", t0-t1, sep=""))
    
    if (estimand == "std") {
        # Define different groups of observations 
        treat_set1 = data_study_augment %>% filter(year == t1 & treat_start_year == t0)
        control_set1 = data_study_augment %>% filter(year == t1 & treat_start_year == 9999) 
        invalid_set1 = data_study_augment %>% filter(year == t1 & !treat_start_year %in% c(t0, 9999))
        
        treat_set2 = data_study_augment %>% filter(time_til == t1-t0) 
        control_set2 = data_study_augment %>% filter(treat_start_year == 9999) 
        invalid_set2 = data_study_augment %>% filter(time_til != t1-t0 & treat_start_year != 9999) 
        
        treat_ideal_idx = treat_set1$X
        control_ideal_idx = control_set1$X
        treat_inv_idx = setdiff(treat_set2$X, treat_set1$X)
        control_inv_idx = setdiff(control_set2$X, control_set1$X)
        
        prior_treat_treat_set3 = data_study_augment %>% 
            filter(state %in% treat_set2$state) %>% 
            filter(treat_status == 1 & time_til < t1-t0)
        
        prior_control_treat_set3 = data_study_augment %>% 
            filter(state %in% treat_set2$state) %>% 
            filter(treat_status == 0 & time_til < t1-t0)
        
        post_treat_treat_set3 = data_study_augment %>% 
            filter(state %in% treat_set2$state) %>% 
            filter(treat_status == 1 & time_til > t1-t0)
        
        post_control_treat_set3 = data_study_augment %>% 
            filter(state %in% treat_set2$state) %>% 
            filter(treat_status == 0 & time_til > t1-t0)
        
        treat_invalid_set3 = data_study_augment %>% 
            filter(!state %in% c(treat_set2$state, control_set2$state)) %>% 
            filter(treat_status == 1)
        
        control_invalid_set3 = data_study_augment %>% 
            filter(!state %in% c(treat_set2$state, control_set2$state)) %>% 
            filter(treat_status == 0)
        
        # Create one column indicating the group for each observation 
        data_study_augment = data_study_augment %>% 
            mutate(group_ind = case_when(X %in% treat_ideal_idx ~ "ideal_treat",
                                         X %in% control_ideal_idx ~ "ideal_control",
                                         X %in% treat_inv_idx ~ "invariance_treat", 
                                         X %in% control_inv_idx ~ "invariance_control", 
                                         X %in% prior_treat_treat_set3$X ~ "prior_treat_of_treat",
                                         X %in% prior_control_treat_set3$X ~ "prior_control_of_treat",
                                         X %in% post_treat_treat_set3$X ~ "post_treat_of_treat", 
                                         X %in% post_control_treat_set3$X ~ "post_control_of_treat", 
                                         X %in% treat_invalid_set3$X ~ "treat_of_invalid", 
                                         X %in% control_invalid_set3$X ~ "control_of_invalid")) %>% 
            mutate(group_ind = factor(group_ind, levels = c("ideal_treat", "ideal_control", 
                                                            "invariance_treat", "invariance_control",
                                                            "prior_treat_of_treat", "prior_control_of_treat", 
                                                            "post_treat_of_treat", "post_control_of_treat", 
                                                            "treat_of_invalid", "control_of_invalid")))
        
        if (metric == "sic") {
            show.r <- sort(order(abs(inf_all$inf), decreasing = TRUE)[1:10])
        } else if (metric == "sic_scaled") {
            show.r <- sort(order(abs(inf_all$inf_scaled), decreasing = TRUE)[1:10])
        } else if (metric == "est_change") {
            show.r <- sort(order(abs(inf_all$est_change), decreasing = TRUE)[1:10])
        }
        labels1 <- rep("",nrow(data_study_augment))
        labels1[show.r] <- as.character(show.r)
    
        plot_data = data_study_augment %>% 
            mutate(sic = inf_all$inf, 
                   sic_scaled = inf_all$inf_scaled, 
                   est_change = inf_all$est_change) %>% 
            rename(female_suicide = asmrs) %>% 
            dplyr::select(X, state, year, female_suicide, treat_start_year, 
                          group_ind, sic, sic_scaled, est_change)
        
        if (metric == "sic") {
            plt = ggplot(plot_data, aes(x = X, xend = X, y = 0, yend = sic, 
                                            state = state, year = year, treat_start_year = treat_start_year, 
                                            female_suicide = female_suicide, 
                                            group = group_ind, color = group_ind)) +
                    geom_segment(linetype = "solid", linewidth = 0.3) + 
                    labs(x = "Observation Index", y = "SIC") +
                    geom_text(aes(x = X, y = sic, label = labels1), size = 3, vjust = -1) +
                    scale_color_manual(values = ess_obs_colors) +
                    theme_bw()
            plt
            
        } else if (metric == "sic_scaled") {
            plt = ggplot(plot_data, aes(x = X, xend = X, y = 0, yend = sic_scaled, 
                                        state = state, year = year, treat_start_year = treat_start_year, 
                                        female_suicide = female_suicide, 
                                        group = group_ind, color = group_ind)) +
                geom_segment(linetype = "solid", linewidth = 0.3) + 
                labs(x = "Observation Index", y = "SIC Scaled") +
                geom_text(aes(x = X, y = sic_scaled, label = labels1), size = 3, vjust = -1) +
                scale_color_manual(values = ess_obs_colors) +
                theme_bw()
            plt
            
        } else {
            plt = ggplot(plot_data, aes(x = X, xend = X, y = 0, yend = est_change, 
                                        state = state, year = year, treat_start_year = treat_start_year, 
                                        female_suicide = female_suicide, 
                                        group = group_ind, color = group_ind)) +
                geom_segment(linetype = "solid", linewidth = 0.3) + 
                labs(x = "Observation Index", y = "Change in Point Estimate") +
                geom_text(aes(x = X, y = est_change, label = labels1), size = 3, vjust = -1) +
                scale_color_manual(values = ess_obs_colors) +
                theme_bw()
            plt
        }
    }
}

GetInfluence = function(estimand, t0, t1) {
    target_name = ifelse(t1 >= t0, paste0("lag_", t1-t0, sep=""), paste0("lead_", t0-t1, sep=""))
    event_study_formula2 = as.formula(
        paste("asmrs ~ ",
              paste(state_covariates,collapse = " + "), " + ",
              paste(year_covariates,collapse = " + "), " + ",
              paste(lead_lag_covariates[-length(lead_lag_covariates)],collapse = " + ")
        )
    )
    
    if (estimand == "std") {
        lmw_att_out2_original = lmw::lmw(event_study_formula2, 
                                         data = data_study_augment, 
                                         estimand = "ATT", method = "URI", 
                                         treat = target_name)
        weight_original = lmw_att_out2_original$weights
        act_index = which(data_study_augment$time_til==t1-t0)
        non_index = which(data_study_augment$time_til!=t1-t0)
        w_act = weight_original[act_index]
        w_non = weight_original[non_index]
        y_act = data_study_augment$asmrs[act_index]
        y_non = data_study_augment$asmrs[non_index]
        T_F = sum(y_act*w_act) - sum(y_non*w_non)
        
        ComputeInfluence_i = function(idx, t0, t1, data, target_name, formula, estimand = "ATT", method = "URI") {
            data_noi = data %>% filter(X != idx) 
            act_index_noi = data_noi$X[which(data_noi$time_til==t1-t0)]
            non_index_noi = data_noi$X[which(data_noi$time_til!=t1-t0)]
            
            lmw_att_out2_noi = lmw::lmw(formula, data = data_noi, estimand, method, treat = target_name)
            w_noi_act = lmw_att_out2_noi$weights[which(data_noi$X %in% act_index_noi)]
            w_noi_non = lmw_att_out2_noi$weights[which(data_noi$X %in% non_index_noi)]
            y_noi_act = data_noi$asmrs[which(data_noi$X %in% act_index_noi)]
            y_noi_non = data_noi$asmrs[which(data_noi$X %in% non_index_noi)]
            T_F_noi = sum(y_noi_act*w_noi_act) - sum(y_noi_non*w_noi_non)
            inf_i = nrow(data_noi)*(T_F_noi-T_F)
            est_change_i = T_F_noi-T_F
            return(c(inf_i, est_change_i))
        }
        
        inf_all = parallel::mcmapply(ComputeInfluence_i, data_study_augment$X,  
                                     MoreArgs = list(t0 = t0, t1 = t1, 
                                                     data = data_study_augment, 
                                                     target_name = target_name, formula = event_study_formula2), 
                                     mc.cores = 5)
        inf_all = data.frame(t(inf_all))
        colnames(inf_all) = c("inf", "est_change")
        inf_all$inf_scaled = abs(inf_all$inf)/max(abs(inf_all$inf))
    }
    return(inf_all)
}

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
                                
                                checkboxInput("invariance", 'Time Shift Invariance. ', FALSE), 
                                #withMathJax(),
                                #helpText('\\( E \\{ Y_{it_1}(t_0) - Y_{it_1}(\\infty) \\mid X_i \\}\\)'),
                                #helpText('\\( = E \\{ Y_{it_1+l}(t_0+l) - Y_{it_1+l}(\\infty) \\mid X_i \\} \\)'),
                                #helpText('for all \\( t_0, t_1. \\)'), 
                                
                                
                                conditionalPanel(
                                    condition = "input.invariance",
                                    checkboxInput("anticipation", 'Limited Treatment Anticipation.', FALSE), 
                                    checkboxInput("washout", 'Washout of Treatment Effect.', FALSE)
                                ), 
                                
                                conditionalPanel(
                                    condition = "input.invariance",
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
                                )
                            ), 
                            mainPanel(plotOutput("implied"))
                        )
               ), 
               
               tabPanel("Estimation",
                        
                        sidebarLayout(
                            sidebarPanel(
                                h4("Estimand"), 
                                selectInput("estimand_3", "Choose the estimand:",
                                            c("ATE_{t0, >t0, t1}" = "new3",
                                              "ATE_{t0, inf, t1}" = "std3")), 
                                withMathJax(),
                                helpText('Estimand can be ATE\\(_{t_0, >t_0, t_1}\\), which is the average causal 
                      effect the target population observed at time \\(t_1\\) from adopting the treatment 
                      for the first time at time \\(t_0\\) to adopting the treatment sometime after time \\(t_0\\); '), 
                                helpText('A special case of the estimand is ATE\\(_{t_0, \\infty, t_1}\\), which is is the average causal 
                      effect in the target population observed at time \\(t_1\\) from adopting the treatment for the first 
                     time at time \\(t_0\\) to never adopting the treatment.'), 
                                
                                sliderInput("t0_3",
                                            "What is the t0 in the target estimand? ", 
                                            value = 1980, 
                                            min = 1964,
                                            max = 1996), 
                                withMathJax(),
                                helpText('\\(t_0\\) denotes the treatment initiation time.'),
                                
                                sliderInput("t1_3",
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
                                
                                checkboxInput("invariance_3", 'Time Shift Invariance. ', FALSE), 
                                
                                conditionalPanel(
                                    condition = "input.invariance_3",
                                    checkboxInput("anticipation_3", 'Limited Treatment Anticipation.', FALSE), 
                                    
                                    checkboxInput("washout_3", 'Washout of Treatment Effect.', FALSE)
                                ), 
                                
                                conditionalPanel(
                                    condition = "input.invariance_3",
                                    helpText('Represent groups of observation corresponding to the assumptions above.'),
                                    checkboxInput("highlight_3", 'Show observations in the four groups with different colors.', FALSE)
                                )
                            ),
                            mainPanel(
                                plotOutput("estimatePlot"), 
                                tableOutput("estimates"))
                        )
               ), 
               
               tabPanel("Effective Sample Size",
                        sidebarLayout(
                            sidebarPanel(
                                h4("Estimand"), 
                                selectInput("estimand_4", "Choose the estimand:",
                                            c(
                                             #"ATE_{t0, >t0, t1}" = "new3",
                                              "ATE_{t0, inf, t1}" = "std4")), 
                                withMathJax(),
                                helpText('Estimand is ATE\\(_{t_0, \\infty, t_1}\\), which is is the average causal 
                      effect in the target population observed at time \\(t_1\\) from adopting the treatment for the first 
                     time at time \\(t_0\\) to never adopting the treatment.'), 
                                sliderInput("t0_4",
                                            "What is the t0 in the target estimand? ", 
                                            value = 1980, 
                                            min = 1964,
                                            max = 1996), 
                                withMathJax(),
                                helpText('\\(t_0\\) denotes the treatment initiation time.'),
                                
                                sliderInput("t1_4",
                                            "What is the t1 in the target estimand? ", 
                                            value = 1985, 
                                            min = 1964,
                                            max = 1996), 
                                withMathJax(),
                                helpText('\\(t_1\\) denotes the time of observation.'),
                                
                                # Horizontal line 
                                tags$hr(), 
                                checkboxInput("twfe", 'Compute ESS for TWFE approach.', FALSE), 
                                checkboxInput("balance1", 'Compute ESS for balancing Approach under ideal experiment setup.', FALSE), 
                                checkboxInput("balance2", 'Compute ESS for balancing Approach under time shift invariance.', FALSE), 
                                checkboxInput("balance3", 'Compute ESS for balancing Approach using all observations.', FALSE)
                                
                            ),
                            mainPanel(
                                plotOutput("ess_plot"), 
                                conditionalPanel(
                                    condition = "input.twfe",
                                    tableOutput("twfe_ess")
                                ),
                                conditionalPanel(
                                    condition = "input.balance1",
                                    tableOutput("balance1_ess")
                                ),
                                conditionalPanel(
                                    condition = "input.balance2",
                                    tableOutput("balance2_ess")
                                ),
                                conditionalPanel(
                                    condition = "input.balance3",
                                    tableOutput("balance3_ess")
                                ))
                        )
               ), 
               
               tabPanel("Sample Influence",
                        sidebarLayout(
                            sidebarPanel(
                                 h4("Estimand"), 
                                selectInput("estimand_5", "Choose the estimand:",
                                            c(
                                                #"ATE_{t0, >t0, t1}" = "new5",
                                                "ATE_{t0, inf, t1}" = "std5")),
                                withMathJax(),
                                helpText('Estimand is ATE\\(_{t_0, \\infty, t_1}\\), which is is the average causal
                      effect in the target population observed at time \\(t_1\\) from adopting the treatment for the first
                     time at time \\(t_0\\) to never adopting the treatment.'),
                                sliderInput("t0_5",
                                            "What is the t0 in the target estimand? ",
                                            value = 1980,
                                            min = 1964,
                                            max = 1996),
                                withMathJax(),
                                helpText('\\(t_0\\) denotes the treatment initiation time.'),

                                sliderInput("t1_5",
                                            "What is the t1 in the target estimand? ",
                                            value = 1985,
                                            min = 1964,
                                            max = 1996),
                                withMathJax(),
                                helpText('\\(t_1\\) denotes the time of observation.'),

                                # Horizontal line
                                tags$hr(),
                                
                                selectInput("influence_metric", "What influence metric to show:",
                                            c("Change of point estimate due to each observation." = "est_change", 
                                              "Exact SIC of each observation." = "sic",
                                              "Scaled SIC of each observation." = "sic_scaled"))

                            ), 
                             mainPanel(
                                plotOutput("infuence_plot_panel"), 
                                plotlyOutput("infuence_plot")
                             )
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
            
            if (input$invariance) {
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
                                   aes(x = year, y = state), position = "identity")  +
                    geom_tile(aes(fill = four_group_ind), colour="gray90", size=0.1, stat="identity") + 
                    scale_fill_manual(values=obs_colors) + 
                    geom_point(aes(x = year, y = state, shape = treat_status)) + 
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
                                   aes(x = year, y = state), position = "identity")  +
                    geom_tile(aes(fill = four_group_ind), colour="gray90", size=0.1, stat="identity") + 
                    scale_fill_manual(values=obs_colors) + 
                    geom_point(aes(x = year, y = state, shape = treat_status)) + 
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
            ctrl_years <- as.numeric(data_study$year[data_study$time_til==input$t1_1-input$t0_1])
            
            if (input$invariance) {
                panelview_data = panelview_data %>% 
                    mutate(four_group_ind = 
                               ifelse(time_til == input$t1_1-input$t0_1, "treat", 
                                      ifelse(between(year, ctrl_year_lbound, ctrl_year_ubound) & 
                                                 time_til == -999 &
                                                 year %in% ctrl_years, "control", "invalid")), 
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
                               aes(x = year, y = state), position = "identity")  +
                geom_tile(aes(fill = four_group_ind), colour="gray90", size=0.1, stat="identity") +
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
    
    
    output$estimatePlot = renderPlot({
        
        # define the valid control lead lag range
        time_til_ubound <- input$t1_3-input$t0_3
        time_til_lbound <- input$t1_3-max(data_study$year)
        
        # generate visualization of the panel data 
        panelview_data = data.frame(cbind(rep(unique(data_study$state), each = length(unique(data_study$year))), 
                                          rep(unique(data_study$year), length(unique(data_study$state)))))
        colnames(panelview_data) = c("state", "year")
        
        if(input$estimand_3 == "new3") {
            panelview_data = panelview_data %>% 
                mutate(year = as.numeric(year)) %>% 
                left_join(., data_study, by = c("state", "year")) %>% 
                dplyr::select(state, year, treat_status, time_til, treat_start_year) %>% 
                mutate(state = factor(state, level = unique(state[order(treat_start_year, decreasing = T)]) ), 
                       year = factor(year, level = min(year):max(year)), 
                       treat_status = factor(treat_status, level = c(1, 0)))
            
            if (input$invariance_3) {
                panelview_data = panelview_data %>% 
                    mutate(four_group_ind = 
                               ifelse(time_til == input$t1_3-input$t0_3, "treat", 
                                      ifelse(between(time_til, time_til_lbound, time_til_ubound) 
                                             | time_til==-999, "control", "invalid")))
                if (input$anticipation_3) {
                    if (input$washout) {
                        panelview_data = panelview_data %>% 
                            mutate(four_group_ind = ifelse(time_til == input$t1_3-input$t0_3, "treat", 
                                                           ifelse(time_til > max(time_til_lbound:time_til_ubound), "control", 
                                                                  ifelse(between(time_til, time_til_lbound, time_til_ubound) 
                                                                         | time_til==-999, "control", "control"))))
                    } else {
                        panelview_data = panelview_data %>% 
                            mutate(four_group_ind = ifelse(time_til == input$t1_3-input$t0_3, "treat", 
                                                           ifelse(time_til > max(time_til_lbound:time_til_ubound), "invalid", 
                                                                  ifelse(between(time_til, time_til_lbound, time_til_ubound) 
                                                                         | time_til==-999, "control", "control"))))
                    }
                } else {
                    if (input$washout) {
                        panelview_data = panelview_data %>% 
                            mutate(four_group_ind = ifelse(time_til == input$t1_3-input$t0_3, "treat", 
                                                           ifelse(time_til > max(time_til_lbound:time_til_ubound), "control", 
                                                                  ifelse(between(time_til, time_til_lbound, time_til_ubound) 
                                                                         | time_til==-999, "control", "invalid"))))
                    } else {
                        panelview_data = panelview_data %>% 
                            mutate(four_group_ind = ifelse(time_til == input$t1_3-input$t0_3, "treat", 
                                                           ifelse(time_til > max(time_til_lbound:time_til_ubound), "invalid", 
                                                                  ifelse(between(time_til, time_til_lbound, time_til_ubound) 
                                                                         | time_til==-999, "control", "invalid"))))
                    }
                }
            } else {
                panelview_data = panelview_data %>% 
                    mutate(
                        four_group_ind = 
                            ifelse(treat_start_year == input$t0_3 & year == input$t1_3, "treat", 
                                   ifelse(year == input$t1_3 & (between(time_til, time_til_lbound, time_til_ubound) 
                                                                | time_til==-999), "control", "invalid")))
            }
            
            if (!input$highlight_3) {
                panelview = ggplot(panelview_data, 
                                   aes(x = year, y = state), position = "identity")  +
                    geom_tile(aes(fill = four_group_ind), colour="gray90", size=0.1, stat="identity") + 
                    scale_fill_manual(values=obs_colors) + 
                    geom_point(data = panelview_data, aes(x = year, y = state, shape = treat_status)) + 
                    scale_shape_manual(values = c(19, 1)) + 
                    labs(x = "Year", y = "State", title = "", 
                         fill = "Observation Group", 
                         shape = "Treatment Status") + panelview_plot_tyle
            } else {
                panelview_data = panelview_data %>%
                    mutate(
                        four_group_ind = ifelse(time_til == input$t1_3-input$t0_3, "treat",
                                                ifelse(time_til > max(time_til_lbound:time_til_ubound), "early",
                                                       ifelse(between(time_til, time_til_lbound, time_til_ubound)
                                                              | time_til==-999, "control", "late"))))
                
                panelview = ggplot(panelview_data, 
                                   aes(x = year, y = state), position = "identity")  +
                    geom_tile(aes(fill = four_group_ind), colour="gray90", size=0.1, stat="identity") + 
                    scale_fill_manual(values=obs_colors) + 
                    geom_point(data = panelview_data, aes(x = year, y = state, shape = treat_status)) + 
                    scale_shape_manual(values = c(19, 1)) + 
                    labs(x = "Year", y = "State", title = "", 
                         fill = "Observation Group", 
                         shape = "Treatment Status") + panelview_plot_tyle
            } 
        } else {
            
            panelview_data = panelview_data %>% 
                mutate(year = as.numeric(year)) %>% 
                left_join(., data_study, by = c("state", "year")) %>% 
                dplyr::select(state, year, treat_status, time_til, treat_start_year) %>% 
                mutate(state = factor(state, level = unique(state[order(treat_start_year, decreasing = T)]) ), 
                       treat_status = factor(treat_status, level = c(1, 0)))
            
            ctrl_year_ubound <- max(as.numeric(data_study$year[data_study$time_til==input$t1_3-input$t0_3]))
            ctrl_year_lbound <- min(as.numeric(data_study$year[data_study$time_til==input$t1_3-input$t0_3]))
            ctrl_years <- as.numeric(data_study$year[data_study$time_til==input$t1_3-input$t0_3])
            
            if (input$invariance_3) {
                # panelview_data = panelview_data %>% 
                #     mutate(four_group_ind = 
                #                ifelse(time_til == input$t1_3-input$t0_3, "treat", 
                #                       ifelse(between(year, ctrl_year_lbound, ctrl_year_ubound) & 
                #                                  time_til == -999 &
                #                                  year %in% ctrl_years, "control", "invalid")), 
                #            year = factor(year, level = min(year):max(year))
                #     ) 
                
                #SZ 20231212: change to stronger version for balancing appraoch to find solution 
                panelview_data = panelview_data %>% 
                    mutate(four_group_ind = 
                               ifelse(time_til == input$t1_3-input$t0_3, "treat", 
                                      ifelse(time_til == -999, "control", "invalid")), 
                           year = factor(year, level = min(year):max(year))
                    ) 
                
            } else {
                panelview_data = panelview_data %>% 
                    mutate(four_group_ind = 
                               ifelse(treat_start_year == input$t0_3 & year == input$t1_3, "treat", 
                                      ifelse(year == input$t1_3 & time_til==-999, "control", "invalid")), 
                           year = factor(year, level = min(year):max(year)))
            }
            
            panelview = ggplot(panelview_data,
                               aes(x = year, y = state), position = "identity")  +
                geom_tile(aes(fill = four_group_ind), colour="gray90", size=0.1, stat="identity") +
                scale_fill_manual(values=obs_colors) +
                geom_point(data = panelview_data, aes(x = year, y = state, shape = treat_status)) +
                scale_shape_manual(values = c(19, 1)) +
                labs(x = "Year", y = "State", title = "",
                     fill = "Observation Group",
                     shape = "Treatment Status") + panelview_plot_tyle
        }
        panelview}, 
        height = 400, width = 630)
    
    output$estimates = renderTable({
        
        EstimateLeadLag = function(estimand, t0, t1) {
            
            # define the valid control lead lag range
            time_til_ubound <- t1-t0
            time_til_lbound <- t1-max(data_study$year)
            
            if (estimand == "std") {
                # Simple Average 
                # minimal assumption 
                data_study_treat1 = data_study_augment %>% 
                    filter(year == t1 & treat_start_year == t0)
                data_study_control1 = data_study_augment %>% 
                    filter(year == t1 & treat_start_year == 9999)
                est1 = mean(data_study_treat1$asmrs) - mean(data_study_control1$asmrs)
                
                # minimal assumption + invariance to time shifts 
                data_study_treat2 = data_study_augment %>% filter(time_til == t1-t0) 
                data_study_control2 = data_study_augment %>% 
                    filter(year %in% data_study_treat2$year & treat_start_year == 9999)
                est2 = mean(data_study_treat2$asmrs) - mean(data_study_control2$asmrs)
                
                # all observations
                data_study_treat3 = data_study_augment %>% filter(time_til == t1-t0) 
                data_study_control3 = data_study_augment %>% 
                    filter(! X %in% data_study_treat3$X)
                est3 = mean(data_study_treat3$asmrs) - mean(data_study_control3$asmrs)
                
                mean_diff = round(c(est1, est2, est3), 3)
                
                # TWFE
                event_study_formula <- as.formula(
                    paste("asmrs ~ ",
                          paste(
                              paste(paste("lead_", 1:lead_min, sep = ""), collapse = " + "),
                              paste(paste("lag_", 0:lag_max, sep = ""), collapse = " + "), sep = " + "),
                          "| year + state | 0 | state"
                    ),
                )
                # all observations
                data_study_augment_w = data_study_augment %>% 
                    mutate(w = ifelse(X %in% c(data_study_treat3$X, data_study_control3$X), 1, 0)) %>% pull(w)
                event_study_reg_w <- felm(event_study_formula, data = data_study_augment, weights = data_study_augment_w)
                coeff_w3 = coef(event_study_reg_w)
                names(coeff_w3) = c(paste("lead_", 1:lead_min, sep = ""), paste("lag_", 0:lag_max, sep = ""))
                target_name = ifelse(t1 >= t0, paste0("lag_", t1-t0, sep=""), paste0("lead_", t0-t1, sep=""))
                twfe_est = c("-", "-", round(coeff_w3[target_name], 3))
                
                # Balance 
                # define moment covariates
                bal = list()
                bal$bal_cov = baseline_covariates
                # Set tolerances
                bal$bal_tol = 0.02
                bal$bal_std = "group"
                bal$bal_alg = FALSE
                
                data_study_augment_subset = data_study_augment %>% 
                    filter(year == t1) %>% 
                    filter(treat_start_year %in% c(t0, 9999)) %>% 
                    mutate(treat_ind = ifelse(treat_start_year == t0, 1, 0)) %>% 
                    dplyr::select(X, treat_ind, any_of(baseline_covariates), asmrs)
                
                df = as.data.frame(cbind(data_study_augment_subset$X, 
                                         data_study_augment_subset$treat_ind, 
                                         data_study_augment_subset[, baseline_covariates], 
                                         data_study_augment_subset$asmrs))
                names(df) = c("idx", "t", baseline_covariates, "Y")
                sbwatt_object_1 = tryCatch(sbw(dat = df, ind = "t", out = "Y", bal = bal), error=function(e) {NULL})
                
                data_study_augment_subset = data_study_augment %>% 
                    filter(X %in% c(data_study_treat2$X, data_study_control2$X)) %>% 
                    mutate(treat_ind = ifelse(time_til == t1-t0, 1, 0)) %>% 
                    dplyr::select(X, treat_ind, any_of(baseline_covariates), asmrs)
                df = as.data.frame(cbind(data_study_augment_subset$X, 
                                         data_study_augment_subset$treat_ind, 
                                         data_study_augment_subset[, baseline_covariates], 
                                         data_study_augment_subset$asmrs))
                names(df) = c("idx", "t", baseline_covariates, "Y")
                sbwatt_object_2 = tryCatch(sbw(dat = df, ind = "t", out = "Y", bal = bal), error=function(e) {NULL})
                
                data_study_augment_subset = data_study_augment %>% 
                    filter(X %in% c(data_study_treat3$X, data_study_control3$X)) %>% 
                    mutate(treat_ind = ifelse(time_til == t1-t0, 1, 0)) %>% 
                    dplyr::select(X, treat_ind, any_of(baseline_covariates), asmrs)
                df = as.data.frame(cbind(data_study_augment_subset$X, 
                                         data_study_augment_subset$treat_ind, 
                                         data_study_augment_subset[, baseline_covariates], 
                                         data_study_augment_subset$asmrs))
                names(df) = c("idx", "t", baseline_covariates, "Y")
                sbwatt_object_3 = tryCatch(sbw(dat = df, ind = "t", out = "Y", bal = bal), error=function(e) {NULL})
                
                balance_est1 = tryCatch(round(as.numeric(estimate(sbwatt_object_1)[[1]][1]), 3), error=function(e) NULL)
                balance_est2 = tryCatch(round(as.numeric(estimate(sbwatt_object_2)[[1]][1]), 3), error=function(e) NULL)
                balance_est3 = tryCatch(round(as.numeric(estimate(sbwatt_object_3)[[1]][1]), 3), error=function(e) NULL)
                
                if (is.null(balance_est1)) {
                    balance_est1 = "-"
                } 
                
                if (is.null(balance_est2)) {
                    balance_est2 = "-"
                }
                
                if (is.null(balance_est3)) {
                    balance_est3 = "-"
                }
                
                balance_est = c(balance_est1, balance_est2, balance_est3)
                
                summary = rbind(mean_diff, twfe_est, balance_est)
                rownames(summary) = c("Mean Difference", "TWFE", "Balance")
                colnames(summary) = c("Ideal Experiment", "Time Shift Invariance", "All Observations")
                summary
                
            } else if (estimand == "new") {
                
                treat_idx = data_study$X[data_study$time_til == t1-t0]
                
                control_idx = data_study$X[(between(data_study$time_til, time_til_lbound, time_til_ubound) 
                                            | data_study$time_til==-999) & data_study$time_til != t1-t0]
                
                early_idx = data_study$X[data_study$time_til > max(time_til_lbound:time_til_ubound)]
                
                late_idx = data_study$X[data_study$time_til < min(time_til_lbound:time_til_ubound) & 
                                            data_study$time_til != -999]
                
                # Simple Average 
                # minimal assumption 
                data_study_treat1 = data_study_augment %>% 
                    filter(year == t1 & treat_start_year == t0)
                data_study_control1 = data_study_augment %>% 
                    filter(year == t1 & treat_start_year > t0)
                est1 = mean(data_study_treat1$asmrs) - mean(data_study_control1$asmrs)
                
                # minimal assumption + invariance to time shifts 
                data_study_treat2 = data_study_augment %>% filter(time_til == t1-t0) 
                data_study_control2 = data_study_augment %>% filter(X %in% control_idx)
                est2 = mean(data_study_treat2$asmrs) - mean(data_study_control2$asmrs)
                
                # all observations
                data_study_treat3 = data_study_augment %>% filter(time_til == t1-t0) 
                data_study_control3 = data_study_augment %>% filter(! X %in% data_study_treat3$X)
                est3 = mean(data_study_treat3$asmrs) - mean(data_study_control3$asmrs)
                mean_diff = round(c(est1, est2, est3), 3)
                
                # TWFE
                event_study_formula <- as.formula(
                    paste("asmrs ~ ",
                          paste(
                              # paste(baseline_covariates, collapse = " + "),
                              paste(paste("lead_", 1:lead_min, sep = ""), collapse = " + "),
                              paste(paste("lag_", 0:lag_max, sep = ""), collapse = " + "), sep = " + "),
                          "| year + state | 0 | state"
                    )
                )
                
                # all observations
                data_study_augment_w = data_study_augment %>%
                    mutate(w = ifelse(X %in% c(data_study_treat3$X, data_study_control3$X), 1, 0)) %>% pull(w)
                event_study_reg_w <- felm(event_study_formula, data = data_study_augment, weights = data_study_augment_w)
                coeff_w3 = coef(event_study_reg_w)
                names(coeff_w3) = c(paste("lead_", 1:lead_min, sep = ""), paste("lag_", 0:lag_max, sep = ""))
                target_name = ifelse(t1 >= t0, paste0("lag_", t1-t0, sep=""), paste0("lead_", t0-t1, sep=""))
                twfe_est = c("-", "-", round(coeff_w3[target_name], 3))

                if (t1 >= t0) {
                    correction_lags = paste("lag_", (t1-t0+1):lag_max, sep = "")
                } else {
                    correction_lags = c(paste("lead_", (t0-t1+1):lead_min, sep = ""),
                                        paste("lag_", 0:lag_max, sep = ""))
                }

                correction = 0
                for (lag in correction_lags) {
                    prop_lag = sum(data_study_augment[, lag] == 1) / nrow(data_study_augment)
                    correction = correction + as.numeric(prop_lag * coeff_w3[lag])
                }
                coeff_w3[target_name] = coeff_w3[target_name] - correction
                twfe_est = c("-", "-", round(coeff_w3["lag_5"], 3))
                
                # Balance 
                data_study_augment_subset = data_study_augment %>% 
                    filter(year == t1) %>% 
                    filter(treat_start_year >= t0) %>% 
                    mutate(treat_ind = ifelse(treat_start_year == t0, 1, 0)) %>% 
                    dplyr::select(X, treat_ind, any_of(baseline_covariates), asmrs)
                
                df = as.data.frame(cbind(data_study_augment_subset$X, 
                                         data_study_augment_subset$treat_ind, 
                                         data_study_augment_subset[, baseline_covariates], 
                                         data_study_augment_subset$asmrs))
                names(df) = c("idx", "t", baseline_covariates, "Y")
                sbwatt_object_1 = tryCatch(sbw(dat = df, ind = "t", out = "Y", bal = bal), error = function(e) NULL)
                
                data_study_augment_subset = data_study_augment %>% 
                    filter(X %in% c(data_study_treat2$X, data_study_control2$X)) %>% 
                    mutate(treat_ind = ifelse(time_til == t1-t0, 1, 0)) %>% 
                    dplyr::select(X, treat_ind, any_of(baseline_covariates), asmrs)
                df = as.data.frame(cbind(data_study_augment_subset$X, 
                                         data_study_augment_subset$treat_ind, 
                                         data_study_augment_subset[, baseline_covariates], 
                                         data_study_augment_subset$asmrs))
                names(df) = c("idx", "t", baseline_covariates, "Y")
                sbwatt_object_2 = tryCatch(sbw(dat = df, ind = "t", out = "Y", bal = bal), error=function(e) NULL)
                
                data_study_augment_subset = data_study_augment %>% 
                    filter(X %in% c(data_study_treat3$X, data_study_control3$X)) %>% 
                    mutate(treat_ind = ifelse(time_til == t1-t0, 1, 0)) %>% 
                    dplyr::select(X, treat_ind, any_of(baseline_covariates), asmrs)
                df = as.data.frame(cbind(data_study_augment_subset$X, 
                                         data_study_augment_subset$treat_ind, 
                                         data_study_augment_subset[, baseline_covariates], 
                                         data_study_augment_subset$asmrs))
                names(df) = c("idx", "t", baseline_covariates, "Y")
                sbwatt_object_3 = tryCatch(sbw(dat = df, ind = "t", out = "Y", bal = bal), error=function(e) NULL)
                
                
                balance_est1 = tryCatch(round(as.numeric(estimate(sbwatt_object_1)[[1]][1]), 3), error=function(e) NULL)
                balance_est2 = tryCatch(round(as.numeric(estimate(sbwatt_object_2)[[1]][1]), 3), error=function(e) NULL)
                balance_est3 = tryCatch(round(as.numeric(estimate(sbwatt_object_3)[[1]][1]), 3), error=function(e) NULL)
                
                if (is.null(balance_est1)) {
                    balance_est1 = "-"
                } 
                
                if (is.null(balance_est2)) {
                    balance_est2 = "-"
                }
                
                if (is.null(balance_est3)) {
                    balance_est3 = "-"
                }
                
                balance_est = c(balance_est1, balance_est2, balance_est3)
                summary = rbind(mean_diff, twfe_est, balance_est)
                rownames(summary) = c("Mean Difference", "TWFE", "Balance")
                colnames(summary) = c("Ideal Experiment", "Time Shift Invariance", "All Observations")
                summary
            }
        }
        
        if (input$estimand_3 == "std3") {
            EstimateLeadLag(estimand = "std", t0 = input$t0_3, t1 = input$t1_3)
        } else if (input$estimand_3 == "new3") {
            EstimateLeadLag(estimand = "new", t0 = input$t0_3, t1 = input$t1_3)
        }
    }, 
    striped = TRUE,
    spacing = "l",
    digits = 3,
    rownames = TRUE, 
    width = "100%",
    caption = ""
    )
    
    output$ess_plot = renderPlot({
        if (input$estimand_4 == "std4") {
            PlotESS(estimand = "std", t0 = input$t0_4, t1 = input$t1_4)
        } 
    })
    
    output$twfe_ess = renderTable({
        if (input$estimand_4 == "std4") {
            GetESS(estimand = "std", t0 = input$t0_4, t1 = input$t1_4, method = "twfe")
        } 
    }, 
    striped = TRUE,
    spacing = "l",
    digits = 3,
    rownames = TRUE, 
    width = "100%",
    caption = ""
    )
    
    output$balance1_ess = renderTable({
        if (input$estimand_4 == "std4") {
            GetESS(estimand = "std", t0 = input$t0_4, t1 = input$t1_4, method = "balance1")
        } 
    }, 
    striped = TRUE,
    spacing = "l",
    digits = 3,
    rownames = TRUE, 
    width = "100%",
    caption = ""
    )
    
    output$balance2_ess = renderTable({
        if (input$estimand_4 == "std4") {
            GetESS(estimand = "std", t0 = input$t0_4, t1 = input$t1_4, method = "balance2")
        } 
    }, 
    striped = TRUE,
    spacing = "l",
    digits = 3,
    rownames = TRUE, 
    width = "100%",
    caption = ""
    )
    
    output$balance3_ess = renderTable({
        if (input$estimand_4 == "std4") {
            GetESS(estimand = "std", t0 = input$t0_4, t1 = input$t1_4, method = "balance3")
        } 
    }, 
    striped = TRUE,
    spacing = "l",
    digits = 3,
    rownames = TRUE, 
    width = "100%",
    caption = ""
    )
    
    inf_res = reactive({ 
        GetInfluence(estimand = "std", t0 = input$t0_5, t1 = input$t1_5) 
    })
    
    output$infuence_plot_panel = renderPlot({
        PlotESS(estimand = "std", t0 = input$t0_5, t1 = input$t1_5)
    })
    
    output$infuence_plot = renderPlotly({
        if (input$influence_metric == "sic") {
            change_plt = PlotInfluence(inf_res(), estimand = "std", input$t0_5, input$t1_5, metric = "sic")
            ggplotly(change_plt, tooltip = c("group", "state", "year", "treat_start_year","female_suicide","sic"))
        } else if (input$influence_metric == "sic_scaled") {
            change_plt = PlotInfluence(inf_res(), estimand = "std", input$t0_5, input$t1_5, metric = "sic_scaled")
            ggplotly(change_plt, tooltip = c("group", "state", "year", "treat_start_year","female_suicide","sic_scaled"))
        } else {
            change_plt = PlotInfluence(inf_res(), estimand = "std", input$t0_5, input$t1_5, metric = "est_change")
            ggplotly(change_plt, tooltip = c("group", "state", "year", "treat_start_year","female_suicide","est_change"))
        }
    })

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
                        geom_text(data = map_plot_data, aes(c_long, c_lat, label = state),
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
    
    
    
    
    
    # output$subplot = renderPlot({
    #     
    #     time_til_ubound <- input$t1_3-input$t0_3
    #     time_til_lbound <- input$t1_3-max(data_study$year)
    #     
    #     # generate visualization of the panel data 
    #     panelview_data = data.frame(cbind(rep(unique(data_study$state), each = length(unique(data_study$year))), 
    #                                       rep(unique(data_study$year), length(unique(data_study$state)))))
    #     colnames(panelview_data) = c("state", "year")
    #     
    #     if(input$estimand == "new3") {
    #         panelview_data = panelview_data %>% 
    #             mutate(year = as.numeric(year)) %>% 
    #             left_join(., data_study, by = c("state", "year")) %>% 
    #             dplyr::select(state, year, treat_status, time_til, treat_start_year) %>% 
    #             mutate(state = factor(state, level = unique(state[order(treat_start_year, decreasing = T)]) ), 
    #                    year = factor(year, level = min(year):max(year)), 
    #                    treat_status = factor(treat_status, level = c(1, 0)))
    #             
    #         if (input$invariance_3) {
    #             panelview_data = panelview_data %>% 
    #                 mutate(four_group_ind = 
    #                            ifelse(time_til == input$t1_3-input$t0_3, "treat", 
    #                                   ifelse(between(time_til, time_til_lbound, time_til_ubound) 
    #                                          | time_til==-999, "control", "invalid")))
    #             if (input$anticipation_3) {
    #                 if (input$washout_3) {
    #                     panelview_data = panelview_data %>% 
    #                         mutate(four_group_ind = ifelse(time_til == input$t1_3-input$t0_3, "treat", 
    #                                                        ifelse(time_til > max(time_til_lbound:time_til_ubound), "control", 
    #                                                               ifelse(between(time_til, time_til_lbound, time_til_ubound) 
    #                                                                      | time_til==-999, "control", "control"))))
    #                 } else {
    #                     panelview_data = panelview_data %>% 
    #                         mutate(four_group_ind = ifelse(time_til == input$t1_3-input$t0_3, "treat", 
    #                                                        ifelse(time_til > max(time_til_lbound:time_til_ubound), "invalid", 
    #                                                               ifelse(between(time_til, time_til_lbound, time_til_ubound) 
    #                                                                      | time_til==-999, "control", "control"))))
    #                 }
    #             } else {
    #                 if (input$washout_3) {
    #                     panelview_data = panelview_data %>% 
    #                         mutate(four_group_ind = ifelse(time_til == input$t1_3-input$t0_3, "treat", 
    #                                                        ifelse(time_til > max(time_til_lbound:time_til_ubound), "control", 
    #                                                               ifelse(between(time_til, time_til_lbound, time_til_ubound) 
    #                                                                      | time_til==-999, "control", "invalid"))))
    #                 } else {
    #                     panelview_data = panelview_data %>% 
    #                         mutate(four_group_ind = ifelse(time_til == input$t1_3-input$t0_3, "treat", 
    #                                                        ifelse(time_til > max(time_til_lbound:time_til_ubound), "invalid", 
    #                                                               ifelse(between(time_til, time_til_lbound, time_til_ubound) 
    #                                                                      | time_til==-999, "control", "invalid"))))
    #                 }
    #             }
    #         } else {
    #             panelview_data = panelview_data %>% 
    #                 mutate(
    #                     four_group_ind = 
    #                         ifelse(treat_start_year == input$t0_3 & year == input$t1_3, "treat", 
    #                                ifelse(year == input$t1_3 & (between(time_til, time_til_lbound, time_til_ubound) 
    #                                                             | time_til==-999), "control", "invalid")))
    #         }
    #         
    #         panelview = ggplot(panelview_data, aes(x = year, y = state, fill = four_group_ind), position = "identity")  +
    #             geom_tile(colour="gray90", size=0.1, stat="identity") + 
    #             scale_fill_manual(values=obs_colors) + 
    #             geom_point(data = panelview_data, aes(x = year, y = state, shape = treat_status)) + 
    #             scale_shape_manual(values = c(19, 1)) + 
    #             labs(x = "Year", y = "State", title = "", 
    #                  fill = "Observation Group", 
    #                  shape = "Treatment Status") + 
    #             panelview_plot_tyle
    #         
    #         } else {
    #             
    #             panelview_data = panelview_data %>% 
    #                 mutate(year = as.numeric(year)) %>% 
    #                 left_join(., data_study, by = c("state", "year")) %>% 
    #                 dplyr::select(state, year, treat_status, time_til, treat_start_year) %>% 
    #                 mutate(state = factor(state, level = unique(state[order(treat_start_year, decreasing = T)]) ), 
    #                        treat_status = factor(treat_status, level = c(1, 0)))
    #             
    #             ctrl_year_ubound <- max(as.numeric(data_study$year[data_study$time_til==input$t1_3-input$t0_3]))
    #             ctrl_year_lbound <- min(as.numeric(data_study$year[data_study$time_til==input$t1_3-input$t0_3]))
    #             ctrl_years <- as.numeric(data_study$year[data_study$time_til==input$t1_3-input$t0_3])
    #             
    #             if (input$invariance_3) {
    #                 panelview_data = panelview_data %>% 
    #                     mutate(four_group_ind = 
    #                                ifelse(time_til == input$t1_3-input$t0_3, "treat", 
    #                                       ifelse(between(year, ctrl_year_lbound, ctrl_year_ubound) & 
    #                                                  time_til == -999 &
    #                                                  year %in% ctrl_years, "control", "invalid")), 
    #                            year = factor(year, level = min(year):max(year))
    #                     ) 
    #                 
    #             } else {
    #                 panelview_data = panelview_data %>% 
    #                     mutate(four_group_ind = 
    #                                ifelse(treat_start_year == input$t0_3 & year == input$t1_3, "treat", 
    #                                       ifelse(year == input$t1_3 & time_til==-999, "control", "invalid")), 
    #                            year = factor(year, level = min(year):max(year)))
    #             }
    #             
    #             panelview = ggplot(panelview_data,
    #                                aes(x = year, y = state, fill = four_group_ind), position = "identity")  +
    #                 geom_tile(colour="gray90", size=0.1, stat="identity") +
    #                 scale_fill_manual(values=obs_colors) +
    #                 geom_point(data = panelview_data, aes(x = year, y = state, shape = treat_status)) +
    #                 scale_shape_manual(values = c(19, 1)) +
    #                 labs(x = "Year", y = "State", title = "",
    #                      fill = "Observation Group",
    #                      shape = "Treatment Status"
    #                 ) + panelview_plot_tyle
    #         }
    #     panelview})
        
    
}

# Run the application 
shinyApp(ui = ui, server = server)
