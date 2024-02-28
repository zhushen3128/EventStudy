

ESS <- function(w) {
  num <- sum(abs(w))^2
  den <- sum(w^2)
  ess <- num/den
  return(ess)
}

# GetESS function takes in the augmented data with dummy columns and the lmw_weights column - from the TWFE model specified by the user 
## arg: 
## 

GetESS = function(data, bal_covariates, method = c("twfe", "experiment", "invariance", "anticipate", "delay", "dissipate")) {
  bal_covariates_names = paste("X", 1:length(bal_covariates), sep = "")
  # Now compute the Balancing weights 
  bal = list()
  bal$bal_cov = bal_covariates_names
  ### Set tolerances
  bal$bal_std = "group"
  # Use adaptive algorithm 
  bal$bal_alg = TRUE
  bal$bal_gri = c(0.0001, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)
  
  if (method == "experiment") {
    # Setting 1: ideal experiment
    df = data %>% 
      filter(group_ind %in% c("IdealTreat", "IdealControl")) %>% 
      mutate(treat_ind = ifelse(group_ind == "IdealTreat", 1, 0)) %>% 
      mutate(treat_ind = factor(treat_ind)) %>% 
      dplyr::select(X, treat_ind, any_of(bal_covariates_names), Outcome) %>% as.data.frame()
    names(df) = c("X", "t", bal_covariates_names, "Y")
    
    # Solve for the Average Treatment Effect on the Treated, ATT (default)
    sbwatt_object_1 = tryCatch(sbw(dat = df, ind = "t", out = "Y", bal = bal), error=function(e) {NULL})
    if (!is.null(sbwatt_object_1)) {
      sbwatt_weight_set1 = sbwatt_object_1$dat_weights %>% dplyr::select(X, sbw_weights) 
      data = data %>% left_join(sbwatt_weight_set1, by = "X")
      data$sbw_weights = base::replace(data$sbw_weights, is.na(data$sbw_weights), 0)
    } 
  } else if (method == "invariance") {
    # Setting 2: ideal experiment + invariance to time shifts 
    df = data %>% 
      filter(group_ind %in% c("IdealTreat", "IdealControl", "InvTreat", "InvControl")) %>% 
      mutate(treat_ind = ifelse(group_ind %in% c("IdealTreat", "InvTreat"), 1, 0)) %>% 
      mutate(treat_ind = factor(treat_ind)) %>% 
      dplyr::select(X, treat_ind, any_of(bal_covariates_names), Outcome) %>% as.data.frame()
    names(df) = c("X", "t", bal_covariates_names, "Y")
    
    # Solve for the Average Treatment Effect on the Treated, ATT (default)
    sbwatt_object_2 = tryCatch(sbw(dat = df, ind = "t", out = "Y", bal = bal), error=function(e) {NULL})
    if (!is.null(sbwatt_object_2)) {
      sbwatt_weight_set2 = sbwatt_object_2$dat_weights %>% dplyr::select(X, sbw_weights)  
      data = data %>% left_join(sbwatt_weight_set2, by = "X")
      data$sbw_weights = base::replace(data$sbw_weights, is.na(data$sbw_weights), 0)
    } 
  } else if (method == "anticipate") {
    # Setting 3: after invoking the limited anticipation assumption 
    df = data %>% 
      filter(group_ind %in% c("IdealTreat", "IdealControl", "InvTreat", "InvControl", "LimitAnticipate")) %>% 
      mutate(treat_ind = ifelse(group_ind %in% c("IdealTreat", "InvTreat"), 1, 0)) %>% 
      mutate(treat_ind = factor(treat_ind)) %>% 
      dplyr::select(X, treat_ind, any_of(bal_covariates_names), Outcome) %>% as.data.frame()
    names(df) = c("X", "t", bal_covariates_names, "Y")
    sbwatt_object_3 = tryCatch(sbw(dat = df, ind = "t", out = "Y", bal = bal), error=function(e) {NULL})
    
    # Solve for the Average Treatment Effect on the Treated, ATT (default)
    if (!is.null(sbwatt_object_3)) {
      sbwatt_weight_set3 = sbwatt_object_3$dat_weights %>% dplyr::select(X, sbw_weights)
      data = data %>% left_join(sbwatt_weight_set3, by = "X")
      data$sbw_weights = base::replace(data$sbw_weights, is.na(data$sbw_weights), 0)
    } 
  } else if (method == "delay") {
    # Setting 4: after invoking the delayed treatment onset assumption 
    df = data %>% 
      filter(group_ind != "Dissipate") %>% 
      mutate(treat_ind = ifelse(group_ind %in% c("IdealTreat", "InvTreat"), 1, 0)) %>% 
      mutate(treat_ind = factor(treat_ind)) %>% 
      dplyr::select(X, treat_ind, any_of(bal_covariates_names), Outcome) %>% as.data.frame()
    names(df) = c("X", "t", bal_covariates_names, "Y")
    sbwatt_object_4 = tryCatch(sbw(dat = df, ind = "t", out = "Y", bal = bal), error=function(e) {NULL})
    
    # Solve for the Average Treatment Effect on the Treated, ATT (default)
    if (!is.null(sbwatt_object_4)) {
      sbwatt_weight_set4 = sbwatt_object_4$dat_weights %>% dplyr::select(X, sbw_weights)
      data = data %>% left_join(sbwatt_weight_set4, by = "X")
      data$sbw_weights = base::replace(data$sbw_weights, is.na(data$sbw_weights), 0)
    } 
  } else if (method == "dissipate") {
    # Setting 5: after invoking the treatment effect dissipation assumption 
    df = data %>% 
      mutate(treat_ind = ifelse(group_ind %in% c("IdealTreat", "InvTreat"), 1, 0)) %>% 
      dplyr::select(X, treat_ind, any_of(bal_covariates_names), Outcome) %>% as.data.frame()
    names(df) = c("X", "t", bal_covariates_names, "Y")
    sbwatt_object_5 = tryCatch(sbw(dat = df, ind = "t", out = "Y", bal = bal), error=function(e) {NULL})
    
    # Solve for the Average Treatment Effect on the Treated, ATT (default)
    if (!is.null(sbwatt_object_5)) {
      sbwatt_weight_set5 = sbwatt_object_5$dat_weights %>% dplyr::select(X, sbw_weights)
      data = data %>% left_join(sbwatt_weight_set5, by = "X")
      data$sbw_weights = base::replace(data$sbw_weights, is.na(data$sbw_weights), 0)
    } 
  }
  
  if (method == "twfe") {
    data_ess = data %>% 
      group_by(group_ind) %>% 
      summarise(
        ss = n(), 
        ess = ESS(lmw_weight),
        ratio_ess_to_ss = ess / ss) %>% ungroup() %>% 
      mutate(percent_of_each = ess / sum(ess, na.rm = T))
    data_ess =  rbind(data_ess, data.frame(group_ind='total', t(colSums(data_ess[, -1]))))
    data_ess[nrow(data_ess), "ess"] = NA
    
    rownames(data_ess) = NULL 
    colnames(data_ess) = c("Group of Observation", 
                           "Actual Sample Size (N)", 
                           "Effective Sample Size (ESS)", 
                           "ESS/N", 
                           "Proportion of Information Borrowing")
    data_ess
  } else {
    if (!is.null(data$sbw_weights)) {
      data_ess = data %>% 
        group_by(group_ind) %>% 
        summarise(
          ss = n(), 
          ess = ESS(sbw_weights),
          ratio_ess_to_ss = ess / ss) %>% ungroup() %>% 
        mutate(percent_of_each = ess / sum(ess, na.rm = T))
      data_ess =  rbind(data_ess, data.frame(group_ind='total', t(colSums(data_ess[, -1], na.rm = T))))
      data_ess[nrow(data_ess), "ess"] = NA
      
      rownames(data_ess) = NULL 
      colnames(data_ess) = c("Group of Observation", 
                             "Actual Sample Size (N)", 
                             "Effective Sample Size (ESS)", 
                             "ESS/N", 
                             "Proportion of Information Borrowing")
      data_ess
    }
  }
}

PlotInfluence = function(data, inf_all, t0, t1, metric) {
  if (metric == "sic") {
    show.r <- sort(order(abs(inf_all$inf), decreasing = TRUE)[1:10])
  } else if (metric == "sic_scaled") {
    show.r <- sort(order(abs(inf_all$inf_scaled), decreasing = TRUE)[1:10])
  } else if (metric == "est_change") {
    show.r <- sort(order(abs(inf_all$est_change), decreasing = TRUE)[1:10])
  }
  labels1 <- rep("",nrow(data))
  labels1[show.r] <- as.character(show.r)
  
  plot_data = data %>% 
    mutate(SIC = inf_all$inf, 
           SIC.Scaled = inf_all$inf_scaled, 
           Est.Change = inf_all$est_change) %>% 
    dplyr::select(X, Unit, Time, Outcome, TreatStartTime, group_ind, SIC, SIC.Scaled, Est.Change)
  
  if (metric == "sic") {
    plt = ggplot(plot_data, aes(x = X, xend = X, y = 0, yend = SIC, 
                                label = Unit, label2 = Time, label3 = TreatStartTime, 
                                label4 = Outcome, label5 = group_ind, color = group_ind)) +
      geom_segment(linetype = "solid", linewidth = 0.3) + 
      labs(x = "Observation Index", y = "SIC", color = "Observation Group") +
      geom_text(aes(x = X, y = SIC, label = labels1), size = 3, vjust = -1) +
      scale_color_manual(values = c(IdealTreat = "#3293e3", 
                                    IdealControl = "#69b334", 
                                    InvTreat = "#c5dbed", 
                                    InvControl = "#ddedd1", 
                                    LimitAnticipate = "#ffc001", 
                                    DelayOnset = "#D2691E", 
                                    Dissipate = "#f56c5b",
                                    InvalidTreat = "#e1e1e1", 
                                    InvalidControl = "#fafafa")) +
      theme_bw()
    plt
    
  } else if (metric == "sic_scaled") {
    plt = ggplot(plot_data, aes(x = X, xend = X, y = 0, yend = SIC.Scaled, 
                                label = Unit, label2 = Time, label3 = TreatStartTime, 
                                label4 = Outcome, label5 = group_ind, color = group_ind)) +
      geom_segment(linetype = "solid", linewidth = 0.3) + 
      labs(x = "Observation Index", y = "SIC Scaled", color = "Observation Group") +
      geom_text(aes(x = X, y = SIC.Scaled, label = labels1), size = 3, vjust = -1) +
      scale_color_manual(values = c(IdealTreat = "#3293e3", 
                                    IdealControl = "#69b334", 
                                    InvTreat = "#c5dbed", 
                                    InvControl = "#ddedd1", 
                                    LimitAnticipate = "#ffc001", 
                                    DelayOnset = "#D2691E", 
                                    Dissipate = "#f56c5b",
                                    InvalidTreat = "#e1e1e1", 
                                    InvalidControl = "#fafafa")) +
      theme_bw()
    plt
    
  } else {
    plt = ggplot(plot_data, aes(x = X, xend = X, y = 0, yend = Est.Change, 
                                label = Unit, label2 = Time, label3 = TreatStartTime, 
                                label4 = Outcome, label5 = group_ind, color = group_ind)) +
      geom_segment(linetype = "solid", linewidth = 0.3) + 
      labs(x = "Observation Index", y = "Change in Point Estimate", color = "Observation Group") +
      geom_text(aes(x = X, y = Est.Change, label = labels1), size = 3, vjust = -1) +
      scale_color_manual(values = c(IdealTreat = "#3293e3", 
                                    IdealControl = "#69b334", 
                                    InvTreat = "#c5dbed", 
                                    InvControl = "#ddedd1", 
                                    LimitAnticipate = "#ffc001", 
                                    DelayOnset = "#D2691E", 
                                    Dissipate = "#f56c5b",
                                    InvalidTreat = "#e1e1e1", 
                                    InvalidControl = "#fafafa")) +
      theme_bw()
    plt
  }
}

GetInfluence = function(data, t0, t1, l_min, l_max) {
  
  unit_covariates <- c(paste("Unit", unique(data$Unit), sep = "_"))
  time_covariates <- c(paste("Time", unique(data$Time), sep = "_"))
  lead_lag_covariates <- c(paste("lead", seq(l_min,1), sep = "_"), 
                           paste("lag", seq(0,l_max), sep = "_"), 
                           "lag_inf")
  twfe_formula <- as.formula(
    paste("Outcome ~ ",
          paste(unit_covariates,collapse = " + "), " + ",
          paste(time_covariates,collapse = " + "), " + ",
          paste(lead_lag_covariates[-length(lead_lag_covariates)],collapse = " + ")
    )
  )
  target_name = ifelse(t1 >= t0, paste0("lag_", t1-t0, sep=""), paste0("lead_", t0-t1, sep=""))
  
  w_act = data$lmw_weight[data$treatment_component == 1]
  w_non = data$lmw_weight[data$treatment_component == 0]
  y_act = data$Outcome[data$treatment_component == 1]
  y_non = data$Outcome[data$treatment_component == 0]
  T_F = sum(y_act*w_act) - sum(y_non*w_non)
  
  ComputeInfluence_i = function(idx, t0, t1, data, target_name, formula, estimand = "ATT", method = "URI") {
    data_noi = data %>% filter(X != idx) 
    lmw_att_out_noi = lmw::lmw(formula, data = data_noi, estimand, method, treat = target_name)
    
    act_index_noi = data_noi$X[which(data_noi$time_til==t1-t0)]
    non_index_noi = data_noi$X[which(data_noi$time_til!=t1-t0)]
    w_noi_act = lmw_att_out_noi$weights[which(data_noi$X %in% act_index_noi)] / length(act_index_noi)
    w_noi_non = lmw_att_out_noi$weights[which(data_noi$X %in% non_index_noi)] / length(non_index_noi)
    y_noi_act = data_noi$Outcome[which(data_noi$X %in% act_index_noi)]
    y_noi_non = data_noi$Outcome[which(data_noi$X %in% non_index_noi)]
    T_F_noi = sum(y_noi_act*w_noi_act) - sum(y_noi_non*w_noi_non)
    inf_i = nrow(data_noi)*(T_F_noi-T_F)
    est_change_i = T_F_noi-T_F
    return(c(idx, inf_i, est_change_i))
  }
  
  # Here we use the parallel computing to speed up the influence computation 
  inf_all = parallel::mcmapply(ComputeInfluence_i, data$X,  
                               MoreArgs = list(t0 = t0, t1 = t1, data = data, 
                                               target_name = target_name, formula = twfe_formula), mc.cores = 8)
  inf_all = data.frame(t(inf_all))
  colnames(inf_all) = c("X", "inf", "est_change")
  inf_all$inf_scaled = abs(inf_all$inf)/max(abs(inf_all$inf))
  return(inf_all)
}

# to do, include balance method SIC values as well 
GetInfluenceBalance = function(estimand, method, t0, t1) {
  
  # compute the Balancing weights 
  ### Define moment covariates
  bal = list()
  bal$bal_cov = baseline_covariates
  ### Set tolerances
  bal$bal_std = "group"
  # Use adaptive algorithm 
  bal$bal_alg = TRUE
  bal$bal_gri = c(0.0001, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)
  
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
  weight_original = data_study_augment$sbw_weights
  act_index = which(data_study_augment$time_til==t1-t0)
  non_index = which(data_study_augment$time_til!=t1-t0)
  w_act = weight_original[act_index]
  w_non = weight_original[non_index]
  y_act = data_study_augment$asmrs[act_index]
  y_non = data_study_augment$asmrs[non_index]
  T_F = sum(y_act*w_act) - sum(y_non*w_non)
  
  ComputeInfluence_i = function(idx, t0, t1, data) {
    data_noi = data %>% filter(X != idx) 
    act_index_noi = data_noi$X[which(data_noi$time_til==t1-t0)]
    non_index_noi = data_noi$X[which(data_noi$time_til!=t1-t0)]
    
    data_subset = data %>% 
      filter(time_til == t1-t0 | treat_start_year == 9999) %>% 
      mutate(treat_ind = ifelse(time_til == t1-t0, 1, 0)) %>% 
      dplyr::select(X, treat_ind, any_of(baseline_covariates), asmrs)
    df = as.data.frame(cbind(data_subset$X, data_subset$treat_ind, 
                             data_subset[, baseline_covariates], data_subset$asmrs))
    names(df) = c("X", "t", baseline_covariates, "Y")
    
    # Solve for the Average Treatment Effect on the Treated, ATT (default)
    sbwatt_object = tryCatch(sbw(dat = df, ind = "t", out = "Y", bal = bal), error=function(e) {NULL})
    if (!is.null(sbwatt_object)) {
      sbwatt_weight_set = sbwatt_object$dat_weights %>% dplyr::select(X, sbw_weights)  
      data = data %>% left_join(sbwatt_weight_set, by = "X")
      data$sbw_weights = replace(data$sbw_weights, is.na(data$sbw_weights), 0)
      
      w_noi_act = data$sbw_weights[which(data_noi$X %in% act_index_noi)]
      w_noi_non = data$sbw_weights[which(data_noi$X %in% non_index_noi)]
      y_noi_act = data_noi$asmrs[which(data_noi$X %in% act_index_noi)]
      y_noi_non = data_noi$asmrs[which(data_noi$X %in% non_index_noi)]
      T_F_noi = sum(y_noi_act*w_noi_act) - sum(y_noi_non*w_noi_non)
      inf_i = nrow(data_noi)*(T_F_noi-T_F)
      est_change_i = T_F_noi-T_F
      
    } else {
      inf_i = est_change_i = NULL
    }
    
    return(c(inf_i, est_change_i))
  }
  
  balance_units = data_study_augment_subset$X
  inf_all = parallel::mcmapply(ComputeInfluence_i, balance_units,  
                               MoreArgs = list(t0 = t0, t1 = t1, 
                                               data = data_study_augment), 
                               mc.cores = 5)
  inf_all = data.frame(t(inf_all))
  colnames(inf_all) = c("inf", "est_change")
  inf_all$inf_scaled = abs(inf_all$inf)/max(abs(inf_all$inf))
  return(inf_all)
}

# GetObsGroup function takes in a dataset with 
# It returns an expanded dataset with a new column of observation group
# Observations groups can be: IdealTreat, IdealControl, InvTreat, InvControl, (DelayOnset), LimitAnticipate, Dissipate
# args: 
## data<data.frame>: should include columns of 
## t0<int>: indicate the treatment initiation time 
## t1<int>: indicate the outcome observation time 
## estimand<str>: can choose from std or new 

# returns: 
## data_obs_group<data.frame>: data + new column of `group_ind`
GetObsGroup = function(data, t0, t1, estimand = "std") {
  
  if (estimand == "std") {
    # Define different groups of observations 
    ideal_treat_set = data %>% filter(Time == t1 & TreatStartTime == t0)
    ideal_control_set = data %>% filter(Time == t1 & TreatStartTime == Inf) 
    invalid_set1 = data %>% filter(Time == t1 & !TreatStartTime %in% c(t0, Inf))
    
    inv_treat_set = data %>% filter(time_til == t1-t0) 
    inv_control_set = data %>% filter(TreatStartTime == Inf) 
    invalid_set2 = data %>% filter(time_til != t1-t0 & TreatStartTime != Inf) 
    
    ideal_treat_idx = ideal_treat_set$X
    ideal_control_idx = ideal_control_set$X
    inv_treat_idx = setdiff(inv_treat_set$X, ideal_treat_set$X)
    inv_control_idx = setdiff(inv_control_set$X, ideal_control_set$X)
    
    delay_set = data %>% 
      filter(Unit %in% inv_treat_set$Unit) %>% 
      filter(Treatment == 1 & time_til < t1-t0)
    
    anticipate_set = data %>% 
      filter(Unit %in% inv_treat_set$Unit) %>% 
      filter(Treatment == 0 & time_til < t1-t0)
    
    dissipate_set = data %>% 
      filter(Unit %in% inv_treat_set$Unit) %>% 
      filter(Treatment == 1 & time_til > t1-t0)
    
    treat_invalid_set3 = data %>% 
      filter(!Unit %in% c(inv_treat_set$Unit, inv_control_set$Unit)) %>% 
      filter(Treatment == 1)
    
    control_invalid_set3 = data %>% 
      filter(!Unit %in% c(inv_treat_set$Unit, inv_control_set$Unit)) %>% 
      filter(Treatment == 0)
    
    # Create one column indicating the group for each observation 
    data_obs_group = data %>% 
      mutate(group_ind = case_when(X %in% ideal_treat_idx ~ "IdealTreat",
                                   X %in% ideal_control_idx ~ "IdealControl",
                                   X %in% inv_treat_idx ~ "InvTreat", 
                                   X %in% inv_control_idx ~ "InvControl", 
                                   X %in% delay_set$X ~ "DelayOnset",
                                   X %in% anticipate_set$X ~ "LimitAnticipate",
                                   X %in% dissipate_set$X ~ "Dissipate", 
                                   X %in% treat_invalid_set3$X ~ "InvalidTreat", 
                                   X %in% control_invalid_set3$X ~ "InvalidControl")) %>% 
      mutate(group_ind = factor(group_ind, levels = c("IdealTreat", "IdealControl", 
                                                      "InvTreat", "InvControl",
                                                      "LimitAnticipate", "DelayOnset", "Dissipate", 
                                                      "InvalidTreat", "InvalidControl")))
  } else if (estimand == "new") { 
    
    # define the valid control lead lag range
    time_til_ubound <- t1-t0
    time_til_lbound <- t1-max(data$Time)
    
    # Define different groups of observations 
    ideal_treat_set = data %>% filter(Time == t1 & TreatStartTime == t0)
    ideal_control_set = data %>% filter(Time == t1 & TreatStartTime > t0) 
    invalid_set1 = data %>% filter(Time == t1 & TreatStartTime < t0)
    
    inv_treat_set = data %>% filter(time_til == t1-t0) 
    inv_control_set = data %>% filter( (time_til >= time_til_lbound & time_til <= time_til_ubound) | TreatStartTime == Inf)  
    
    ideal_treat_idx = ideal_treat_set$X
    ideal_control_idx = ideal_control_set$X
    inv_treat_idx = setdiff(inv_treat_set$X, ideal_treat_set$X)
    inv_control_idx = setdiff(inv_control_set$X, ideal_control_set$X)
    
    anticipate_set = data %>% 
      filter(Treatment == 0 & time_til < min(time_til_lbound:time_til_ubound) & TreatStartTime != Inf)
    
    dissipate_set = data %>% 
      filter(Treatment == 1 & time_til > max(time_til_lbound:time_til_ubound))
    
    treat_invalid_set3 = data %>% 
      filter(!X %in% c(inv_treat_set$X, inv_control_set$X, anticipate_set$X, dissipate_set$X)) %>% 
      filter(Treatment == 1)
    
    control_invalid_set3 = data %>% 
      filter(!X %in% c(inv_treat_set$X, inv_control_set$X, anticipate_set$X, dissipate_set$X)) %>% 
      filter(Treatment == 0)
    
    # Create one column indicating the group for each observation 
    data_obs_group = data %>% 
      mutate(group_ind = case_when(X %in% ideal_treat_idx ~ "IdealTreat",
                                   X %in% ideal_control_idx ~ "IdealControl",
                                   X %in% inv_treat_idx ~ "InvTreat", 
                                   X %in% inv_control_idx ~ "InvControl", 
                                   X %in% anticipate_set$X ~ "LimitAnticipate",
                                   X %in% dissipate_set$X ~ "Dissipate", 
                                   X %in% treat_invalid_set3$X ~ "InvalidTreat", 
                                   X %in% control_invalid_set3$X ~ "InvalidControl")) %>% 
      mutate(group_ind = factor(group_ind, levels = c("IdealTreat", "IdealControl", 
                                                      "InvTreat", "InvControl",
                                                      "LimitAnticipate", "Dissipate", 
                                                      "InvalidTreat", "InvalidControl")))
    
  } else {
    stop("Estimand to be specified as `std` or `new`. ")
  }
  return (data_obs_group)
}


GetImpliedWeights = function(data, t0, t1, l_min, l_max) {
  unit_covariates <- c(paste("Unit", unique(data$Unit), sep = "_"))
  time_covariates <- c(paste("Time", unique(data$Time), sep = "_"))
  lead_lag_covariates <- c(paste("lead", seq(l_min,1), sep = "_"), 
                           paste("lag", seq(0,l_max), sep = "_"), 
                           "lag_inf")
  # define the regression formula 
  twfe_formula <- as.formula(
    paste("Outcome ~ ",
          paste(unit_covariates,collapse = " + "), " + ",
          paste(time_covariates,collapse = " + "), " + ",
          paste(lead_lag_covariates[-length(lead_lag_covariates)],collapse = " + ")
    )
  )
  target_name = ifelse(t1 >= t0, paste0("lag_", t1-t0, sep=""), paste0("lead_", t0-t1, sep=""))
  lmw_att_out = lmw::lmw(twfe_formula, data = data, estimand = "ATT", method = "URI", treat = target_name)
  
  weight_original = lmw_att_out$weights
  treatment_component_idx = data$X[data$time_til == t1-t0]
  control_component_idx = data$X[data$time_til != t1-t0]
  
  data_weight = data %>% 
    # 2024/02/15: adapt to the changes of the lmw package
    mutate(treatment_component = ifelse(X %in% treatment_component_idx, 1, 0), 
           treatment_component = factor(treatment_component, levels = c(1, 0)), 
           lmw_weight = ifelse(treatment_component == 1, weight_original/length(treatment_component_idx), weight_original/length(control_component_idx))) %>% 
    mutate(uniform_weight = ifelse(treatment_component == 1, 1/length(treatment_component_idx), 1/length(control_component_idx))) 
  data_weight
}










