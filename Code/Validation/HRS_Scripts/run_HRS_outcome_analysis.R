# Load packages
library(readr)
library(dplyr)
library(survival)
library(nnet)

# Read data
data <- read_csv("../Output/HRS_outcome_data.csv")

#######################################
# Bring in functions
#######################################
run_survival_analysis <- function(data, clocks, outcome_time, outcome_status, control_covariates){
  
  
  all_results <- data.frame(
    Dataset = character(),
    outcome = character(),
    clock = character(),
    hazard_ratio = numeric(),
    lower_ci = numeric(),
    upper_ci = numeric(),
    n_subjects = numeric(),
    n_events = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (clock in clocks){
    
    formula <- as.formula(paste("Surv(", outcome_time, ",", outcome_status, ") ~", 
                                paste(c(clock, control_covariates), collapse = " + ")))
    
    model <- suppressWarnings(coxph(formula, data = data))
    
    n_subjects <- model$n[1]  # The first element of model$n gives the number of subjects
    n_events <- model$nevent  # This gives the number of events in the model
    
    coeffs <- summary(model)[["coefficients"]]
    
    hr <- coeffs[clock, "exp(coef)"]
    
    conf_ints <- confint(model)
    lower_ci <- exp(conf_ints[clock, "2.5 %"])
    upper_ci <- exp(conf_ints[clock, "97.5 %"])
    
    # Add to results dataframe
    clock_results <- data.frame(
      Dataset = "HRS",
      outcome = "Death",
      clock = clock,
      hazard_ratio = hr,
      lower_ci = lower_ci,
      upper_ci = upper_ci,
      n_subjects = n_subjects,
      n_events = n_events,
      stringsAsFactors = FALSE
    )
    
    all_results <- rbind(all_results, clock_results)
    
  }
  
  # Print formatted table at the end
  for(i in 1:nrow(all_results)) {
    cat(sprintf("%-25s HR: %.4f (95%% CI: %.4f-%.4f)\n", 
                all_results$clock[i], 
                all_results$hazard_ratio[i], 
                all_results$lower_ci[i], 
                all_results$upper_ci[i]))
  }
  cat("\n")
  
  return(all_results)
  
}

get_pois_results_revised <- function(data, clocks, control_covariates, conf_level = 0.95){
  
  all_results = data.frame(
    Outcome = character(),
    Clock = character(),
    rate_ratio = numeric(),
    lower_ci = numeric(),
    upper_ci = numeric(),
    n_subjects = numeric(),
    n_events = numeric()
  )
  
  for (clock in clocks){
    
    #################################################
    # Chronic Disease - Conditional Prospective only
    #################################################
    
    chron_CP_formula <- as.formula(paste("chrondxe15", " ~ ", "chrondxe13 + ", 
                                         paste(c(clock, control_covariates), collapse = " + ")))
    
    chron_CP_model <- suppressWarnings(glm(chron_CP_formula, family = poisson(link = "log"), data = data))
    
    # Get coefficient and CI
    chron_CP_coeff <- exp(chron_CP_model[["coefficients"]][[clock]])
    chron_CP_ci <- suppressWarnings(exp(confint(chron_CP_model, parm = clock, level = conf_level)))
    
    
    n_subjects_chron <- nrow(chron_CP_model$model)
    
    complete_data_chron <- chron_CP_model$model
    n_events_chron <- sum(complete_data_chron$chrondxe15 > complete_data_chron$chrondxe13, na.rm = TRUE)
    
    # Create chronic disease row
    clock_CD_results <- data.frame(
      Outcome = "Chronic Disease",
      Clock = clock,
      rate_ratio = chron_CP_coeff,
      lower_ci = chron_CP_ci[1],
      upper_ci = chron_CP_ci[2],
      n_subjects = n_subjects_chron,
      n_events = n_events_chron
    )
    
    all_results <- rbind(all_results, clock_CD_results)
    
    #################################################
    # ADLs - Conditional Prospective only
    #################################################
    
    ADL_CP_formula <- as.formula(paste("adl5a15", " ~ ", "adl5a13 + ", 
                                       paste(c(clock, control_covariates), collapse = " + ")))
    
    ADL_CP_model <- suppressWarnings(glm(ADL_CP_formula, family = poisson(link = "log"), data = data))
    
    # Get coefficient and CI
    ADL_CP_coeff <- exp(ADL_CP_model[["coefficients"]][[clock]])
    ADL_CP_ci <- suppressWarnings(exp(confint(ADL_CP_model, parm = clock, level = conf_level)))
    
    n_subjects_ADL <- nrow(ADL_CP_model$model)
    
    complete_data_ADL <- ADL_CP_model$model
    n_events_ADL <- sum(complete_data_ADL$adl5a15 > complete_data_ADL$adl5a13, na.rm = TRUE)
    
    # Create ADL row
    clock_ADL_results <- data.frame(
      Outcome = "ADLs",
      Clock = clock,
      rate_ratio = ADL_CP_coeff,
      lower_ci = ADL_CP_ci[1],
      upper_ci = ADL_CP_ci[2],
      n_subjects = n_subjects_ADL,
      n_events = n_events_ADL
    )
    
    all_results <- rbind(all_results, clock_ADL_results)
    
    #################################################
    # IADLs - Conditional Prospective only
    #################################################
    
    IADL_CP_formula <- as.formula(paste("iadl5a15", " ~ ", "iadl5a13 + ", 
                                        paste(c(clock, control_covariates), collapse = " + ")))
    
    IADL_CP_model <- suppressWarnings(glm(IADL_CP_formula, family = poisson(link = "log"), data = data))
    
    # Get coefficient and CI
    IADL_CP_coeff <- exp(IADL_CP_model[["coefficients"]][[clock]])
    IADL_CP_ci <- suppressWarnings(exp(confint(IADL_CP_model, parm = clock, level = conf_level)))
    
    
    n_subjects_IADL <- nrow(IADL_CP_model$model)
    
    complete_data_IADL <- IADL_CP_model$model
    n_events_IADL <- sum(complete_data_IADL$iadl5a15 > complete_data_IADL$iadl5a13, na.rm = TRUE)
    
    
    # Create IADL row
    clock_IADL_results <- data.frame(
      Outcome = "IADLs",
      Clock = clock,
      rate_ratio = IADL_CP_coeff,
      lower_ci = IADL_CP_ci[1],
      upper_ci = IADL_CP_ci[2],
      n_subjects = n_subjects_IADL,
      n_events = n_events_IADL
    )
    
    all_results <- rbind(all_results, clock_IADL_results)
  }
  
  return(all_results)
}

multi_logit_restructured <- function(data, clocks, control_covariates){
  
  all_results = data.frame(
    clock = character(),
    outcome = character(),
    hazard_ratio = numeric(),
    lower_ci = numeric(),
    upper_ci = numeric(),
    n_subjects = numeric(),
    n_events = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Get the total number of subjects
  n_subjects <- nrow(data)
  
  for (clock in clocks){
    
    CP_formula <- as.formula(paste("cogfunction2020", " ~ cogfunction2016 + ",
                                   paste(c(clock, control_covariates), collapse = " + ")))
    
    CP_model <- suppressWarnings(multinom(CP_formula, data = data))
    
    model_data <- model.frame(CP_formula, data = data)
    n_subjects <- nrow(model_data)
    
    complete_data <- model_data[complete.cases(model_data), ]
    
    # Count transitions from lower to higher states (worsening)
    # For CIND: count cases where cogfunction2020 is 2 and cogfunction2016 was 1
    n_events_cind <- sum(complete_data$cogfunction2020 == 2 & complete_data$cogfunction2016 == 1, na.rm = TRUE)
    
    # For Dementia: count cases where cogfunction2020 is 3 and cogfunction2016 was 1 or 2
    n_events_dementia <- sum(complete_data$cogfunction2020 == 3 & 
                               (complete_data$cogfunction2016 == 1 | complete_data$cogfunction2016 == 2), 
                             na.rm = TRUE)
    
    # Extract coefficients and SEs for Level 2 vs Level 1
    CP_coeff_2 <- coef(CP_model)[1, clock]  # First row is Level 2 vs reference
    CP_se_2 <- summary(CP_model)$standard.errors[1, clock]
    
    # Calculate odds ratios and confidence intervals
    CP_OR_2 <- exp(CP_coeff_2)
    CP_lower_2 <- exp(CP_coeff_2 - (1.96 * CP_se_2))
    CP_upper_2 <- exp(CP_coeff_2 + (1.96 * CP_se_2))
    
    # Add the Level 2 vs Level 1 results to the dataframe
    level2_results <- data.frame(
      clock = clock,
      outcome = "CIND",
      hazard_ratio = CP_OR_2,
      lower_ci = CP_lower_2,
      upper_ci = CP_upper_2,
      n_subjects = n_subjects,
      n_events = n_events_cind,
      stringsAsFactors = FALSE
    )
    
    all_results <- rbind(all_results, level2_results)
    
    # Extract coefficients and SEs for Level 3 vs Level 1
    CP_coeff_3 <- coef(CP_model)[2, clock]  # Second row is Level 3 vs reference
    CP_se_3 <- summary(CP_model)$standard.errors[2, clock]
    
    # Calculate odds ratios and confidence intervals
    CP_OR_3 <- exp(CP_coeff_3)
    CP_lower_3 <- exp(CP_coeff_3 - (1.96 * CP_se_3))
    CP_upper_3 <- exp(CP_coeff_3 + (1.96 * CP_se_3))
    
    # Add the Level 3 vs Level 1 results to the dataframe
    level3_results <- data.frame(
      clock = clock,
      outcome = "Dementia",
      hazard_ratio = CP_OR_3,
      lower_ci = CP_lower_3,
      upper_ci = CP_upper_3,
      n_subjects = n_subjects,
      n_events = n_events_dementia,
      stringsAsFactors = FALSE
    )
    
    all_results <- rbind(all_results, level3_results)
    
  }
  
  return(all_results)
  
}


##########################################
# Get results (Minimal covariates)
##########################################

clocks <- c("DNAmALResid_scaled", "methALResid_scaled", 
            "methALTResid_scaled",
            "DunedinPACEResid_scaled", "PCGrimAgeResid_scaled", "PCPhenoAgeResid_scaled")

## Minimal covariates
minimal_covariates <- c("PAGE", "Female", "Race")

# Survival results
minimal_surv_results <- run_survival_analysis(
  data = data,
  clocks = clocks,
  outcome_time = "death_time",
  outcome_status = "death_status",
  control_covariates = minimal_covariates
)

# Disease/disability results
minimal_dd_results <- get_pois_results_revised(data, clocks, minimal_covariates)

# Cognition results
minimal_cog_results <- multi_logit_restructured(data, clocks, minimal_covariates)

##########################################
# Get results (Minimal covariates + cells)
##########################################

# Cell covariates
cell_covariates <- c("PAGE", "Female", "Race",
                     "Bas", "Bmem", "Bnv",
                     "CD4mem", "CD4nv", "CD8mem",
                     "CD8nv", "Eos", "Mono",
                     "Neu", "NK", "Treg")

# Survival results
cell_surv_results <- run_survival_analysis(
  data = data,
  clocks = clocks,
  outcome_time = "death_time",
  outcome_status = "death_status",
  control_covariates = cell_covariates
)

# Disease/disability results
cell_dd_results <- get_pois_results_revised(data, clocks, cell_covariates)

# Cognition results
cell_cog_results <- multi_logit_restructured(data, clocks, cell_covariates)

##########################################
# Get results (Minimal covariates + BMI)
##########################################

# BMI covariates
bmi_covariates <- c("PAGE", "Female", "Race", "bmi1213")

# Survival results
bmi_surv_results <- run_survival_analysis(
  data = data,
  clocks = clocks,
  outcome_time = "death_time",
  outcome_status = "death_status",
  control_covariates = bmi_covariates
)

# Disease/disability results
bmi_dd_results <- get_pois_results_revised(data, clocks, bmi_covariates)

# Cognition results
bmi_cog_results <- multi_logit_restructured(data, clocks, bmi_covariates)

##########################################
# Get results (Minimal covariates + cells + bmi)
##########################################

# BMI + cell covariates
bmi_cell_covariates <- c("PAGE", "Female", "Race",
                     "Bas", "Bmem", "Bnv",
                     "CD4mem", "CD4nv", "CD8mem",
                     "CD8nv", "Eos", "Mono",
                     "Neu", "NK", "Treg",
                     "bmi1213")

# Survival results
bmi_cell_surv_results <- run_survival_analysis(
  data = data,
  clocks = clocks,
  outcome_time = "death_time",
  outcome_status = "death_status",
  control_covariates = bmi_cell_covariates
)

# Disease/disability results
bmi_cell_dd_results <- get_pois_results_revised(data, clocks, bmi_cell_covariates)

# Cognition results
bmi_cell_cog_results <- multi_logit_restructured(data, clocks, bmi_cell_covariates)


##########################################
# Merge minimal results
##########################################

# Remove 'Dataset' column from the survival dataset
minimal_surv_results$Dataset <- NULL

# Need to adjust some more column names
minimal_dd_results <- minimal_dd_results %>%
  rename(outcome = Outcome,
         clock = Clock,
         hazard_ratio = rate_ratio)

# Merge datasets
minimal_results <- rbind(minimal_surv_results, minimal_dd_results, minimal_cog_results)

##########################################
# Merge cell results
##########################################

# Remove 'Dataset' column from the survival dataset
cell_surv_results$Dataset <- NULL

# Need to adjust some more column names
cell_dd_results <- cell_dd_results %>%
  rename(outcome = Outcome,
         clock = Clock,
         hazard_ratio = rate_ratio)

# Merge datasets
cell_results <- rbind(cell_surv_results, cell_dd_results, cell_cog_results)


##########################################
# Merge BMI results
##########################################

# Remove 'Dataset' column from the survival dataset
bmi_surv_results$Dataset <- NULL

# Need to adjust some more column names
bmi_dd_results <- bmi_dd_results %>%
  rename(outcome = Outcome,
         clock = Clock,
         hazard_ratio = rate_ratio)

# Merge datasets
bmi_results <- rbind(bmi_surv_results, bmi_dd_results, bmi_cog_results)

##########################################
# Merge bmi + cell results
##########################################

# Remove 'Dataset' column from the survival dataset
bmi_cell_surv_results$Dataset <- NULL

# Need to adjust some more column names
bmi_cell_dd_results <- bmi_cell_dd_results %>%
  rename(outcome = Outcome,
         clock = Clock,
         hazard_ratio = rate_ratio)

# Merge datasets
bmi_cell_results <- rbind(bmi_cell_surv_results, bmi_cell_dd_results, bmi_cell_cog_results)


##########################################
# Write results
##########################################

write.csv(minimal_results, "../Output/HRS_outcome_results.csv")
write.csv(cell_results, "../Output/HRS_outcome_results_w_cells.csv")
write.csv(bmi_results, "../Output/HRS_outcome_results_w_bmi.csv")
write.csv(bmi_cell_results, "../Output/HRS_outcome_results_w_bmi_cells.csv")
