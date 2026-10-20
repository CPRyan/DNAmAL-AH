library(readr)
library(stringr)
library(survival)
library(ggplot2)
library(dplyr)

# I want to bring in the fixed outcomes
fhs_outcomes_data <- read_csv("../../FittingPoA/FHS_Data/adjusted_filtered_outcomes.csv")
fhs_pheno_data <- read_csv("../../FittingPoA/FHS_Data/Initial/framingham_offspring_clean_pheno_dataset.csv")
fhs_cells <- read_csv("../../FittingPoA/FHS_Data/Initial/Framingham_CPR_SALAS_cell_counts.csv")

# Add in the Add Health predictions
predictions <- read_csv("../Output/Framingham_predictions.csv")

# Merge all
pheno_outcomes <- merge(fhs_pheno_data, fhs_outcomes_data, by.x="dbGaP_Subject_ID", by.y="subject_id")
pheno_outcomes_cells <- merge(pheno_outcomes, fhs_cells, by="Barcode")
raw_survival_data <- merge(pheno_outcomes_cells, predictions, by.x="Barcode", by.y="IDATid")

names(raw_survival_data)

survival_data <- raw_survival_data %>%
  select("Barcode", "dbGaP_Subject_ID",
         "Age", "Female", "bmi_ln_adj_sc",
         "Bas", "Bmem", "Bnv",
         "CD4mem", "CD4nv", "CD8mem",
         "CD8nv", "Eos", "Mono",
         "Neu", "NK", "Treg",
         "death_status", "death_date",
         "chd_status", "chd_date",
         "chf_status", "chf_date",
         "cvd_status", "cvd_date",
         "stroke_status", "stroke_date",
         "dementia_status", "dementia_date",
         "AL_mean", "methAL",
         "methALT")


# Make the scaled/resid columns
survival_data$AL_mean_scaled <- scale(survival_data$AL_mean)
survival_data$methAL_scaled <- scale(survival_data$methAL)
survival_data$methALT_scaled <- scale(survival_data$methALT)

survival_data <- survival_data %>%
  mutate(AL_meanResid = resid(lm(AL_mean ~ Age, data = .))) %>%
  mutate(methALResid = resid(lm(methAL ~ Age, data = .))) %>%
  mutate(methALTResid = resid(lm(methALT ~ Age, data = .)))

survival_data$AL_meanResid_scaled <- scale(survival_data$AL_meanResid)
survival_data$methALResid_scaled <- scale(survival_data$methALResid)
survival_data$methALTResid_scaled <- scale(survival_data$methALTResid)


# Bring in other clocks from survival data
prev_survival_data <- read_csv("../../FittingPoA/FHS_Data/complete_survival_data.csv")
other_clocks <- prev_survival_data %>%
  select("Barcode",
         "DunedinPACE", "DunedinPACE_scaled", 
         "DunedinPACEResid", "DunedinPACEResid_scaled",
         "PCPhenoAge", "PCPhenoAge_scaled", 
         "PCPhenoAgeResid", "PCPhenoAgeResid_scaled",
         "PCGrimAge", "PCGrimAge_scaled",
         "PCGrimAgeResid", "PCGrimAgeResid_scaled"
  )

survival_data <- merge(survival_data, other_clocks, by="Barcode")

write.csv(survival_data, "../Output/Framingham_survival_data.csv")

#####################################
# Survival Analyses
#####################################

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
    
    coeffs <- summary(model)$coefficients
    conf_ints <- confint(model)
    
    # Get sample size info
    n_subjects <- nrow(data)
    n_events <- sum(data[[outcome_status]], na.rm = TRUE)
    
    # Create row of results
    clock_results <- data.frame(
      Dataset = "FHS",
      outcome = outcome_status,
      clock = clock,
      hazard_ratio = coeffs[clock, "exp(coef)"],
      lower_ci = exp(conf_ints[clock, "2.5 %"]),
      upper_ci = exp(conf_ints[clock, "97.5 %"]),
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

run_survival_analyses <- function(data, clocks, outcomes, control_covariates){
  
  combined_results <- data.frame(
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
  
  for (outcome in outcomes){
    
    outcome_status = str_c(outcome, "_status")
    outcome_time = str_c(outcome, "_date")
    
    outcome_results <- run_survival_analysis(
      data = data,
      clocks = clocks,
      outcome_time = outcome_time,
      outcome_status = outcome_status,
      control_covariates = control_covariates
    )
    
    combined_results <- rbind(combined_results, outcome_results)
  }
  
  return(combined_results)
}


clocks <- c("AL_meanResid_scaled", "methALResid_scaled", "methALTResid_scaled",
            "DunedinPACEResid_scaled", "PCGrimAgeResid_scaled", "PCPhenoAgeResid_scaled")

control_covariates <- c("Female", "Age")

outcomes <- c("death", "chd", "chf",
              "cvd", "stroke", "dementia")

results_all_nocells <- run_survival_analyses(
  data=survival_data,
  clocks=clocks,
  outcomes=outcomes,
  control_covariates = control_covariates
)

control_covariates_w_bmi <- c("Female", "Age", "bmi_ln_adj_sc")

results_all_bmi <- run_survival_analyses(
  data=survival_data,
  clocks=clocks,
  outcomes=outcomes,
  control_covariates = control_covariates_w_bmi
)

control_covariates_w_cells <- c("Female", "Age", 
                                "Bas", "Bmem", "Bnv",
                                "CD4mem", "CD4nv", "CD8mem",
                                "CD8nv", "Eos", "Mono",
                                "Neu", "NK", "Treg")

results_all_cells <- run_survival_analyses(
  data=survival_data,
  clocks=clocks,
  outcomes=outcomes,
  control_covariates = control_covariates_w_cells
)

control_covariates_w_bmi_and_cells <- c("Female", "Age", 
                                        "Bas", "Bmem", "Bnv",
                                        "CD4mem", "CD4nv", "CD8mem",
                                        "CD8nv", "Eos", "Mono",
                                        "Neu", "NK", "Treg",
                                        "bmi_ln_adj_sc")

results_all_bmi_and_cells <- run_survival_analyses(
  data=survival_data,
  clocks=clocks,
  outcomes=outcomes,
  control_covariates = control_covariates_w_bmi_and_cells
)

write.csv(results_all_nocells, "../Output/fhs_survival_results.csv")
write.csv(results_all_bmi, "../Output/fhs_survival_results_w_bmi.csv")
write.csv(results_all_cells, "../Output/fhs_survival_results_w_cells.csv")
write.csv(results_all_bmi_and_cells, "../Output/fhs_survival_results_w_bmi_cells.csv")
