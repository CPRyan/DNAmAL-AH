# Load packages
library(readr)
library(dplyr)

# Load datasets
predictions <- read_csv("../Output/HRS_predictions.csv")
hrs_outcomes <- read_csv("../../HRS/HRS_Data/HRS_outcomes_and_clocks.csv")
hrs_disease_disability <- read_csv("../../HRS/HRS_Data/HRS_disease_disability_data.csv")
hrs_cognition <- read_csv("../../HRS/HRS_Data/HRS_cognition_data.csv")

# Simplify datasets
pheno <- hrs_cognition %>%
  select(c("id", "IDATid",
           "Race", "Female", "PAGE",
           "bmi1213", "smoker1213",
           "pwavessmoked", "n_smokedata",
           "Bas", "Bmem", "Bnv",
           "CD4mem", "CD4nv", "CD8mem", 
           "CD8nv", "Eos", "Mono", 
           "Neu", "NK", "Treg"))

other_clocks <- hrs_cognition %>%
  select(c("IDATid",
           "DunedinPACE", "DunedinPACE_scaled",
           "DunedinPACEResid", "DunedinPACEResid_scaled",
           "PCGrimAge", "PCGrimAge_scaled",
           "PCGrimAgeResid", "PCGrimAgeResid_scaled",
           "PCPhenoAge", "PCPhenoAge_scaled",
           "PCPhenoAgeResid", "PCPhenoAgeResid_scaled"))


surv_outcomes <- hrs_outcomes %>%
  select(c("IDATid", "death_status", "death_time"))

disease_disability_outcomes <- hrs_disease_disability %>%
  select(c("IDATid",
           "chrondxe13", "chrondxe14", "chrondxe15",
           "adl5a13", "adl5a14", "adl5a15",
           "iadl5a13"), "iadl5a14", "iadl5a15")
  
cog_outcomes <- hrs_cognition %>%
  select(c("IDATid",
           "cogfunction2016", "cogfunction2018", "cogfunction2020"))

# Merge datasets
surv_cog_outcomes <- merge(surv_outcomes, cog_outcomes, by="IDATid")
all_outcomes <- merge(surv_cog_outcomes, disease_disability_outcomes, by="IDATid")

pheno_outcomes <- merge(pheno, all_outcomes, by="IDATid")
HRS_outcome_data <- merge(pheno_outcomes, other_clocks, by="IDATid")

HRS_outcome_data <- merge(HRS_outcome_data, predictions, by="IDATid")

# Make the scaled/resid columns
HRS_outcome_data$DNAmAL_scaled <- scale(HRS_outcome_data$DNAmAL)
HRS_outcome_data$methAL_scaled <- scale(HRS_outcome_data$methAL)
HRS_outcome_data$methALT_scaled <- scale(HRS_outcome_data$methALT)

HRS_outcome_data <- HRS_outcome_data %>%
  mutate(DNAmALResid = resid(lm(DNAmAL ~ PAGE, data = .))) %>%
  mutate(methALResid = resid(lm(methAL ~ PAGE, data = .))) %>%
  mutate(methALTResid = resid(lm(methALT ~ PAGE, data = .)))

HRS_outcome_data$DNAmALResid_scaled <- scale(HRS_outcome_data$DNAmALResid)
HRS_outcome_data$methALResid_scaled <- scale(HRS_outcome_data$methALResid)
HRS_outcome_data$methALTResid_scaled <- scale(HRS_outcome_data$methALTResid)

# Write to CSV
write.csv(HRS_outcome_data, "../Output/HRS_outcome_data.csv")

