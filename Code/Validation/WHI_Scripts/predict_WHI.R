# Load packages
library(dplyr)
library(readr)
library(glmnet)

pacman::p_load(tidyverse, data.table, sjlabelled, caret, glmnet, purrr, readxl, writexl, patchwork, rlang)


# Load data
WHI_betas <- readRDS("../../WHI/WHI_Data/Initial/GRSet_fully_filtered_bmiq.rds")

# Load metadata and cvfit objects
metadata <- readRDS("../Models/enet_metadata_al_msa_probes_AL_mean.rds")
cvfit <- readRDS("../Models/enet_cvfit_al_msa_probes_AL_mean.rds")


########################################################
# Function to make predictions (for legacy)
########################################################

get_predictions <- function(betas, 
                            probe_metadata,
                            probe_cvfit,
                            prediction_name="prediction"){
  
  ## Assess probe missingness:
  print("Data dimensions:")
  dim(betas)
  print("\nFirst few column names:")
  head(colnames(betas))
  
  # Get coefficients for lambda.min
  coef_min <- as.matrix(coef(probe_cvfit, s = "lambda.min"))
  nonzero_coefs <- coef_min[coef_min != 0, , drop = FALSE]
  intercept <- nonzero_coefs["(Intercept)", ]
  feature_coefs <- nonzero_coefs[rownames(nonzero_coefs) != "(Intercept)", , drop = FALSE]
  
  print("Number of non-zero coefficients:")
  print(length(feature_coefs))
  print("Intercept value:")
  print(intercept)
  
  model_features <- rownames(feature_coefs)
  available_features <- intersect(model_features, rownames(betas))
  sample_ids <- colnames(betas)
  
  missing_features <- setdiff(model_features, available_features)
  
  print(paste("Total model features needed:", length(model_features)))
  print(paste("Available features in dataset:", length(available_features)))
  print(paste("Missing features:", length(missing_features)))
  print(paste("Percentage of model features available:", 
              round(length(available_features) / length(model_features) * 100, 2), "%"))
  
  # Calculate importance of missing features
  if(length(missing_features) > 0) {
    missing_coef_sum <- sum(abs(feature_coefs[missing_features, ]))
    total_coef_sum <- sum(abs(feature_coefs))
    print(paste("Sum of absolute coefficients for missing features:", missing_coef_sum))
    print(paste("Percentage of total coefficient weight missing:", 
                round(missing_coef_sum / total_coef_sum * 100, 2), "%"))
  }
  
  ########################################################
  # Manual approach for prediction
  ########################################################
  
  # Get the sample IDs
  sample_ids <- colnames(betas)
  
  # Extract only the features we need and their coefficients
  available_coefs <- feature_coefs[available_features, , drop = FALSE]
  
  # Convert to ensure it's a numeric matrix
  available_coefs <- as.matrix(available_coefs)
  
  # Extract the corresponding rows from the beta matrix
  beta_subset <- betas[available_features, , drop = FALSE]
  
  # Ensure it's a numeric matrix
  beta_subset <- as.matrix(beta_subset)
  
  # Replace NAs with zeros in one step
  beta_subset[is.na(beta_subset)] <- 0
  
  # Calculate predictions using matrix multiplication
  # Let's debug what's happening
  print("Dimensions of available_coefs:")
  print(dim(available_coefs))
  print("Class of available_coefs:")
  print(class(available_coefs))
  
  print("Dimensions of beta_subset:")
  print(dim(beta_subset))
  print("Class of beta_subset:")
  print(class(beta_subset))
  
  # Mmatrix multiplication
  predictions <- rep(intercept, ncol(beta_subset))
  for (i in 1:length(available_features)) {
    feature <- available_features[i]
    coef_value <- as.numeric(available_coefs[i, ])
    feature_values <- beta_subset[i, ]
    predictions <- predictions + (coef_value * feature_values)
  }
  names(predictions) <- sample_ids
  
  # Create results dataframe
  result_df <- data.frame(
    IDATid = sample_ids
  )
  result_df[[prediction_name]] <- predictions
  
  return(result_df)
}

########################################################
# Get results
########################################################

results <- get_predictions(WHI_betas,
                                  metadata,
                                  cvfit,
                                  prediction_name = "AL_mean")

########################################################
# Compute Chamberlain Clock
########################################################

source("../Scripts/methAL_calculator.R")

# Read the EPIC version table
meth_al_table <- read_excel("../Data/Supplementary Tables - CpGs.xlsx",
                            sheet = "Table 2.", skip = 1)

# Read the 450k compatible version
meth_alt_table <- read_excel("../Data/Supplementary Tables - CpGs.xlsx",
                             sheet = "Table 4.", skip = 1)

key_probes <- setdiff(
  unique(c(meth_al_table$`CpG name`, meth_alt_table$`CpG name`)),
  "(Intercept)"
)

chamberlain_values <-predict_cpg_scores(betas = WHI_betas, 
                                        model_tables = list(meth_al_table, meth_alt_table), 
                                        model_names = c("methAL", "methALT"), 
                                        id_colname = "Barcode")

########################################################
# Combine and export results
########################################################

merged_results <- merge(results, chamberlain_values, by.x="IDATid", by.y='Barcode')

write.csv(merged_results, "../Output/WHI_predictions.csv")
