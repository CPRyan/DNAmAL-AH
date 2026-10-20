####################
# Set up library and cores:
####################
.libPaths("/ifs/sec/cpc/addhealth/users/cpr/R/x86_64-pc-linux-gnu-library/4.3")
# Where to find my custom package library

# Setup parallel computation - use all cores on our computer. (removed because it was resulting in errors)
# (num_cores = RhpcBLASctl::get_num_cores())
# options(mc.cores = num_cores)
options(expressions = 5e5)

# Setup paths
cpr_home <-"/ifs/sec/cpc/addhealth/users/cpr/AL_clock"
epig_home <-"/ifs/sec/cpc/addhealth/epig"

####################
# Load packages
####################
pacman::p_load(tidyverse, data.table, sjlabelled, caret, glmnet)

####################
# Setup Data
####################

###########
# Phenotype Data
###########
pheno <-read_csv(paste0(cpr_home, "/Output/Data/al_scores_cells_manifest_for_ewas.csv"))

print("pheno_loaded")
###########
# Load EWAS Betas
###########

# Save betas in EWAS and Sugden
betas <-read_csv(here::here("Output/Data", "betas_from_ewas_and_sugden_probes.csv"))
row_names <-betas$probe
betas <-as.data.frame(betas[,-1]); 
rownames(betas)<-row_names


print('betas_loaded and rownames fixed')
###########
# Load EWAS Betas
###########
load(file = here::here("Output/Data", "ewas_sugden_probe_list.Rdata"))

print('probelist loaded')


####################
# Match pheno rows with betas columns
####################
# Ensure all rows in pheno match column names in betas
pheno_filtered <- pheno %>% 
  filter(Barcode %in% colnames(betas)) %>%
  arrange(match(Barcode, colnames(betas)))  # Ensure correct order

# Reorder betas columns to match pheno$Barcode order
betas_reordered <- betas[, pheno_filtered$Barcode]

# Check if the order is now correctly aligned
alignment_check <-all(colnames(betas_reordered) == pheno_filtered$Barcode)
alignment_check
# Must be TRUE

rm(betas)
####################
# Setup models
####################
m1 <-readRDS(here::here("Output/Models", "enet_cvfit_al_msa_probes_AL_mean.rds"))

models <-list(m1)

# Setup probe list names for for loop
for_loop_probe_names <-c("al_msa_probes")

# Setup outcomes for naming predictions
outcomes <- c("AL_mean")

####################
# Run Machine Learning
####################

# Set seed to select folds
set.seed(1)

# Pick rows in pheno that correspond to training set. (should be the same as prior due to setting set)
train_rows <-sample(nrow(pheno_filtered), round(nrow(pheno_filtered)*0.8))

# Pick training and testing Y 
y_train <-pheno_filtered[train_rows,] # 3222 train
y_test  <-pheno_filtered[-train_rows,]# 806 test


pheno_filtered <-pheno_filtered %>%
  mutate(train_test = if_else(Barcode %in% y_train$Barcode, "train", "test"))

print("train and test constructed")


# Loop through models and their corresponding probe sets
for (i in seq_along(models)) {
  
  # Select the model (cvfit) and corresponding probe set
  cvfit <- models[[i]]  
  probe_set_name <- for_loop_probe_names[i]
  outcome_name <- outcomes[i]
  
  # Extract the corresponding probe set
  probe_set <- probes_list[[probe_set_name]]
  
  # Extract and transpose the matching X data
  x <- betas_reordered[rownames(betas_reordered) %in% probe_set, ] %>% t()
  
  # Make predictions using the stored glmnet model
  predicted_values <- as.vector(predict(cvfit, newx = x, s = "lambda.min"))
  
  # Generate column name dynamically
  col_name <- paste0(outcome_name, "_", probe_set_name, "_pred")
  
  # Ensure prediction length matches dataset
  if (length(predicted_values) == nrow(pheno_filtered)) {
    pheno_filtered[[col_name]] <- predicted_values
  } else {
    warning("Skipping ", col_name, " due to length mismatch.")
  }
  
  # Print status update
  message("Predictions saved for: ", col_name)
}

write_csv(pheno_filtered, file = here::here("Output/Data", "AL_phenos_with_predictions.csv"))


##################
# CV Plots
##################
pacman::p_load(glmnet)

for (i in seq_along(models)) {
  
  # Dynamically generate a filename
  file_path <- here::here("Output/Figures", paste0("cvfit_MSE_", for_loop_probe_names[i], "_", outcomes[i], ".pdf"))
  
  # Open PDF device
  pdf(file_path)
  
  # Create the plot
  plot(models[[i]])
  
  # Close the device
  dev.off()
  
  message("Saved plot: ", file_path)
}

