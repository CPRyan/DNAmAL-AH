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
pacman::p_load(tidyverse, data.table, sjlabelled, SuperLearner, caret, glmnet)

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
if(alignment_check)print("Alignment of betas and phenotype is done, data ready to proceed")else("Problems with alignment - end run")


####################
# Run Machine Learning
####################

# Set seed to select folds
set.seed(1)

# Pick rows in pheno that correspond to training set.
train_rows <-sample(nrow(pheno_filtered), round(nrow(pheno_filtered)*0.8))

# Pick training and testing Y
y_train <-pheno_filtered[train_rows,] # 3222 train
y_test  <-pheno_filtered[-train_rows,]# 806 test

# Save train/test split for later
save(y_train, y_test, file = here::here("Output/Models", "y_train_test_data_split.rda"))

# Pick training and testing X
x_train <-betas_reordered[,y_train$Barcode]
x_test <-betas_reordered[,y_test$Barcode]

print("train and test constructed")



outcomes <- c("AL_mean")

# Loop over each probe set in probes_list
for (probe_set_name in names(probes_list)) {
  
  probe_set <- probes_list[[probe_set_name]]  # Extract the actual probe list
  
  # Loop over each outcome
  for (o in outcomes) {
    x <- x_train[rownames(x_train) %in% probe_set, ] %>% t()
    y <- y_train[[o]]
    
    xtest <-x_test[rownames(x_test) %in% probe_set, ] %>% t()
    ytest <-y_test[[o]]
    
    set.seed(1) # Ensure reproducibility across cores
    
    cv_fit <- cv.glmnet(x = x,
                        y = y,
                        alpha = 0.5)
    
    # Dynamically generate a filename with both the probe set name and outcome
    file_path <- here::here("Output/Models", paste0("enet_cvfit_", probe_set_name, "_", o, ".rds"))
    
    # Save the model output
    saveRDS(cv_fit, file = file_path)
    
    # Print status update
    message("Saved model for probe set: ", probe_set_name, " and outcome: ", o, " at ", file_path)
    ###################
    # Metadata
    ################### 
    metadata <- list(
      # Feature names in correct order
      feature_names = colnames(x),
      
      # Model parameters
      lambda_min = cv_fit$lambda.min,
      lambda_1se = cv_fit$lambda.1se,
      # Coefficients for both lambda values
      coef_min = as.matrix(coef(cv_fit, s = 'lambda.min')),
      coef_1se = as.matrix(coef(cv_fit, s = 'lambda.1se')),
      # Training data characteristics
      n_features = ncol(x),
      feature_means = colMeans(x),
      feature_sds = apply(x, 2, sd),
      # Model performance metrics
      r2_min = cor(predict(cv_fit, newx = xtest, s = 'lambda.min'),
                   ytest)^2,
      r2_1se = cor(predict(cv_fit, newx = xtest, s = 'lambda.1se'),
                   ytest)^2,
      # Time stamp
      save_date = Sys.time()
    )
    
    # Dynamically generate a filename with both the probe set name and outcome
    file_path2 <- here::here("Output/Models", paste0("enet_metadata_", probe_set_name, "_", o, ".rds"))
    
    # Save the model output
    saveRDS(metadata, file = file_path2)
    
  }
}

