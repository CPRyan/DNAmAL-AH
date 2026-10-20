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
pacman::p_load(tidyverse, data.table, sjlabelled, caret, glmnet, purrr, writexl, ComplexUpset)

####################
# Setup Data
####################

####################
# Setup models
####################
m1 <-readRDS(here::here("Output/Models", "enet_cvfit_al_msa_probes_AL_mean.rds"))

models <-list(m1)

# Setup probe list names for for loop
for_loop_probe_names <-c("al_msa_probes")

# Setup outcomes for naming predictions
outcomes <- c("AL_mean")


# Updated helper function to extract both summary and coefficients
extract_nonzero_with_weights <- function(model, outcome_name, probe_set_name) {
  coef_min <- coef(model, s = "lambda.min")
  
  # Convert sparse matrix to data frame
  coef_df <- as.data.frame(as.matrix(coef_min))
  coef_df <- tibble::rownames_to_column(coef_df, "feature")
  names(coef_df)[2] <- "coefficient"
  
  # Filter to non-zero coefficients only
  nonzero_coefs <- coef_df %>% filter(coefficient != 0)
  
  # Summary info (excluding intercept)
  n_probes <- sum(nonzero_coefs$feature != "(Intercept)")
  
  list(
    summary = tibble(
      model = probe_set_name,
      outcome = outcome_name,
      n_probes = n_probes
    ),
    coefs = nonzero_coefs
  )
}

# Run over all models
all_results <- pmap(
  list(models, outcomes, for_loop_probe_names),
  extract_nonzero_with_weights
)

# Combine summaries
summary_table <- bind_rows(map(all_results, "summary"))

# Build Excel sheet list
excel_sheets <- list("Summary" = summary_table)

for (i in seq_along(all_results)) {
  sheet_name <- paste0(
    str_replace_all(for_loop_probe_names[[i]], "_", " "),
    " - ",
    outcomes[[i]]
  )
  sheet_name <- str_trunc(sheet_name, 31, ellipsis = "")  # Excel sheet name limit
  excel_sheets[[sheet_name]] <- all_results[[i]]$coefs
}

# Write to Excel
write_xlsx(excel_sheets, path = here::here("Output/Tables", "enet_nonzero_probe_summary_and_weights.xlsx"))


########################
# Upset plot removed - single probe set only
########################

# Note: Upset plot not generated as we only have one probe set (al_msa_probes)
# Coefficient analysis provides complete probe set information

