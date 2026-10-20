predict_cpg_scores <- function(betas, model_tables, model_names = NULL, id_colname = "ID") {
  if (is.null(model_names)) {
    model_names <- paste0("Model", seq_along(model_tables))
  }
  
  # Extract all CpGs used across models
  all_cpgs <- unique(unlist(lapply(model_tables, function(tbl) {
    tbl$`CpG name`[grepl("^cg", tbl$`CpG name`)]
  })))
  
  # Filter betas once
  betas_subset <- betas[rownames(betas) %in% all_cpgs, , drop = FALSE]
  
  # Prepare output
  individual_predictions <- matrix(NA_real_, nrow = ncol(betas), ncol = length(model_tables))
  colnames(individual_predictions) <- model_names
  individual_ids <- colnames(betas)
  
  for (i in seq_along(model_tables)) {
    model_tbl <- as.data.frame(model_tables[[i]])
    model_tbl$Coefficient <- as.numeric(as.character(model_tbl$Coefficient))
    
    intercept <- model_tbl$Coefficient[model_tbl$`CpG name` == "(Intercept)"]
    model_cpgs <- model_tbl[grepl("^cg", model_tbl$`CpG name`), ]
    rownames(model_cpgs) <- model_cpgs$`CpG name`
    
    # Keep CpGs present in both model and data
    overlapping_cpgs <- intersect(rownames(betas_subset), rownames(model_cpgs))
    if (length(overlapping_cpgs) == 0) {
      warning(paste("No CpGs matched for model", model_names[i]))
      next
    }
    
    betas_model <- betas_subset[overlapping_cpgs, , drop = FALSE]
    coefs <- model_cpgs[overlapping_cpgs, "Coefficient"]
    
    # Prediction = weighted sum + intercept
    for (j in seq_len(ncol(betas_model))) {
      individual_predictions[j, i] <- sum(betas_model[, j] * coefs, na.rm = TRUE) + intercept
    }
    
    # Print coverage info
    message(sprintf("Your dataset has %d of %d CpGs for the %s measure.",
                    length(overlapping_cpgs), nrow(model_cpgs), model_names[i]))
  }
  
  # Return as tibble
  output_df <- as_tibble(individual_predictions)
  output_df[[id_colname]] <- individual_ids
  
  return(output_df)
}
