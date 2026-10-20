write_sjplot_xlsx <- function(
    model_sets,
    output_dir,
    xlsx_name = "SLE_DNAm_model_tables.xlsx",
    dv_label = "SLE Score (scaled)",
    pred_labels
) {
  
  required_packages <- c("sjPlot", "broom", "openxlsx")
  
  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
  
  if (length(missing_packages) > 0) {
    stop(
      "Install required packages: ",
      paste(missing_packages, collapse = ", ")
    )
  }
  
  if (is.null(names(model_sets))) {
    stop("model_sets must be a named list.")
  }
  
  # Create and then resolve the output directory once
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  output_dir <- normalizePath(
    output_dir,
    winslash = "/",
    mustWork = TRUE
  )
  
  # Formatting helpers
  format_number <- function(x, digits = 2) {
    ifelse(
      is.na(x),
      "",
      formatC(x, format = "f", digits = digits)
    )
  }
  
  format_p <- function(p) {
    
    stars <- ifelse(
      is.na(p), "",
      ifelse(
        p < 0.001, "***",
        ifelse(
          p < 0.01, "**",
          ifelse(p < 0.05, "*", "")
        )
      )
    )
    
    value <- ifelse(
      is.na(p), "",
      ifelse(
        p < 0.001,
        "<0.001",
        formatC(p, format = "f", digits = 3)
      )
    )
    
    paste0(value, stars)
  }
  
  # Excel workbook
  wb <- openxlsx::createWorkbook()
  
  title_style <- openxlsx::createStyle(
    fontSize = 13,
    textDecoration = "bold",
    halign = "left"
  )
  
  model_header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    fgFill = "#D9E2F3",
    border = "bottom",
    borderStyle = "medium"
  )
  
  column_header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    fgFill = "#E7E6E6",
    border = "bottom",
    borderStyle = "thin"
  )
  
  body_style <- openxlsx::createStyle(
    halign = "center",
    valign = "center"
  )
  
  first_column_style <- openxlsx::createStyle(
    halign = "left",
    valign = "center"
  )
  
  summary_style <- openxlsx::createStyle(
    textDecoration = "bold",
    border = "top",
    borderStyle = "thin"
  )
  
  for (table_title in names(model_sets)) {
    
    models <- model_sets[[table_title]]
    n_models <- length(models)
    
    # Coefficient order across the three models
    term_order <- unique(
      unlist(
        lapply(
          models,
          function(model) names(stats::coef(model))
        )
      )
    )
    
    if (length(pred_labels) != length(term_order)) {
      stop(
        "For ", table_title, ", pred_labels has ",
        length(pred_labels), " labels but the models have ",
        length(term_order), " unique coefficients:\n",
        paste(term_order, collapse = ", ")
      )
    }
    
    # Name the labels so sjPlot matches them to coefficients
    named_pred_labels <- stats::setNames(
      pred_labels,
      term_order
    )
    
    # --------------------------------------------------------
    # HTML table
    # --------------------------------------------------------
    
    file_stub <- gsub(
      "[^A-Za-z0-9]+",
      "_",
      table_title
    )
    
    file_stub <- gsub(
      "^_|_$",
      "",
      file_stub
    )
    
    html_file <- file.path(
      output_dir,
      paste0(file_stub, "_models.html")
    )
    
    sjPlot::tab_model(
      models,
      title = table_title,
      dv.labels = rep(dv_label, n_models),
      pred.labels = named_pred_labels,
      show.ci = 0.95,
      show.se = FALSE,
      show.stat = TRUE,
      show.r2 = TRUE,
      show.obs = TRUE,
      p.style = "numeric_stars",
      digits = 2,
      digits.p = 3,
      file = html_file
    )
    
    # --------------------------------------------------------
    # Excel table
    # --------------------------------------------------------
    
    model_results <- lapply(
      models,
      broom::tidy,
      conf.int = TRUE,
      conf.level = 0.95
    )
    
    n_term_rows <- length(term_order)
    n_summary_rows <- 3
    n_excel_columns <- 1 + n_models * 4
    
    table_body <- matrix(
      "",
      nrow = n_term_rows + n_summary_rows,
      ncol = n_excel_columns
    )
    
    table_body[, 1] <- c(
      unname(named_pred_labels[term_order]),
      "R²",
      "Adjusted R²",
      "Observations"
    )
    
    for (i in seq_along(models)) {
      
      result <- model_results[[i]]
      
      # Align estimates with the common term order
      result <- result[
        match(term_order, result$term),
        ,
        drop = FALSE
      ]
      
      first_model_column <- 2 + (i - 1) * 4
      
      estimate_column <- first_model_column
      ci_column <- first_model_column + 1
      statistic_column <- first_model_column + 2
      p_column <- first_model_column + 3
      
      table_body[
        seq_len(n_term_rows),
        estimate_column
      ] <- format_number(result$estimate)
      
      table_body[
        seq_len(n_term_rows),
        ci_column
      ] <- ifelse(
        is.na(result$conf.low),
        "",
        paste0(
          format_number(result$conf.low),
          " – ",
          format_number(result$conf.high)
        )
      )
      
      table_body[
        seq_len(n_term_rows),
        statistic_column
      ] <- format_number(result$statistic)
      
      table_body[
        seq_len(n_term_rows),
        p_column
      ] <- format_p(result$p.value)
      
      model_glance <- broom::glance(models[[i]])
      
      table_body[
        n_term_rows + 1,
        estimate_column
      ] <- format_number(
        model_glance$r.squared,
        digits = 3
      )
      
      table_body[
        n_term_rows + 2,
        estimate_column
      ] <- format_number(
        model_glance$adj.r.squared,
        digits = 3
      )
      
      table_body[
        n_term_rows + 3,
        estimate_column
      ] <- as.character(stats::nobs(models[[i]]))
    }
    
    sheet_name <- gsub(
      "[:\\\\/?*\\[\\]]",
      "_",
      table_title
    )
    
    sheet_name <- substr(sheet_name, 1, 31)
    
    openxlsx::addWorksheet(
      wb,
      sheetName = sheet_name,
      gridLines = FALSE,
      zoom = 90
    )
    
    # Table title
    openxlsx::writeData(
      wb,
      sheet = sheet_name,
      x = matrix(table_title, nrow = 1),
      startRow = 1,
      startCol = 1,
      colNames = FALSE
    )
    
    openxlsx::mergeCells(
      wb,
      sheet = sheet_name,
      cols = 1:n_excel_columns,
      rows = 1
    )
    
    # Model headers
    model_headers <- c(
      "",
      unlist(
        lapply(
          seq_len(n_models),
          function(x) c(dv_label, "", "", "")
        )
      )
    )
    
    column_headers <- c(
      "Predictor",
      rep(
        c("Estimate", "95% CI", "Statistic", "p"),
        n_models
      )
    )
    
    openxlsx::writeData(
      wb,
      sheet = sheet_name,
      x = matrix(model_headers, nrow = 1),
      startRow = 2,
      startCol = 1,
      colNames = FALSE
    )
    
    openxlsx::writeData(
      wb,
      sheet = sheet_name,
      x = matrix(column_headers, nrow = 1),
      startRow = 3,
      startCol = 1,
      colNames = FALSE
    )
    
    # Merge each dependent-variable heading across four columns
    for (i in seq_len(n_models)) {
      
      first_col <- 2 + (i - 1) * 4
      last_col <- first_col + 3
      
      openxlsx::mergeCells(
        wb,
        sheet = sheet_name,
        cols = first_col:last_col,
        rows = 2
      )
    }
    
    openxlsx::writeData(
      wb,
      sheet = sheet_name,
      x = as.data.frame(
        table_body,
        stringsAsFactors = FALSE
      ),
      startRow = 4,
      startCol = 1,
      colNames = FALSE,
      rowNames = FALSE
    )
    
    final_row <- 3 + nrow(table_body)
    summary_first_row <- 4 + n_term_rows
    
    openxlsx::addStyle(
      wb,
      sheet = sheet_name,
      style = title_style,
      rows = 1,
      cols = 1:n_excel_columns,
      gridExpand = TRUE
    )
    
    openxlsx::addStyle(
      wb,
      sheet = sheet_name,
      style = model_header_style,
      rows = 2,
      cols = 1:n_excel_columns,
      gridExpand = TRUE
    )
    
    openxlsx::addStyle(
      wb,
      sheet = sheet_name,
      style = column_header_style,
      rows = 3,
      cols = 1:n_excel_columns,
      gridExpand = TRUE
    )
    
    openxlsx::addStyle(
      wb,
      sheet = sheet_name,
      style = body_style,
      rows = 4:final_row,
      cols = 2:n_excel_columns,
      gridExpand = TRUE
    )
    
    openxlsx::addStyle(
      wb,
      sheet = sheet_name,
      style = first_column_style,
      rows = 4:final_row,
      cols = 1
    )
    
    openxlsx::addStyle(
      wb,
      sheet = sheet_name,
      style = summary_style,
      rows = summary_first_row:final_row,
      cols = 1:n_excel_columns,
      gridExpand = TRUE,
      stack = TRUE
    )
    
    openxlsx::setColWidths(
      wb,
      sheet = sheet_name,
      cols = 1,
      widths = 24
    )
    
    openxlsx::setColWidths(
      wb,
      sheet = sheet_name,
      cols = 2:n_excel_columns,
      widths = 13
    )
    
    openxlsx::freezePane(
      wb,
      sheet = sheet_name,
      firstActiveRow = 4,
      firstActiveCol = 2
    )
  }
  
  xlsx_file <- file.path(
    output_dir,
    xlsx_name
  )
  
  openxlsx::saveWorkbook(
    wb,
    file = xlsx_file,
    overwrite = TRUE
  )
  
  message("Tables written to: ", output_dir)
  
  invisible(xlsx_file)
}