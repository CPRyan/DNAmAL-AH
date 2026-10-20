# ============================================================
# Plotting function
# ============================================================

plot_sle_boxplot <- function(
    data,
    sle_var,
    y_var,
    sle_label,
    y_label
) {
  
  data %>%
    filter(
      !is.na(.data[[sle_var]]),
      !is.na(.data[[y_var]])
    ) %>%
    ggplot(
      aes(
        x = factor(.data[[sle_var]]),
        y = .data[[y_var]],
        fill = .data[[sle_var]]
      )
    ) +
    geom_boxplot(
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.10,
      shape = 21,
      colour = "black",
      alpha = 0.20
    ) +
    scale_fill_gradientn(
      colours = grDevices::hcl.colors(
        n = 10,
        palette = "Viridis",
        rev = TRUE
      )
    ) +
    labs(
      x = sle_label,
      y = y_label
    ) +
    theme_minimal() +
    theme(
      legend.position = "none"
    )
}