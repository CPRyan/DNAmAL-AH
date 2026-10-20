## Load packages
library(readr)
library(dplyr)
library(ggplot2)

## Load data
data <- read_csv("../Output/fhs_survival_results.csv")

## Plotting Function

outcomes <- c("death_status",
              "chd_status",
              "chf_status",
              "cvd_status",
              "stroke_status",
              "dementia_status")

outcome_names <- c("Death",
                   "CHD",
                   "CHF",
                   "CVD",
                   "Stroke",
                   "Dementia")

## Create a named vector for easy lookup
outcome_labels <- setNames(outcome_names, outcomes)

# Define all available clocks with their colors
all_clock_colors <- c(
  "AL_meanResid_scaled" = "navy", 
  "methALResid_scaled" = "blue", 
  "methALTResid_scaled" = "#4c8db3", 
  "DunedinPACEResid_scaled" = "#a13cc9", 
  "PCPhenoAgeResid_scaled" = "#f37577", 
  "PCGrimAgeResid_scaled" = "#b31012"
)

all_clock_labels <- c(
  "AL_meanResid_scaled" = "DNAmAL",
  "methALResid_scaled" = "Chamberlain MethAL",
  "methALTResid_scaled" = "Chamberlain MethALT",
  "DunedinPACEResid_scaled" = "DunedinPACE",
  "PCPhenoAgeResid_scaled" = "PCPhenoAge",
  "PCGrimAgeResid_scaled" = "PCGrimAge"
)

plot_all_survival_results_dot <- function(results_data,
                                          selected_clocks = NULL,
                                          title = NULL) {
  
  if (is.null(title)){
    title = "Hazard Ratios by Outcome and Clock in FHS"
  }
  
  # If no clocks specified, use all of them
  if (is.null(selected_clocks)) {
    selected_clocks <- names(all_clock_labels)
  }
  
  # Subset the clock labels and colors to only include selected clocks
  clock_labels <- all_clock_labels[selected_clocks]
  clock_colors <- all_clock_colors[selected_clocks]
  
  # Ensure the order is preserved
  clock_labels <- clock_labels[order(match(names(clock_labels), selected_clocks))]
  clock_colors <- clock_colors[order(match(names(clock_colors), selected_clocks))]
  
  # First, prepare the data by mapping outcome and clock names
  plot_data <- results_data %>%
    # Make sure outcome values are in the outcomes vector
    filter(outcome %in% outcomes) %>%
    # Filter to only include the selected clocks
    filter(clock %in% selected_clocks) %>%
    # Create display names for outcomes and clocks
    mutate(
      outcome_name = outcome_labels[outcome],
      clock_name = clock_labels[clock],
      # Order outcomes for consistent display
      outcome_name = factor(outcome_name, levels = outcome_names),
      # Make clocks a factor for consistent ordering based on selected_clocks order
      clock_name = factor(clock_name, levels = unname(clock_labels[selected_clocks]))
    )
  
  # Get event counts for each outcome to display in panels
  event_counts <- plot_data %>%
    group_by(outcome, outcome_name) %>%
    summarise(n_events = first(n_events), .groups = "drop") %>%
    mutate(panel_label = paste0(outcome_name, "\n(n = ", n_events, " events)"))
  
  # Join the event counts back to the plot data
  plot_data <- plot_data %>%
    left_join(event_counts %>% dplyr::select(outcome, panel_label), by = "outcome")
  
  # Create the plot
  p <- ggplot(plot_data) +
    # Reference line at HR = 1 (moved to back)
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    # Add error bars with matching colors and increased thickness
    geom_errorbar(aes(x = clock_name, y = hazard_ratio, 
                      ymin = lower_ci, ymax = upper_ci, 
                      color = clock_name), 
                  width = 0.15, linewidth = 0.8) +
    # Use solid colored points with no outline
    geom_point(aes(x = clock_name, y = hazard_ratio, color = clock_name), size = 4) +
    # Set colors manually - match by position instead of by name
    scale_color_manual(values = unname(clock_colors)) +
    # Set y-axis limits based on data range with a bit of padding
    scale_y_continuous(limits = c(0.9, 2.0)) +
    # Create the faceted layout using panel labels with event counts
    facet_wrap(~ panel_label, nrow = 2, ncol = 3) +
    # Labels
    labs(
      title = title,
      x = NULL,
      y = "Hazard Ratio (95% CI)"
    ) +
    # Theme
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 9),
      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80", linewidth = 0.5),
      panel.spacing = unit(1, "lines")
    )
  
  return(p)
}

# Example usage:

# Create plot with all clocks (default behavior)
all_survival_plot <- plot_all_survival_results_dot(data)
print(all_survival_plot)


# Save plots
ggsave(
  filename = "../Output/fhs_survival_plot.pdf", 
  plot = all_survival_plot,
  width = 12,
  height = 8,
  dpi = 600,
  units = "in"
)

saveRDS(all_survival_plot, "../Output/fhs_survival_plot.rds")

