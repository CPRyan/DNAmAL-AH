# Load packages
library(readr)

## Load in all datasets

# Survival
survival_minimal <- read_csv("../Output/whi_survival_results.csv")
survival_cells <- read_csv("../Output/whi_survival_results_w_cells.csv")
survival_bmi <- read_csv("../Output/whi_survival_results_w_bmi.csv")
survival_bmi_cells <- read_csv("../Output/whi_survival_results_w_bmi_cells.csv")

###########################################
# Get clock-specific results
###########################################

# Filter each dataframe to only include AL_meanResid_scaled entries
baseline_data <- survival_minimal %>% 
  filter(clock == "AL_meanResid_scaled") %>%
  mutate(model = "Baseline")  # Add a new column for the model type

cells_data <- survival_cells %>%
  filter(clock == "AL_meanResid_scaled") %>%
  mutate(model = "With Cells")

bmi_data <- survival_bmi %>%
  filter(clock == "AL_meanResid_scaled") %>%
  mutate(model = "With BMI")

bmi_cells_data <- survival_bmi_cells %>%
  filter(clock == "AL_meanResid_scaled") %>%
  mutate(model = "With BMI + Cells")

# Combine all filtered data into one dataframe
sensitivity_results <- rbind(baseline_data, cells_data, bmi_data, bmi_cells_data)

# Reorder columns to make 'model' appear earlier for better readability
sensitivity_results <- sensitivity_results %>%
  dplyr::select(outcome, model, hazard_ratio, lower_ci, upper_ci, n_subjects, n_events, clock)

DNAmAL_sensitivity_results <- sensitivity_results

################
################
################

# Filter each dataframe to only include methALResid_scaled entries
baseline_data <- survival_minimal %>% 
  filter(clock == "methALResid_scaled") %>%
  mutate(model = "Baseline")  # Add a new column for the model type

cells_data <- survival_cells %>%
  filter(clock == "methALResid_scaled") %>%
  mutate(model = "With Cells")

bmi_data <- survival_bmi %>%
  filter(clock == "methALResid_scaled") %>%
  mutate(model = "With BMI")

bmi_cells_data <- survival_bmi_cells %>%
  filter(clock == "methALResid_scaled") %>%
  mutate(model = "With BMI + Cells")

# Combine all filtered data into one dataframe
sensitivity_results <- rbind(baseline_data, cells_data, bmi_data, bmi_cells_data)

# Reorder columns to make 'model' appear earlier for better readability
sensitivity_results <- sensitivity_results %>%
  dplyr::select(outcome, model, hazard_ratio, lower_ci, upper_ci, n_subjects, n_events, clock)

methAL_sensitivity_results <- sensitivity_results

################
################
################

# Filter each dataframe to only include methALTResid_scaled entries
baseline_data <- survival_minimal %>% 
  filter(clock == "methALTResid_scaled") %>%
  mutate(model = "Baseline")  # Add a new column for the model type

cells_data <- survival_cells %>%
  filter(clock == "methALTResid_scaled") %>%
  mutate(model = "With Cells")

bmi_data <- survival_bmi %>%
  filter(clock == "methALTResid_scaled") %>%
  mutate(model = "With BMI")

bmi_cells_data <- survival_bmi_cells %>%
  filter(clock == "methALTResid_scaled") %>%
  mutate(model = "With BMI + Cells")

# Combine all filtered data into one dataframe
sensitivity_results <- rbind(baseline_data, cells_data, bmi_data, bmi_cells_data)

# Reorder columns to make 'model' appear earlier for better readability
sensitivity_results <- sensitivity_results %>%
  dplyr::select(outcome, model, hazard_ratio, lower_ci, upper_ci, n_subjects, n_events, clock)

methALT_sensitivity_results <- sensitivity_results

###########################################
# Plotting stuff
###########################################

# Plot function
outcomes <- c("DEATH",
              "ANYCANCER",
              "CHD",
              "STROKE")

outcome_names <- c("DEATH",
                   "ANYCANCER",
                   "CHD",
                   "STROKE")

## Create a named vector for easy lookup
outcome_labels <- setNames(outcome_names, outcomes)

model_colors <- c("#096192", "#1171ba", "#1399c6", "#24aae2")

## Create better display names for the clocks
model_labels <- c(
  "Baseline" = "Baseline",
  "With Cells" = "With Cells",
  "With BMI" = "With BMI",
  "With BMI + Cells" = "With BMI + Cells"
)


plot_sensitivity_analyses_dot <- function(results_data, title=NULL) {
  
  if (is.null(title)){
    title = "Sensitivity Analysis in HRS"
  }
  
  
  # Hardcoded outcome order
  sorted_outcomes <- c("DEATH", "CHD", "STROKE", "ANYCANCER")
  
  # First, prepare the data by mapping outcome and clock names
  plot_data <- results_data %>%
    # Make sure outcome values are in the outcomes vector
    filter(outcome %in% outcomes) %>%
    # Filter to only include the clocks we want
    filter(model %in% names(model_labels)) %>%
    # Create display names for outcomes and clocks
    mutate(
      outcome_name = outcome_labels[outcome],
      model_name = model_labels[model],
      # Order outcomes based on hardcoded order
      outcome = factor(outcome, levels = sorted_outcomes),
      outcome_name = factor(outcome_name, levels = outcome_labels[sorted_outcomes]),
      # Make clocks a factor for consistent ordering
      model_name = factor(model_name, levels = unname(model_labels))
    )
  
  # Get event counts for each outcome to display in panels
  event_counts <- plot_data %>%
    group_by(outcome, outcome_name) %>%
    summarise(n_events = first(n_events), .groups = "drop")
  
  # Create panel labels conditionally based on whether n_events is available
  event_counts <- event_counts %>%
    mutate(panel_label = case_when(
      is.na(n_events) ~ as.character(outcome_name),
      TRUE ~ paste0(outcome_name, "\n(n = ", n_events, " events)")
    ))
  
  # Order the panel labels
  event_counts$panel_label <- factor(event_counts$panel_label, 
                                     levels = event_counts$panel_label[order(event_counts$outcome)])
  
  # Join the event counts back to the plot data
  plot_data <- plot_data %>%
    left_join(event_counts %>% dplyr::select(outcome, panel_label), by = "outcome")
  
  # Create the plot
  p <- ggplot(plot_data) +
    # Reference line at HR = 1 (moved to back)
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    # Add error bars with matching colors and increased thickness
    geom_errorbar(aes(x = model_name, y = hazard_ratio, 
                      ymin = lower_ci, ymax = upper_ci, 
                      color = model_name), 
                  width = 0.15, linewidth = 0.8) +
    # Use solid colored points with no outline
    geom_point(aes(x = model_name, y = hazard_ratio, color = model_name), size = 4) +
    # Set colors manually 
    scale_color_manual(values = model_colors) +
    # Set y-axis limits based on data range with a bit of padding
    scale_y_continuous(limits = c(0.8, 1.7)) +
    # Create the faceted layout using panel labels with event counts
    facet_wrap(~ panel_label, nrow = 2, ncol = 2) +
    # Labels
    labs(
      title = title,
      x = NULL,
      y = "Effect Size (95% CI)"
    ) +
    # Theme
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24, margin = margin(b = 15)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.text.y = element_text(size = 9),
      axis.title.y = element_text(size = 18, margin = margin(r = 20, l = 10)),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(size = 16, face = "bold"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80", linewidth = 0.5),
      panel.spacing = unit(1, "lines")
    )
  
  return(p)
}

###########################################
# Make the plots
###########################################

DNAmAL_sensitivity_plots <- plot_sensitivity_analyses_dot(DNAmAL_sensitivity_results,
                                                          title="DNAmAL Sensitivity Analysis in FHS")
print(DNAmAL_sensitivity_plots)

methAL_sensitivity_plots <- plot_sensitivity_analyses_dot(methAL_sensitivity_results,
                                                          title="Chamberlain methAL Sensitivity Analysis in FHS")
print(methAL_sensitivity_plots)

methALT_sensitivity_plots <- plot_sensitivity_analyses_dot(methALT_sensitivity_results,
                                                           title="Chamberlain methALT Sensitivity Analysis in FHS")
print(methALT_sensitivity_plots)

###########################################
# Save the plots
###########################################
ggsave(
  filename = "../Output/whi_DNAmAL_sensitivity_results.pdf", 
  plot = DNAmAL_sensitivity_plots,
  width = 10,  # Width in inches
  height = 6, # Height in inches
  dpi = 600,   # Resolution
  units = "in"
)

saveRDS(DNAmAL_sensitivity_plots, "../Output/whi_DNAmAL_sensitivity_results.rds")

ggsave(
  filename = "../Output/whi_methAL_sensitivity_results.pdf", 
  plot = methAL_sensitivity_plots,
  width = 10,  # Width in inches
  height = 10, # Height in inches
  dpi = 600,   # Resolution
  units = "in"
)

saveRDS(methAL_sensitivity_plots, "../Output/whi_methAL_sensitivity_results.rds")

ggsave(
  filename = "../Output/whi_methALT_sensitivity_results.pdf", 
  plot = methALT_sensitivity_plots,
  width = 10,  # Width in inches
  height = 10, # Height in inches
  dpi = 600,   # Resolution
  units = "in"
)

saveRDS(methALT_sensitivity_plots, "../Output/whi_methALT_sensitivity_results.rds")

