# Load packages
library(readr)
library(dplyr)
library(ggplot2)
library(ggforce)
library(patchwork)
library(gridExtra)

# Load datasets
minimal_results <- read_csv("../Output/HRS_outcome_results.csv")

##################################################################
# Plotting details
##################################################################

# Plot function
outcomes <- c("Death",
              "Dementia",
              "ADLs",
              "IADLs",
              "CIND",
              "Chronic Disease")

outcome_names <- c("Death",
                   "Dementia",
                   "ADLs",
                   "IADLs",
                   "CIND",
                   "Chronic Disease")

## Create a named vector for easy lookup
outcome_labels <- setNames(outcome_names, outcomes)

# Define all available clocks with their colors
all_clock_colors <- c(
  "DNAmALResid_scaled" = "navy", 
  "methALResid_scaled" = "blue", 
  "methALTResid_scaled" = "#4c8db3", 
  "DunedinPACEResid_scaled" = "#a13cc9", 
  "PCPhenoAgeResid_scaled" = "#f37577", 
  "PCGrimAgeResid_scaled" = "#b31012"
)

all_clock_labels <- c(
  "DNAmALResid_scaled" = "DNAmAL",
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
    title = "Effect Sizes by Outcome and Clock in HRS"
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
  
  # Hardcoded outcome order
  sorted_outcomes <- c("Death",
                       "CIND",
                       "Dementia",
                       "ADLs",
                       "IADLs",
                       "Chronic Disease")
  
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
      # Order outcomes based on hardcoded order
      outcome = factor(outcome, levels = sorted_outcomes),
      outcome_name = factor(outcome_name, levels = outcome_labels[sorted_outcomes]),
      # Make clocks a factor for consistent ordering based on selected_clocks order
      clock_name = factor(clock_name, levels = unname(clock_labels[selected_clocks]))
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
  
  # Create the common plot structure with NO y-axis label
  p_base <- ggplot(plot_data) +
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
    # Remove y-axis label from individual plots
    labs(x = NULL, y = NULL) +
    # Common theme elements
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(size = 16, face = "bold"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80", linewidth = 0.5),
      panel.spacing = unit(1, "lines"),
      # Explicitly remove y-axis title
      axis.title.y = element_blank()
    )
  
  # Split plot data by outcome position
  first_row_outcomes <- sorted_outcomes[1:3]
  second_row_outcomes <- sorted_outcomes[4:6]
  
  # First row plots (with original scale and NO x-axis text)
  first_row_data <- plot_data %>% filter(outcome %in% first_row_outcomes)
  
  first_row_plot <- p_base %+% first_row_data +
    scale_y_continuous(limits = c(0.8, 2.1)) +
    facet_wrap(~ panel_label, nrow = 1, ncol = 3) +
    # Add the title to the top plot only
    labs(title = title) +
    theme(
      # Hide x-axis text and ticks for top row
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      # Title styling
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24, margin = margin(b = 15)),
      # Y-axis styling
      axis.text.y = element_text(size = 9)
    )
  
  # Second row plots (with different scale and WITH x-axis text)
  second_row_data <- plot_data %>% filter(outcome %in% second_row_outcomes)
  
  second_row_plot <- p_base %+% second_row_data +
    scale_y_continuous(limits = c(0.85, 1.5)) + # Custom scale for second row
    facet_wrap(~ panel_label, nrow = 1, ncol = 3) +
    theme(
      # Show x-axis text for bottom row
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      # Y-axis styling
      axis.text.y = element_text(size = 9)
    )
  
  # Combine using patchwork and add a common y-axis label
  combined_plot <- first_row_plot / second_row_plot +
    plot_layout(heights = c(1, 1.2)) +  # Slightly more height for bottom row with labels
    # Add the common y-axis label with patchwork's annotation system
    plot_annotation(
      title = NULL,
      caption = NULL,
      theme = theme(plot.margin = margin(5.5, 5.5, 5.5, 35))  # Add extra left margin
    ) &
    # This adds a common y-axis label to all plots
    theme(
      plot.tag = element_text(size = 0)  # Hide tags if any
    )
  
  # Return the combined plot with an invisible label using draw_label
  # We wrap this in a call to a function that adds the label
  add_shared_label <- function(plot) {
    # Create a viewport for the y-axis label
    vp <- viewport(x = 0.02, y = 0.5, just = c("left", "center"))
    
    # Return a function that, when called, will draw the plot and then add the label
    function() {
      # Draw the plot
      print(plot)
      
      # Push the viewport
      pushViewport(vp)
      
      # Add the y-axis label
      grid.text("Effect Size (95% CI)", rot = 90, gp = gpar(fontsize = 18))
      
      # Pop the viewport
      popViewport()
    }
  }
  
  # This is the simplest way to create the common y-axis label with minimum dependencies
  return(combined_plot + 
           theme(plot.margin = margin(5.5, 5.5, 5.5, 50)) +  # Make room for y-axis label
           labs(tag = "Effect Size (95% CI)") +  # Use tag for y-axis label
           theme(plot.tag = element_text(size = 22, angle = 90, vjust = 0.5, hjust = 0.5),
                 plot.tag.position = c(-0.07, 1.0)))  # Position in the left margin, vertically centered
}

##################################################################
# Make plots
##################################################################

# Create plot with all clocks (default behavior)
all_outcomes_plot <- plot_all_survival_results_dot(minimal_results)
print(all_outcomes_plot)

##################################################################
# Save plots
##################################################################

ggsave(
  filename = "../Output/HRS_overall_outcome_results.pdf", 
  plot = all_outcomes_plot,
  width = 12,  # Width in inches
  height = 8, # Height in inches
  dpi = 600,   # Resolution
  units = "in"
)

saveRDS(all_outcomes_plot, "../Output/HRS_overall_outcome_results.rds")

