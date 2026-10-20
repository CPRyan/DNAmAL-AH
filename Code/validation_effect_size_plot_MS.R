pacman::p_load(tidyverse)

# ---------------------------
# Read + clean
# ---------------------------
fhs <- read_csv("Output/Validation/DNAmAL_ResultsFeb2026/fhs_survival_results.csv") %>% 
  filter(clock == "AL_meanResid_scaled") %>% 
  mutate(
    clock   = "DNAmAL",
    Dataset = "FHS",
    outcome = fct_recode(
      outcome,
      Death    = "death_status",
      CHD      = "chd_status",
      CHF      = "chf_status",
      CVD      = "cvd_status",
      Stroke   = "stroke_status",
      Dementia = "dementia_status"
    )
  ) %>% 
  select(-`...1`)

hrs <- read_csv("Output/Validation/DNAmAL_ResultsFeb2026/HRS_outcome_results.csv") %>% 
  filter(clock == "DNAmALResid_scaled") %>% 
  mutate(
    clock   = "DNAmAL",
    Dataset = "HRS",
    outcome = fct_recode(
      outcome,
      "Chronic\nDisease" = "Chronic Disease"
    )
  ) %>% 
  select(-`...1`) 

whi <- read_csv("Output/Validation/DNAmAL_ResultsFeb2026/whi_survival_results.csv") %>% 
  filter(clock == "AL_meanResid_scaled") %>%
  filter(outcome != "ANYFX") %>% 
  mutate(
    clock   = "DNAmAL",
    Dataset = "WHI",
    outcome = fct_recode(
      outcome,
      Stroke       = "STROKE",
      "Any\nCancer" = "ANYCANCER",
      Death        = "DEATH"
    )
  ) %>% 
  select(-`...1`, -std_error)

df <- bind_rows(fhs, whi, hrs)

# ---------------------------
# Shared ordering + domains
# Top to bottom:
# CVD-related -> cognitive/functional -> other -> death
# ---------------------------
outcome_order_top_to_bottom <- c(
  "CHD",
  "CHF",
  "CVD",
  "Stroke",
  "CIND",
  "Dementia",
  "ADLs",
  "IADLs",
  "Chronic\nDisease",
  "Any\nCancer",
  "Death"
)

df <- df %>%
  mutate(
    Dataset = factor(Dataset, levels = c("FHS", "WHI", "HRS")),
    Domain = case_when(
      outcome %in% c("CHD", "CHF", "CVD", "Stroke") ~ "Cardiovascular",
      outcome %in% c("CIND", "Dementia", "ADLs", "IADLs") ~ "Cognitive / functional",
      outcome %in% c("Chronic\nDisease") ~ "Chronic disease",
      outcome %in% c("Any\nCancer") ~ "Cancer",
      outcome %in% c("Death") ~ "Mortality",
      TRUE ~ "Other"
    ),
    # reverse because ggplot places first factor level at the bottom on y-axis
    outcome = factor(outcome, levels = rev(outcome_order_top_to_bottom))
  )

# ---------------------------
# Colors
# ---------------------------
domain_cols <- c(
  "Cardiovascular"         = "#2C7FB8",  # blue
  "Cognitive / functional" = "#8E5EA2",  # purple
  "Chronic disease"        = "#C47A2C",  # warm rust/gold-brown
  "Cancer"                 = "#1B9E77",  # green/teal
  "Mortality"              = "#333333"   # dark gray
)

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(df, aes(x = hazard_ratio, y = outcome, color = Domain)) +
  geom_vline(
    xintercept = 1,
    color = "grey55",
    linewidth = 0.7,
    linetype = "22"
  ) +
  geom_segment(
    aes(x = lower_ci, xend = upper_ci, y = outcome, yend = outcome),
    linewidth = 1.05,
    alpha = 0.75,
    lineend = "round"
  ) +
  geom_point(size = 3.8) +
  facet_wrap(~Dataset, nrow = 1, scales = "free_y") +
  scale_color_manual(values = domain_cols) +
  scale_x_continuous(
    breaks = seq(1.0, 2.2, by = 0.2),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  coord_cartesian(xlim = c(0.95, 2.3)) +
  labs(
    x = "Hazard Ratio (95% CI)",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 8)),
    axis.text.x  = element_text(size = 11, color = "grey20"),
    axis.text.y  = element_text(size = 11, color = "grey20"),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.5),
    
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    
    panel.spacing = grid::unit(1.2, "lines"),
    plot.margin = margin(8, 12, 8, 8)
  )

p

ggsave(
  here::here("Output/Figures", "Hazard_ratios_DNAmAL_only.pdf"),
  plot = p,
  width = 11.5,
  height = 4.6
)
