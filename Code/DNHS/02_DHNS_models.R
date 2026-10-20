# GEE

# 1. Variables
covars <- c("Female", "Age", "focal_smoking")

stress_vars <- c(
  "focal_trauma",
  "focal_stress"
)

clock_vars <- c(
  "DNAmAL_r_scaled",
  "methAL_r_scaled",
  "methALT_r_scaled",
  "DunedinPACE_r_scaled",
  "PCGrimAge_r_scaled",
  "PCPhenoAge_r_scaled"
)

clock_labels <- c(
  DNAmAL_r_scaled = "DNAmAL",
  methAL_r_scaled = "methAL",
  methALT_r_scaled = "methALT",
  DunedinPACE_r_scaled = "DunedinPACE",
  PCGrimAge_r_scaled = "PCGrimAge",
  PCPhenoAge_r_scaled = "PCPhenoAge"
)

stress_labels <- c(
  focal_trauma = "Trauma Score",
  focal_stress = "Stress Score"
)



# 2. Model grid
model_grid <- crossing(
  outcome = clock_vars,
  stress = stress_vars,
  adjust_clock = c("None", clock_vars)
) %>%
  filter(adjust_clock == "None" | adjust_clock != outcome)

# 3. Fit GEE models
fit_clock_model <- function(outcome, stress, adjust_clock) {
  rhs <- c(covars, stress, if (adjust_clock != "None") adjust_clock)
  
  dat <- df_focal %>%
    select(RESP, all_of(c(outcome, rhs))) %>%
    drop_na() %>%
    arrange(RESP)
  
  geeglm(
    reformulate(rhs, response = outcome),
    id = RESP,
    data = dat,
    family = gaussian,
    corstr = "exchangeable"
  )
}

models <- model_grid %>%
  mutate(
    model = pmap(list(outcome, stress, adjust_clock), fit_clock_model),
    n_obs = map_int(model, nobs),
    n_people = map_int(model, ~ n_distinct(.x$id))
  )


# 4. Tidy results with 95% CI
results <- models %>%
  mutate(tidy = map(model, ~ broom::tidy(.x) %>%
                      mutate(
                        conf.low = estimate - 1.96 * std.error,
                        conf.high = estimate + 1.96 * std.error
                      ))) %>%
  select(outcome, stress, adjust_clock, n_obs, n_people, tidy) %>%
  unnest(tidy)

# Keep the coefficient for the focal stress/trauma measure
stress_results <- results %>%
  filter(term == stress) %>%
  mutate(
    outcome = recode(outcome, !!!clock_labels),
    stress = recode(stress, !!!stress_labels),
    adjust_clock = recode(adjust_clock, !!!clock_labels)
  )

trauma_results <- stress_results %>% filter(stress == "Trauma Score")
focal_stress_results <- stress_results %>% filter(stress == "Stress Score")





dir.create("Output/DNHS/Tables", recursive = TRUE, showWarnings = FALSE)

wb <- createWorkbook()

addWorksheet(wb, "Trauma")
writeData(wb, "Trauma", trauma_results)

addWorksheet(wb, "Stress")
writeData(wb, "Stress", focal_stress_results)

saveWorkbook(
  wb,
  "Output/DNHS/Tables/DNHS_stress_DNAm_GEE_models.xlsx",
  overwrite = TRUE
)








plot_dat <- stress_results %>%
  mutate(
    outcome = factor(outcome, levels = rev(c(
      "DNAmAL", "methAL", "methALT",
      "DunedinPACE", "PCGrimAge", "PCPhenoAge"
    ))),
    adjust_clock = recode(adjust_clock, "None" = "Base model (DNAmAL only)")
  )

clock_cols <- c(
  "DNAmAL" = "#8C2D04",
  "methAL" = "#C75B12",
  "methALT" = "#E68613",
  "DunedinPACE" = "#F3A738",
  "PCGrimAge" = "#F7C95E",
  "PCPhenoAge" = "#FCE39C"
)

p1 <- plot_dat %>%
  filter(adjust_clock == "Base model (DNAmAL only)") %>%
  ggplot(aes(x = estimate, y = outcome, color = outcome)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15) +
  geom_point(size = 3) +
  facet_wrap(~stress) +
  scale_color_manual(values = clock_cols) +
  labs(
    x = "Association with epigenetic biomarker",
    y = NULL,
    title = "A. Association of stress exposures with epigenetic biomarkers"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )





adjust_cols <- c(
  "Base model (DNAmAL only)" = "#8C2D04",
  "methAL" = "#1B9E77",
  "methALT" = "#66A61E",
  "DunedinPACE" = "#7570B3",
  "PCGrimAge" = "#E7298A",
  "PCPhenoAge" = "#1F78B4"
)

p2_dat <- plot_dat %>%
  filter(outcome == "DNAmAL") %>%
  mutate(adjust_clock = factor(
    adjust_clock,
    levels = rev(c(
      "Base model (DNAmAL only)", "methAL", "methALT",
      "DunedinPACE", "PCGrimAge", "PCPhenoAge"
    ))
  ))

p2 <- ggplot(p2_dat, aes(x = estimate, y = adjust_clock, colour = adjust_clock)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15) +
  geom_point(size = 3) +
  facet_wrap(~stress) +
  scale_colour_manual(values = adjust_cols) +
  labs(
    x = "Stress-exposure coefficient predicting DNAmAL",
    y = NULL,
    title = "B. DNAmAL association after adjustment for other biomarkers"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )



p_1and2 <- p1 / p2 &
  theme(plot.title = element_text(face = "bold", size = 12))

dir.create("Output/DNHS/Figures", recursive = TRUE, showWarnings = FALSE)

ggsave(
  "Output/DNHS/Figures/DNHS_stress_DNAm_GEE_models.pdf",
  p_1and2, width = 10, height = 9, device = cairo_pdf
)
