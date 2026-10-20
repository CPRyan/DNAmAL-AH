model_results  <-sle_results %>%
  filter(sle == "SLE score 2") %>%
  mutate(
    outcome = factor(outcome, levels = rev(c(
      "DNAmAL-AH", "methAL", "methALT",
      "DunedinPACE", "PCGrimAge", "PCPhenoAge"
    ))))


library(tidyverse)

# Replace model_results with the name of your results data frame
clock_order <- c(
  "DNAmAL-AH",
  "DunedinPACE",
  "PCGrimAge",
  "PCPhenoAge",
  "methAL",
  "methALT"
)

adjust_order <- c(
  "None",
  "DNAmAL-AH",
  "DunedinPACE",
  "PCGrimAge",
  "PCPhenoAge",
  "methAL",
  "methALT"
)

heatmap_dat <- model_results %>%
  filter(
    sle == "SLE score 2",
    term == "sle_score2_scaled"
  ) %>%
  mutate(
    outcome = factor(outcome, levels = rev(clock_order)),
    adjust_clock = factor(adjust_clock, levels = adjust_order),
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    cell_label = paste0(
      sprintf("%.3f", estimate),
      significance
    )
  ) %>%
  complete(
    outcome = factor(rev(clock_order), levels = rev(clock_order)),
    adjust_clock = factor(adjust_order, levels = adjust_order)
  )




estimate_limit <- max(abs(heatmap_dat$estimate), na.rm = TRUE)

ggplot(
  heatmap_dat,
  aes(
    x = adjust_clock,
    y = outcome,
    fill = estimate
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = cell_label),
    size = 3.5,
    na.rm = TRUE
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-estimate_limit, estimate_limit),
    na.value = "grey90",
    name = expression(beta)
  ) +
  labs(
    x = "Additional epigenetic clock included in model",
    y = "Outcome clock",
    caption = paste(
      "Cells show the SLE coefficient.",
      "*P < 0.05; **P < 0.01; ***P < 0.001.",
      "Grey diagonal cells represent self-adjustment models that were not estimated."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    plot.caption = element_text(
      hjust = 0,
      size = 9
    )
  )


ggsave("Output/HRS/Figures/HRS_models_adjusted_heatmap.pdf", height = 7, width = 7.5)
