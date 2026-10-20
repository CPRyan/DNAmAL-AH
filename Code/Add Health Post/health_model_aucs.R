pacman::p_load(readr, tidyverse, readxl, sjlabelled)

My_Theme = theme(
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 16),
  axis.title.y = element_text(size = 18), 
  axis.text.y = element_text(size = 14))


health_aucs <- read_excel("Output/Tables/Addhealth_Tables.xlsx", 
                               sheet = "auc_health_models")



clock_labels <- c(
  "AL_r_scaled"          = "Allostatic\n Load",
  "DNAmAL_r_scaled"      = "DNAmAL-AH",
  "DunedinPACE_r_scaled" = "DunedinPACE",
  "PCGrimAge_r_scaled"   = "PC GrimAge",
  "PCPhenoAge_r_scaled"   = "PC PhenoAge", 
  "methAL_r_scaled"      = "methAL",
  "methALT_r_scaled"     = "methALT"
)


outcome_labels <- c(
  "diab_joint_w6_home"   = "Type II Diabetes (n=2,762; events=399)",
  "hypert_joint_w6_home" = "Hypertension (n=2,762; events=1034)"
)

# Optional but recommended: make the ordering and palette matching explicit
health_aucs_plot <- health_aucs |>
  mutate(
    model = factor(model, levels = c("Model 1", "Model 2", "Model 3")),
    clock = factor(
      clock,
      levels = c(
        "AL_r_scaled",
        "DNAmAL_r_scaled",
        "DunedinPACE_r_scaled",
        "PCGrimAge_r_scaled",
        "PCPhenoAge_r_scaled", 
        "methAL_r_scaled",
        "methALT_r_scaled"
      )
    ),
    outcome = factor(
      outcome,
      levels = c(
        "diab_joint_w6_home",
        "hypert_joint_w6_home"
      )
    )
  )

pd_auc <- position_dodge(width = 0.80)

ggplot(
  health_aucs_plot,
  aes(
    x = model,
    y = AUC,
    fill = clock,
    group = clock
  )
) +
  geom_col(
    position = pd_auc,
    width = 0.70,
    alpha = 0.85,
    colour = "white",
    linewidth = 0.25
  ) +
  facet_wrap(
    ~ outcome,
    ncol = 1,
    labeller = labeller(outcome = outcome_labels)
  ) +
  paletteer::scale_fill_paletteer_d("MetBrewer::Nizami", labels = clock_labels) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    x = "Model x Clock",
    y = "Area under the ROC curve (AUC)",
    fill = "Clock"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    legend.position = "right",
    legend.background = element_rect(
      fill = scales::alpha("white", 0.85),
      colour = "grey50"
    ),
    legend.key = element_blank()
  ) +
  My_Theme


ggsave("Output/Figures/health_model_AUCs_barchart.pdf", height = 5, width = 7)

