pacman::p_load(readr, tidvyerse, readxl)


hypert_table <- read_excel("Output/Tables/Addhealth_Tables.xlsx", 
                              sheet = "Hypertension Cross Tab")

diab_table <- read_excel("Output/Tables/Addhealth_Tables.xlsx",
                             sheet = "Diabetes Cross Tab")

library(tidyverse)
library(scales)

prep_transition <- function(tab, outcome_name) {
  
  tab %>%
    # Assumes columns are ordered: Wave 5 status, W6 No, W6 Yes, W6 Missing
    setNames(c("w5", "No", "Yes", "Missing")) %>%
    mutate(
      outcome = outcome_name,
      w5 = recode(as.character(w5),
                  `0` = "No",
                  `1` = "Yes")
    ) %>%
    pivot_longer(
      cols = c(No, Yes, Missing),
      names_to = "w6",
      values_to = "n"
    ) %>%
    group_by(outcome, w5) %>%
    mutate(
      total = sum(n),
      prop = n / total
    ) %>%
    ungroup()
}


transition_data <- bind_rows(
  prep_transition(hypert_table, "Hypertension"),
  prep_transition(diab_table, "Type 2 diabetes")
) %>%
  mutate(
    outcome = factor(
      outcome,
      levels = c("Hypertension", "Type 2 diabetes")
    ),
    w5 = factor(w5, levels = c("No", "Yes")),
    w6 = factor(w6, levels = c("No", "Yes", "Missing")),
    text_color = if_else(w6 == "No", "black", "white")
  )


group_totals <- transition_data %>%
  distinct(outcome, w5, total)


transition_plot <- ggplot(
  transition_data,
  aes(x = w5, y = prop, fill = w6)
) +
  geom_col(
    width = 0.70,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    aes(
      label = percent(prop, accuracy = 0.1),
      color = text_color
    ),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    size = 3.5
  ) +
  geom_text(
    data = group_totals,
    aes(
      x = w5,
      y = 1.02,
      label = paste0("n = ", comma(total))
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.3
  ) +
  facet_wrap(~outcome, ncol = 1) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "No" = "grey80",
      "Yes" = "#38678D",
      "Missing" = "grey45"
    ),
    name = "Wave 6 status"
  ) +
  scale_color_identity() +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = percent_format(),
    limits = c(0, 1.13),
    expand = c(0, 0)
  ) +
  labs(
    x = "Wave 5 status",
    y = "Percentage of Wave 5 group"
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.spacing = unit(1, "lines")
  )

transition_plot

ggsave("Output/Figures/health_transition_percentages.pdf")
