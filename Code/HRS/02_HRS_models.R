library(tidyverse)
library(broom)
library(openxlsx)


# Is the association between stressful life events and DNAmAL robust to adjustment for established epigenetic aging measures, and conversely, do those measures retain independent associations with stress after adjustment for DNAmAL?

# Model A: Does DNAmAL contain stress information beyond e.g. GrimAge?
  
# Model B: Does e.g. GrimAge contain stress information beyond DNAmAL?


# 1. Define the variables
covars <- c("Female", "PAGE", "smoker1213")

sle_vars <- c(
  "sle_score1_scaled",
  "sle_score2_scaled",
  "sle_score3_scaled"
)

clock_vars <- c(
  "DNAmAL_r_scaled",
  "DunedinPACE_r_scaled",
  "PCGrimAge_r_scaled",
  "PCPhenoAge_r_scaled", 
  "methAL_r_scaled",
  "methALT_r_scaled"
)


# 2. Create every model I want. 

# This creates the base model plus models adjusting for each of the other clocks.
model_grid <- crossing(
  outcome = clock_vars,
  sle = sle_vars,
  adjust_clock = c("None", clock_vars)
) %>%
  filter(adjust_clock == "None" | adjust_clock != outcome)

# 3. Fit all models. 
fit_clock_model <- function(outcome, sle, adjust_clock) {
  rhs <- c(covars, sle, if (adjust_clock != "None") adjust_clock)
  lm(reformulate(rhs, response = outcome), data = hrs)
}

models <- model_grid %>%
  mutate(
    model = pmap(list(outcome, sle, adjust_clock), fit_clock_model),
    n = map_int(model, nobs),
    r2 = map_dbl(model, ~ summary(.x)$r.squared)
  )

# 4. Make a tidy table.
results <- models %>%
  mutate(tidy = map(model, ~ tidy(.x, conf.int = TRUE))) %>%
  select(outcome, sle, adjust_clock, n, r2, tidy) %>%
  unnest(tidy)

# 5. Make names readable. 
results <- results %>%
  mutate(
    outcome = recode(outcome,
                     DNAmAL_r_scaled = "DNAmAL-AH",
                     DunedinPACE_r_scaled = "DunedinPACE",
                     PCGrimAge_r_scaled = "PCGrimAge",
                     PCPhenoAge_r_scaled = "PCPhenoAge", 
                     methAL_r_scaled = "methAL",
                     methALT_r_scaled = "methALT"
    ),
    sle = recode(sle,
                 sle_score1_scaled = "SLE score 1",
                 sle_score2_scaled = "SLE score 2",
                 sle_score3_scaled = "SLE score 3"
    ),
    adjust_clock = recode(adjust_clock,
                          DNAmAL_r_scaled = "DNAmAL-AH",
                          DunedinPACE_r_scaled = "DunedinPACE",
                          PCGrimAge_r_scaled = "PCGrimAge",
                          PCPhenoAge_r_scaled = "PCPhenoAge", 
                          methAL_r_scaled = "methAL",
                          methALT_r_scaled = "methALT"
    )
  )




sle_results <- results %>%
  filter(
    (sle == "SLE score 1" & term == "sle_score1_scaled") |
      (sle == "SLE score 2" & term == "sle_score2_scaled") |
      (sle == "SLE score 3" & term == "sle_score3_scaled")
  )


sle_score1 <- sle_results %>%
  filter(
    sle == "SLE score 1"  
    )

sle_score2 <- sle_results %>%
  filter(
    sle == "SLE score 2"  
  )

sle_score3 <- sle_results %>%
  filter(
    sle == "SLE score 3"  
  )



m1 <-models %>% 
  filter(sle == "sle_score2_scaled" & outcome == "DNAmAL_r_scaled" & adjust_clock == "None") %>% 
  pull(model)

tab_model(m1)


# Find a way to add R2 for each model, then save as one large file for SLE 2. 

wb <-openxlsx::loadWorkbook("Output/HRS/Tables/SLE_DNAm_model_tables.xlsx")

sheets <- list(
  sle_score1_clocks = sle_score1,
  sle_score2_clocks = sle_score2,
  sle_score3_clocks = sle_score3
)

for (s in names(sheets)) {
  if (s %in% names(wb)) openxlsx::removeWorksheet(wb, s)
  openxlsx::addWorksheet(wb, s)
  openxlsx::writeData(wb, s, sheets[[s]])
}

openxlsx::saveWorkbook(
  wb,
  "Output/HRS/Tables/SLE_DNAm_model_tables.xlsx",
  overwrite = TRUE
)




# SLE score 2 results only
plot_dat <- sle_results %>%
  filter(sle == "SLE score 2") %>%
  mutate(
    outcome = factor(outcome, levels = rev(c(
      "DNAmAL-AH", "DunedinPACE", "PCGrimAge", "PCPhenoAge", "methAL", "methALT"
    ))),
    adjust_clock = recode(adjust_clock, "None" = "Base model (DNAmAL-AH)")
  )





# A. Association of SLE with each biomarker
p1 <- plot_dat %>%
  filter(adjust_clock == "Base model (DNAmAL-AH)") %>%
  ggplot(aes(x = estimate, y = outcome, color = outcome)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15, linewidth = 1) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("#5773C0FF", "#7DA7EAFF", "#8CC8BCFF", "#EDB144FF", "#C8570DFF","#B83326FF")) +
  labs(
    x = "Assoc. w/ Stressful Life Events (standardized β)",
    y = "Epigenetic measure assoc. with SLE") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )+
  theme(legend.position = "none")


  

p1

ggsave(filename = "Output/HRS/Figures/HRS_SLE_predict_clocks.pdf",  device = cairo_pdf, width = 5, height =4)












adjust_cols <- c(
  "Base model\n (DNAmAL-AH only)" = "#B83326FF",
  "DunedinPACE"        = "#C64A3EFF",
  "PCGrimAge"       = "#D46157FF",
  "PCPhenoAge"   = "#DF7A70FF",
  "methAL"     = "#EAA39CFF",
  "methALT"    = "#F3CCC8FF"
)




clock_order <- c(
  "methALT", 
  "methAL",
  "PCPhenoAge", 
  "PCGrimAge", 
  "DunedinPACE", 
  "Base model (DNAmAL-AH)")

p2_dat <- plot_dat %>%
  filter(outcome == "DNAmAL-AH") %>%
  mutate(
    adjust_clock = fct_relevel(adjust_clock, clock_order), 
    adjust_clock = recode(adjust_clock, "Base model (DNAmAL-AH)" = "Base model\n (DNAmAL-AH only)")
  )


p2 <- ggplot(p2_dat, aes(x = estimate, y = adjust_clock, colour = adjust_clock)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15, linewidth =1) +
  geom_point(size = 3.5) +
  labs(
    x = "Assoc. of SLE with DNAmAL-AH (standardized β)",
    y = "Additional epigenetic measure in Base Model"  ) +
  scale_colour_manual(values = adjust_cols) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none", 
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )


ggsave(filename = "Output/HRS/Figures/HRS_SLE_DNAmAL_after_other_biomarkers.pdf",  device = cairo_pdf)


p_1and2 <-p1 + p2 +
  plot_layout(widths = c(1, 1)) &
  theme(plot.title = element_text(face = "bold", size = 12))


ggsave(plot = p_1and2, filename =  "Output/HRS/Figures/HRS_SLE_predict_clocks_and_DNAmAL_after_other_biomarkers.pdf", width = 10, height = 5,  device = cairo_pdf)

