pacman::p_load(readr, tidyverse, readxl, sjlabelled)

My_Theme = theme(
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 16),
  axis.title.y = element_text(size = 18), 
  axis.text.y = element_text(size = 14))


m1 <- read_excel("Output/Tables/Addhealth_Tables.xlsx", 
                               sheet = "health6_m1") %>% 
  filter(term == 'clock')%>% 
  mutate(model = "m1")



m2 <- read_excel("Output/Tables/Addhealth_Tables.xlsx", 
                 sheet = "health6_m2_health5") %>% 
  filter(term == 'clock')%>% 
  mutate(model = "m2")



m3 <- read_excel("Output/Tables/Addhealth_Tables.xlsx", 
                 sheet = "health6_m3_health5_cells") %>% 
  filter(term == 'clock') %>% 
  mutate(model = "m3")



all_mods <-bind_rows(m1, m2, m3) %>% 
  as_factor(clock) %>% 
  mutate(clock = fct_relevel(clock, "AL_r_scaled", "DNAmAL_r_scaled", "DunedinPACE_r_scaled", "PCGrimAge_r_scaled", "PCPhenoAge_r_scaled", "methAL_r_scaled", "methALT_r_scaled")) 


health_tab_clocks <-all_mods %>% 
  rename(Outcome = outcome, 
         'Epigenetic Measure' = clock, 
         Estimate = estimate) %>% 
  mutate(Outcome = fct_recode(Outcome,   
                              "Type II Diabetes" = "diab_joint_w6_home",
                              "Hypertension" = "hypert_joint_w6_home"),
         `Epigenetic Measure` = fct_recode(`Epigenetic Measure`, 
                                           "Allostatic Load" = "AL_r_scaled",         
                                           "DNAmAL-AH" =  "DNAmAL_r_scaled",    
                                           "DunedinPACE" = "DunedinPACE_r_scaled", 
                                           "PCGrimAge" = "PCGrimAge_r_scaled", 
                                           "PCPhenoAge" = "PCPhenoAge_r_scaled", 
                                           "methAL" =  "methAL_r_scaled",     
                                           "methALT" = "methALT_r_scaled")) %>% 
  select( model, Outcome, `Epigenetic Measure`, Estimate, CI_lower, CI_upper,  N, Events)
  

write_csv(health_tab_clocks, "Output/Tables/add_health_tables_all_models.csv")


pd <- position_dodge(width = 0.65)


clock_labels <- c(
  "AL_r_scaled"          = "Allostatic\n Load",
  "DNAmAL_r_scaled"      = "DNAmAL-AH",
  "DunedinPACE_r_scaled" = "DunedinPACE",
  "PCGrimAge_r_scaled"   = "PCGrimAge",
  "PCPhenoAge_r_scaled"   = "PCPhenoAge", 
  "methAL_r_scaled"      = "methAL",
  "methALT_r_scaled"     = "methALT"
)

outcome_labels <- c(
  "diab_joint_w6_home"   = "Type II Diabetes",
  "hypert_joint_w6_home" = "Hypertension"
)

ggplot(
  all_mods,
  aes(
    x = clock,
    y = estimate,
    colour = model,
    group = model
  )
) +
  geom_hline(yintercept = 1 , linetype = "dashed")+
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.4,
    position = pd
  ) +
  geom_point(
    size = 4,
    position = pd,
    alpha = 0.7
  ) +
  facet_wrap(~ outcome, ncol = 1, 
             labeller = labeller(outcome = outcome_labels)) +
  scale_x_discrete(labels = clock_labels) +
  labs(
    x = "Clock",
    y = "Estimate",
    colour = "Model"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )+
  scale_color_manual(values = c("#0D1D3CFF", "#4F5F28FF", "#871C0FFF"))+
  scale_fill_manual(values = c("#0D1D3CFF", "#4F5F28FF", "#871C0FFF"))+
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.85),
      colour = "grey50"
    ),
    legend.key = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_rect(fill = "white", colour = "black")
  )+My_Theme


ggsave("Output/Figures/health_model_estimates.pdf", height = 6, width = 7)
