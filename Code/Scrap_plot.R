pacman::p_load(tidyverse, patchwork)


fhs <- read_csv("Output/Validation/DNAmAL_ResultsFeb2026/fhs_survival_results.csv") %>% filter(clock == "AL_meanResid_scaled") %>% 
  mutate(clock = "DNAmAL", 
         outcome = fct_recode(outcome, Death = 'death_status', CHD = 'chd_status', CHF = 'chf_status', CVD = 'cvd_status', Stroke = 'stroke_status', Dementia = 'dementia_status'), 
         outcome = fct_relevel(outcome, c('CHD', 'CHF', 'CVD', 'Stroke', 'Dementia', 'Death'))) %>% 
  select(-`...1`)

hrs <- read_csv("Output/Validation/DNAmAL_ResultsFeb2026/HRS_outcome_results.csv") %>% filter(clock == 'DNAmALResid_scaled') %>% 
  mutate(clock = "DNAmAL", 
         Dataset = "HRS", 
         outcome = fct_recode(outcome, 'Chronic\nDisease' = "Chronic Disease"),
         outcome = fct_relevel(outcome, c('ADLs', 'IADLs', 'CIND', 'Dementia', 'Chronic\nDisease', 'Death'))) %>% 
  select(-`...1`) 

whi <- read_csv("Output/Validation/DNAmAL_ResultsFeb2026/whi_survival_results.csv") %>% filter(clock == "AL_meanResid_scaled") %>%
  filter(outcome != 'ANYFX') %>% 
  mutate(clock = "DNAmAL", 
         outcome = fct_recode(outcome, Stroke = 'STROKE', `Any\nCancer` = "ANYCANCER", Death = "DEATH"), 
         outcome = fct_relevel(outcome, c("CHD", "CHF", "Stroke", "Any\nCancer", "Death"))) %>% 
  select(-`...1`, -std_error)

df <-bind_rows(bind_rows(fhs, hrs), whi)


p1 <-fhs %>% 
  ggplot(., aes(x = outcome, y = hazard_ratio))+
  geom_point(aes(y = hazard_ratio), size = 2.4)+
  geom_errorbar(color = 'navy', aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, linewidth = 1.05)+
  ylim(c(0.95, 2.3))+
  facet_wrap(~Dataset, nrow = 1)+
  geom_hline(yintercept = 1, linetype = 'dashed')+
  theme_bw()+
  ylab("Hazard Ratio (95% CI)")+
  xlab("Outcome")+
  theme( axis.title.x = element_blank())
  

p2 <-hrs %>% 
ggplot(., aes(x = outcome, y = hazard_ratio))+
  geom_point(aes(y = hazard_ratio), size = 2.4)+
  geom_errorbar(color = 'navy', aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, linewidth = 1.05)+
  ylim(c(0.95, 2.3))+
  facet_wrap(~Dataset, nrow = 1)+
  geom_hline(yintercept = 1, linetype = 'dashed')+
  theme_bw()+
  ylab("Hazard Ratio (95% CI)")+
  xlab("Outcome")+
  theme(axis.title.y = element_blank())

p3 <-whi %>% 
ggplot(., aes(x = outcome, y = hazard_ratio))+
  geom_point(aes(y = hazard_ratio), size = 2.4)+
  geom_errorbar(color = 'navy', aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, linewidth = 1.05)+
  ylim(c(0.95, 2.3))+
  facet_wrap(~Dataset, nrow = 1)+
  geom_hline(yintercept = 1, linetype = 'dashed')+
  theme_bw()+
  ylab("Hazard Ratio (95% CI)")+
  xlab("Outcome")+
  theme( axis.title.x = element_blank(), axis.title.y = element_blank())

p1 + p2 + p3

ggsave(here::here("Output/Figures", "Hazard_ratios_DNAmAL_only.pdf"), width = 12, height = 3.5)

