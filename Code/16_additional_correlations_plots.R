####################
# Set up library and cores:
####################
.libPaths("/ifs/sec/cpc/addhealth/users/cpr/R/x86_64-pc-linux-gnu-library/4.4")
# Where to find my custom package library
library(haven)
library(tidyverse)
library(janitor)
library(openxlsx)
library(gtsummary)
library(broom)




###############
# Load data
###############

pheno <-read_csv(here::here("Output/Data", "AL_phenos_with_predictions.csv"))%>% 
  sjlabelled::as_character(aid) %>%
  # rename something suitable
  rename(DNAmAL = "AL_mean_al_msa_probes_pred", 
         AL = "AL_mean")


# Nobody in this dataset is missing AL or DNAmAL, so we don't need to clean it. 

# Scale
pheno_scaled <-pheno %>%
  mutate(across(
    .cols = c(sbp, dbp, pr, tc, bmi, waist, hba1c_log, tg_log, crp_log, il6_log, tnf_log, hdl_inv, AL, DNAmAL), 
    .fns = ~as.numeric(scale(.)),
    .names = "{.col}"
  ))


biomarkers <-c("sbp", "dbp", "pr", "tc", "bmi", "waist", "hba1c_log", "tg_log", "crp_log", "il6_log", "tnf_log", "hdl_inv")
scores <-c("AL", "DNAmAL")

# r2_results <-crossing(biomarker = biomarkers, score = scores) %>%
#   pmap_dfr(function(biomarker, score) {
#     x <-pheno_scaled[[biomarker]]
#     y <-pheno_scaled[[score]]
#     keep <-complete.cases(x, y)
#     r <-cor(x[keep], y[keep])
#     tibble(biomarker, score, r, R2 = r^2, N = sum(keep))
#   })

lm_results <-crossing(biomarker = biomarkers, score = scores) %>%
  pmap_dfr(function(biomarker, score) {
    fit <-lm(reformulate(score, biomarker), data = pheno_scaled)
    tidy(fit, conf.int = TRUE) %>%
      filter(term == score) %>%
      transmute(biomarker, score, estimate, 
                CI_lower = conf.low, 
                CI_upper = conf.high, 
                N = nobs(fit))
  }) %>%
  sjlabelled::as_factor(score) %>%
  mutate(score = fct_relevel(score, "DNAmAL"))



ggplot(lm_results, aes(x = score, y = estimate, fill = score))+
  geom_col(width =0.7)+
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width =0.2)+
  facet_wrap(~biomarker, 
             labeller = as_labeller(c("bmi" = "BMI", 
                                      "crp_log" = "CRP (log)", 
                                      "dbp" = "DBP", 
                                      "hba1c_log" = "HbA1c (log)", 
                                      "hdl_inv" = "HDL (inv)", 
                                      "il6_log" = "IL6 (log)", 
                                      "pr" = "PR", 
                                      "sbp" = "SBP", 
                                      "tc" = "Tcells", 
                                      "tg_log" = "Trigs (log)", 
                                      "tnf_log" = "TNFa (log)", 
                                      "waist" = "Waist Circ.")))+
  scale_fill_manual(values = c("AL" = "goldenrod3", "DNAmAL" = "dodgerblue4"))+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "top")+
  ylab("Estimate (95% CI)")+xlab(NULL)
  
ggsave(here::here("Output/Figures", "barchart_items_AL_and_DNAmAL_facetted_horiz.pdf"), height = 6, width=7)
 


