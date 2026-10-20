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


# Setup paths
health_home <-"/ifs/sec/cpc/addhealth/users/share/ForCalen"
epig_home <-"/ifs/sec/cpc/addhealth/epig"


# Load health data
health5 <-read_sas(paste0(health_home, "/", "w5health.sas7bdat"))
health6 <-read_sas(paste0(health_home, "/", "w6health.sas7bdat"))

# Load allostatic load predictors
pheno <-read_csv(here::here("Output/Data", "AL_phenos_with_predictions.csv"))%>% 
  sjlabelled::as_character(aid) %>%
  # rename something suitable
  rename(DNAmAL = "AL_mean_al_msa_probes_pred", 
         AL = "AL_mean")
  
  
# Load Clock Data (already has)
clocks <-read_csv(paste0(epig_home, "/clocks/clocks.nobatchcorrection.postdetp.fordissemination.092525.csv"))%>% dplyr::select(aid, DunedinPACE, PCGrimAge, PCPhenoAge) %>% sjlabelled::as_character(aid)

# Load methAL and methALT
meth <- read_csv("Output/Data/AddHealth_MethAL_ALT.csv")

############
# Clean up health measures and rename.
############
health5 <-health5 %>% 
  rename(hypert_w5_surv = 'H5ID6C', 
         hypert_age_w5_surv = 'H5ID6CA',
         diab_w5_surv = 'H5ID6D', 
         diab_age_w5_surv = 'H5ID6DA',
         diag_depress_w5_surv = 'H5ID6G', 
         diag_depress_age_w5_surv = 'H5ID6GA', 
         depress_scale_w5_surv = cesd,
         diag_hbp_w5_home = 'H5Q045A', 
         diag_diab_w5_home = 'H5Q045B', 
         hypert_joint_w5_home = 'H5HTENJC', 
         diab_joint_w5_home = 'H5DIABJC'
         ) %>%
  mutate(across(where(is.numeric), ~na_if(.x, c(997))),
         across(where(is.numeric), ~na_if(.x, c(94))),
         across(where(is.numeric), ~na_if(.x, c(95))),
         across(where(is.numeric), ~na_if(.x, c(96))))


health6 <-health6 %>% 
  rename(diag_diab_w6_surv = 'H6ID7D', 
         diag_diab_age_w6_surv = 'H6ID7DA',
         curr_diab_w6_surv = 'H6ID7DC', 
         diag_hypert_w6_surv = 'H6ID7C', 
         diag_hypert_age_w6_surv = 'H6ID7CA', 
         curr_hypert_w6_surv = 'H6ID7CC', 
         diag_depress_w6_surv = 'H6ID7G',
         diag_depress_age_w6_surv = 'H6ID7GA',
         curr_depress_w6_surv = 'H6ID7GC',
         depress_scale_w6_surv = cesd,
         diag_diab_w6_home = 'H6CQ47B', 
         diag_hypert_w6_home = 'H6CQ47A',
         hypert_joint_w6_home = 'H6HTENJC', 
         hypert_joint_w6_home_alt = 'H6HTENJC2',
         diab_joint_w6_home = 'H6DIABJC') %>%
  mutate(across(where(is.numeric), ~na_if(.x, c(-9998))),
         across(where(is.numeric), ~na_if(.x, c(-9997))),
         across(where(is.numeric), ~na_if(.x, c(-9996))))

health5 <-health5 %>% select(AID, hypert_joint_w5_home, diab_joint_w5_home)
health6 <-health6 %>% select(AID, hypert_joint_w6_home, diab_joint_w6_home)


################
# Merge together.
################
df <-pheno %>%
  left_join(., clocks, by = 'aid') %>%
  left_join(., meth, by = "Barcode") %>%
  left_join(health5, by = c("aid" = "AID")) %>%
  left_join(health6, by = c("aid" = "AID"))


################
# Residualize clocks and scale
################
clock_vars <-c("DunedinPACE", "PCGrimAge", "PCPhenoAge", "AL", "DNAmAL", "methAL", "methALT")

# residualize on Age
df <-df %>% 
  filter(!if_any(c(DunedinPACE, PCGrimAge, PCPhenoAge, AL, DNAmAL, methAL, methALT), is.na)) %>%
  mutate(across(
    all_of(clock_vars), 
    ~resid(lm(.x ~ age)),
    .names = "{.col}_r"
  ))

# scale
df <-df %>%
  mutate(across(
    .cols = c(DunedinPACE_r, PCGrimAge_r, PCPhenoAge_r, AL_r, DNAmAL_r, methAL_r, methALT_r), 
    .fns = ~as.numeric(scale(.)),
    .names = "{.col}_scaled"
  ))




hyper_table5 <-tabyl(df, hypert_joint_w5_home)
hyper_table6 <-tabyl(df, hypert_joint_w6_home)
hyper_xtable <-tabyl(df, hypert_joint_w5_home, hypert_joint_w6_home)

diab_table5 <-tabyl(df, diab_joint_w5_home)
diab_table6 <-tabyl(df, diab_joint_w6_home)
diab_xtable <-tabyl(df, diab_joint_w5_home, diab_joint_w6_home)



wb <-loadWorkbook(file = "Output/Tables/Addhealth_Tables.xlsx")
addWorksheet(wb, sheetName = "Hypertension Ns W5")
writeData(wb, sheet = "Hypertension Ns W5", x = hyper_table5)
addWorksheet(wb, sheetName = "Hypertension Ns W6")
writeData(wb, sheet = "Hypertension Ns W6", x = hyper_table6)
addWorksheet(wb, sheetName = "Hypertension Cross Tab")
writeData(wb, sheet = "Hypertension Cross Tab", x = hyper_xtable)

addWorksheet(wb, sheetName = "Diabetes Ns W5")
writeData(wb, sheet = "Diabetes Ns W5", x = diab_table5)
addWorksheet(wb, sheetName = "Diabetes Ns W6")
writeData(wb, sheet = "Diabetes Ns W6", x = diab_table6)
addWorksheet(wb, sheetName = "Diabetes Cross Tab")
writeData(wb, sheet = "Diabetes Cross Tab", x = diab_xtable)


saveWorkbook(wb, "Output/Tables/Addhealth_Tables.xlsx", overwrite = T)



write.csv(
  df %>%
    filter(!if_all(.cols = c(hypert_joint_w5_home, hypert_joint_w6_home, diab_joint_w5_home, diab_joint_w6_home), is.na))%>%
    select(aid, age, sex, AL, DNAmAL:diab_joint_w6_home) %>%
    summary(), 
  here::here('Output/Tables', 'add_health_health_characteristics_summary.csv'), 
  row.names = TRUE
)


#########################
# Modelling the ability of AL (and DNAmAL and other clocks) to predict health in W6 - with and without health in W5.
#########################
# Setup model core variables
clock_vars_r_s <-c("AL_r_scaled", "DNAmAL_r_scaled", "DunedinPACE_r_scaled", "PCGrimAge_r_scaled", "PCPhenoAge_r_scaled", "methAL_r_scaled", "methALT_r_scaled")

outcomes <-c("hypert_joint_w6_home", "diab_joint_w6_home")

cells <-c("CD4T", "CD8T", "Bcell", "NK", "Mono", "Gran")

# Create the grid
model_grid <-crossing(outcome = outcomes, clock = clock_vars_r_s)

# Build and map models.
model1 <-model_grid %>%
  pmap_dfr(\(outcome, clock) {
    fit <-glm(reformulate(c("age", "sex", clock), outcome), 
              data = df, family = binomial)
   
     events <-sum(model.response(model.frame(fit)) == 1)
    
    tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
      mutate(term = if_else(term == clock, "clock", term)) %>%
      transmute(outcome, clock, term, estimate, 
                CI_lower = conf.low, CI_upper = conf.high, p.value, 
                Events = events,
                N = nobs(fit))
  })


# Build and map models.
model2 <-model_grid %>%
  mutate(w5_health = if_else(outcome == "hypert_joint_w6_home", 
                             "hypert_joint_w5_home", "diab_joint_w5_home")) %>%
  pmap_dfr(\(outcome, clock, w5_health) {
    fit <-glm(reformulate(c("age", "sex", w5_health, clock), outcome), 
              data = df, family = binomial)
    
    events <-sum(model.response(model.frame(fit)) == 1)
    
    tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
      mutate(term = if_else(term == clock, "clock", term)) %>%
      transmute(outcome, clock, term, estimate, 
                CI_lower = conf.low, CI_upper = conf.high, p.value, 
                Events = events,
                N = nobs(fit))
  })

model3 <-model_grid %>%
  mutate(w5_health = if_else(outcome == "hypert_joint_w6_home", 
                             "hypert_joint_w5_home", "diab_joint_w5_home")) %>%
  pmap_dfr(\(outcome, clock, w5_health) {
    fit <-glm(reformulate(c("age", "sex", w5_health, cells, clock), outcome), 
              data = df, family = binomial)
    
    events <-sum(model.response(model.frame(fit)) == 1)
    
    tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
      mutate(term = if_else(term == clock, "clock", term)) %>%
      transmute(outcome, clock, term, estimate, 
                CI_lower = conf.low, CI_upper = conf.high, p.value, 
                Events = events,
                N = nobs(fit))
  })



wb <-loadWorkbook(file = "Output/Tables/Addhealth_Tables.xlsx")
addWorksheet(wb, sheetName = "health6_m1")
writeData(wb, sheet = "health6_m1", x = model1)
addWorksheet(wb, sheetName = "health6_m2_health5")
writeData(wb, sheet = "health6_m2_health5", x = model2)
addWorksheet(wb, sheetName = "health6_m3_health5_cells")
writeData(wb, sheet = "health6_m3_health5_cells", x = model3)

saveWorkbook(wb, "Output/Tables/Addhealth_Tables.xlsx", overwrite = T)


#####################
# Get AUC
#####################

get_auc <-function(outcome, clock, additional_vars = NULL, model) {
  fit <-glm(reformulate(c("age", "sex", additional_vars, clock), outcome), 
            data = df, family = binomial)
  
  observed <-model.response(model.frame(fit))
  predicted <-fitted(fit)
  
  tibble(
    model, 
    outcome, 
    clock,
    AUC = as.numeric(pROC::auc(observed, predicted)), 
    N = nobs(fit),
    Events = sum(observed == 1)
  )
}

# AUC for model 1
auc_model1 <-model_grid %>%
  pmap_dfr(function(outcome, clock) {
    get_auc(outcome, clock, model = "Model 1")
  })


# AUC for model 2
auc_model2 <-model_grid %>%
  mutate(w5_health = if_else(outcome == "hypert_joint_w6_home", 
                             "hypert_joint_w5_home", "diab_joint_w5_home")) %>%
  pmap_dfr(function(outcome, clock, w5_health) {
    get_auc(outcome, clock, additional_vars = w5_health, model = "Model 2")
  })

# AUC for model 3
auc_model3 <-model_grid %>%
  mutate(w5_health = if_else(outcome == "hypert_joint_w6_home", 
                             "hypert_joint_w5_home", "diab_joint_w5_home")) %>%
  pmap_dfr(function(outcome, clock, w5_health) {
    get_auc(outcome, clock, additional_vars = c(w5_health, cells), model = "Model 3")
  })


auc_results <-bind_rows(auc_model1, auc_model2, auc_model3)


wb <-loadWorkbook(file = "Output/Tables/Addhealth_Tables.xlsx")
addWorksheet(wb, sheetName = "auc_health_models")
writeData(wb, sheet = "auc_health_models", x = auc_results)
saveWorkbook(wb, "Output/Tables/Addhealth_Tables.xlsx", overwrite = T)
