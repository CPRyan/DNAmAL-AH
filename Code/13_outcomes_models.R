####################
# Set up library and cores:
####################
.libPaths("/ifs/sec/cpc/addhealth/users/cpr/R/x86_64-pc-linux-gnu-library/4.3")
# Where to find my custom package library

# Setup parallel computation - use all cores on our computer. (removed because it was resulting in errors)
# (num_cores = RhpcBLASctl::get_num_cores())
# options(mc.cores = num_cores)
options(expressions = 5e5)

# Setup paths
cpr_home <-"/ifs/sec/cpc/addhealth/users/cpr/AL_clock"
epig_home <-"/ifs/sec/cpc/addhealth/epig"

####################
# Load packages
####################
pacman::p_load(tidyverse, data.table, sjlabelled, caret, glmnet, openxlsx, GGally, patchwork, paletteer, purrr, broom, sjPlot)


###################
# Load data and merge
####################

# Load allostatic load predictors
al_preds <-read_csv(here::here("Output/Data", "AL_phenos_with_predictions.csv"))%>% sjlabelled::as_character(aid)

# Load Clock Data (already has)
clocks <-read_csv(paste0(epig_home, "/clocks/clocks.nobatchcorrection.postdetp.fordissemination.092525.csv"))%>% dplyr::select(aid, DunedinPACE, PCGrimAge, PCPhenoAge) %>% sjlabelled::as_character(aid)

# Load Chamberlain
chamberlain <-read_csv(here::here("Output/Data", "AddHealth_MethAL_ALT.csv"))

# Load stress measures:
stress <-read_csv("/ifs/sec/cpc/addhealth/users/share/aid_sle.csv")%>% sjlabelled::as_character(aid)

#######
# Merge
#######
df <-left_join(al_preds, clocks, by = c('aid')) %>%
  left_join(., chamberlain, by = 'Barcode') %>%
  left_join(., stress, by = 'aid') 

#######
# Scale mean 0 SD 1
#######
df <- df %>%
  mutate(across(
    .cols = c(AL_mean, AL_mean_al_msa_probes_pred,  methALT),   # replace with your variable list
    .fns  = ~ as.numeric(scale(.)), 
    .names = "{.col}_scaled"
  ))


#######
# Models
#######

# the X candidates you want to swap in:
x_vars <- names(df)[52:89] # <- replace with your X names

# names of the two outcomes
y1 <- "AL_mean_scaled"
y2 <- "AL_mean_al_msa_probes_pred_scaled"

# shared covariates (keep train_test*sle as an interaction term)
covars <- c('age', "sex","CD4T","CD8T","NK","Bcell","Mono","Gran","train_test*sle")

build_formula <- function(y, x, covars) {
  reformulate(termlabels = c(x, covars), response = y)
}

fit_and_tidy <- function(data, y, x, covars, model_label) {
  fm <- build_formula(y, x, covars)
  m  <- lm(fm, data = data)
  tidy(m, conf.int = TRUE) %>%
    filter(term == x) %>%                   # assumes X is continuous (see note below)
    mutate(model = model_label, outcome = y, X = x, .before = 1)
}

# SAFE wrapper: if lm/tidy fails, return a single NA row for this (model, X)
safe_fit_and_tidy <- function(data, y, x, covars, model_label) {
  tryCatch(
    fit_and_tidy(data, y, x, covars, model_label),
    error = function(e) {
      tibble(
        model = model_label, outcome = y, X = x, term = x,
        estimate = NA_real_, std.error = NA_real_, statistic = NA_real_,
        p.value = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
        .error = conditionMessage(e)
      )
    }
  )
}

# run all models for each X across both outcomes
results <- map_dfr(
  x_vars,
  ~ bind_rows(
    safe_fit_and_tidy(df, y1, .x, covars, "m1"),
    safe_fit_and_tidy(df, y2, .x, covars, "m2")
  )
) %>% mutate(
  model = dplyr::recode(model, m1 = "AL", m2 = "DNAmAL"),
  model = factor(model, levels = c("AL","DNAmAL"))
)




wb <-loadWorkbook(file = "Output/Tables/Addhealth_Tables.xlsx")
addWorksheet(wb, sheetName = "SLE_outcomes_models")
writeData(wb, sheet = "SLE_outcomes_models", x = results)
saveWorkbook(wb, "Output/Tables/Addhealth_Tables.xlsx", overwrite = T)

results %>% filter(p.value < 0.05) %>% View()



results %>% filter(X %in% c("sle_earlylife", "sle_adult", "sle"))


results %>% filter(X == 'sle') %>%
ggplot(., aes(x = model, y = estimate, fill = model)) +
  geom_col() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  facet_wrap(~ X, scales = "free_y") +
  scale_fill_manual(values = c(AL = "darkgoldenrod4", DNAmAL = "darkorchid4"))+
  labs(x = NULL, y = "Estimated effect of stressor on AL measure", title = "Effect sizes with 95% CI (upper whisker)") +
  theme_bw(base_size = 12)+theme(legend.position = 'none')

ggsave("Output/Figures/SLE_effect_estimates_AL_DNAmAL.pdf", width = 6, height = 6)

######################
# Actual models for SLE
######################

df <-df %>%
  mutate(smoker = fct_relevel(smoker, c('never', 'former', 'current')))

fm <- build_formula("AL_mean_scaled", "sle", covars)
fm2 <- build_formula("AL_mean_al_msa_probes_pred_scaled", "sle", covars)

m1  <- lm(fm, data = df)
m2  <- lm(fm2, data = df)

tab_model(m1, m2, pred.labels = c("Intercept", "SLE", "Age", "Sex [M]", "CD4T", "CD8T", "NK", "Bcell", "Mono", "Gran",  "Test/Train [Training]", "SLE*Test/Train"), dv.labels = c("Allostatic Load", "DNAmAL"), file = here::here("Output/Tables/SLE_models_comparison.html"))

