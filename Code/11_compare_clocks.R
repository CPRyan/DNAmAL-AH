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
pacman::p_load(tidyverse, data.table, sjlabelled, caret, glmnet, openxlsx, GGally, ggrastr, sjPlot, ppcor)

# Function
source(paste0(cpr_home, "/../", "helper_functions_addhealth.R"))
####################
# Setup Data
####################
# Load allostatic load predictors
al_preds <-read_csv(here::here("Output/Data", "AL_phenos_with_predictions.csv"))%>% sjlabelled::as_character(aid)

# Load Clock Data (already has)
clocks <-read_csv(paste0(epig_home, "/clocks/clocks.nobatchcorrection.postdetp.fordissemination.092525.csv"))%>% dplyr::select(aid, DunedinPACE, PCGrimAge, PCPhenoAge) %>% sjlabelled::as_character(aid)

# Load Chamberlain
chamberlain <-read_csv(here::here("Output/Data", "AddHealth_MethAL_ALT.csv")) 

df <-left_join(al_preds, clocks, by = c('aid')) %>%
  left_join(., chamberlain, by = 'Barcode') 

df <-df %>% 
  filter(!if_any(c(DunedinPACE, PCGrimAge, PCPhenoAge, AL_mean, methAL, methALT), is.na)) %>%
  mutate(DunedinPACE_r = resid(lm(DunedinPACE ~ age)), 
              PCGrimAge_r = resid(lm(PCGrimAge ~ age)), 
              PCPhenoAge_r = resid(lm(PCPhenoAge ~ age)),
              AL_mean_al_msa_probes_pred_r = resid(lm(AL_mean_al_msa_probes_pred ~ age)), 
              AL_mean_r = resid(lm(AL_mean ~ age)),
              methAL_r = resid(lm(methAL ~ age)),
              methALT_r = resid(lm(methALT ~ age)))


df2 <-df %>% filter(train_test == "test")

# 804

#####################


#####################
# Examine associations between clocks
#####################
library(GGally)
library(tidyverse)

corr_matrix <-  df2 %>%
  dplyr::select(AL_mean_r,  "AL_mean_al_msa_probes_pred_r", 
         DunedinPACE_r, PCPhenoAge_r, PCGrimAge_r, methAL_r, methALT_r) |>
  na.omit()|>
  cor() |>
  round(3)


wb <-loadWorkbook(file = "Output/Tables/Addhealth_Tables.xlsx")
addWorksheet(wb, sheetName = "al_and_clocks_corrmat")
writeData(wb, sheet = "al_and_clocks_corrmat", x = corr_matrix)
saveWorkbook(wb, "Output/Tables/Addhealth_Tables.xlsx", overwrite = T)






custom_colors <- c(train = "deepskyblue4", test = "darkorange4")

p <- ggpairs(
  df %>%
    dplyr::select(train_test, AL_mean_r,  "AL_mean_al_msa_probes_pred_r", 
          # Obsolete predictions removed - only AL_mean_al_msa_probes_pred used
           DunedinPACE_r, PCPhenoAge_r, PCGrimAge_r, methAL_r, methALT_r) %>%
  rename("Allostatic Load" = AL_mean_r, 
         "DNAmAL" = AL_mean_al_msa_probes_pred_r, 
         "Train/Test" = train_test),
  mapping = aes(color = `Train/Test`, fill = `Train/Test`),
  lower = list(continuous = function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_point_rast(alpha = 0.2, size = 1) +
      geom_smooth(method = "lm", se = FALSE)
  }),
  upper = list(continuous = wrap("cor", size = 3, method = "pearson")),
  diag = list(continuous = wrap("barDiag", alpha = 0.5))
)

# Add manual color + fill scales globally
p + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal()

ggsave(here::here("Output/Figures", "Clocks_residualized_AL_measures_scatterplot_matrix.pdf"), height = 8, width = 8)


############
# NOTE: REMOVE age by taking residuals from clocks
# 
# Compare clocks in hold out data (not all data)
# 
############

lm1 <-lm(AL_mean_r ~AL_mean_al_msa_probes_pred_r + age + sex, data = df2)
lm2 <-update(lm1, . ~ . + DunedinPACE_r + PCPhenoAge_r + PCGrimAge_r)
lm3 <-update(lm1, . ~ . + methAL_r)
lm4 <-update(lm1, . ~ . + methALT_r)

summary(lm1)
summary(lm2)
summary(lm3)
summary(lm4)

tab_model(lm1, lm2, pred.labels = c("Intercept", "DNAmAL", "Age", "Sex [M]", "DunedinPACE", "PCPhenoAge", "PCGrimAge"), dv.labels = c("Basic Model", "Epigenetic Clock Model"), file = here::here("Output/Tables/clocks_models_comparison.html"), show.est = F, show.std = TRUE)

tab_model(lm1, lm3, lm4, pred.labels = c("Intercept", "DNAmAL", "Age", "Sex [M]", "MethAL", "MethALT"), dv.labels = c("Basic Model", "MethAL", "MethALT"), file = here::here("Output/Tables/chamberlain_models_comparison.html"), show.est = F, show.std = TRUE)

