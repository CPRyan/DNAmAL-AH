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
pacman::p_load(tidyverse, data.table, sjlabelled, caret, glmnet, purrr, writexl, patchwork, rlang, openxlsx, GGally)

####################
# Setup Data
####################
al_preds <-read_csv(here::here("Output/Data", "AL_phenos_with_predictions.csv"))



# Single probe set evaluation
x_vars <- c("AL_mean")

y_vars <- c("AL_mean_al_msa_probes_pred")

plot_list <- vector("list", length(x_vars))

for (i in seq_along(x_vars)) {
  x_sym <- sym(x_vars[i])
  y_sym <- sym(y_vars[i])
  
  p <- ggplot(al_preds, aes(x = !!x_sym, y = !!y_sym, group = train_test)) +
    geom_point(alpha = 0.2, pch = 21, color = 'black', aes(fill = train_test)) +
    geom_smooth(method = "lm", se = FALSE, aes(color = train_test)) +
    labs(
      x = x_vars[i],
      y = y_vars[i],
      title = paste(x_vars[i], "vs", y_vars[i])
    ) +
    scale_fill_manual(values = c('navyblue', 'firebrick4'))+
    scale_color_manual(values = c('navyblue', 'firebrick2'))+
    theme_minimal()+
    theme(legend.position='bottom')
  
  plot_list[[i]] <- p
  
  ggsave(
    filename = paste0(here::here("Output/Figures/"), "scatter_", x_vars[i], "_vs_", y_vars[i], ".pdf"),
    plot = p,
    width = 6, height = 4, dpi = 300
  )
}

# Combine and save final plot
combined_plot <- wrap_plots(plotlist = plot_list, ncol = 3)
ggsave(here::here("Output/Figures", "combined_scatter_plots.pdf"), plot = combined_plot, width = 12, height = 9, dpi = 300)





################################
# R2 for all models
################################

# Initialize list to collect results
r2_results <- vector("list", length(x_vars))

for (i in seq_along(x_vars)) {
  x_col <- x_vars[i]
  y_col <- y_vars[i]
  
  # Filter data for train and test sets
  train_data <- al_preds %>% filter(train_test == "train")
  test_data  <- al_preds %>% filter(train_test == "test")
  
  # Compute R² for train
  lm_train <- lm(as.formula(paste(y_col, "~", x_col)), data = train_data)
  r2_train <- summary(lm_train)$r.squared
  
  # Compute R² for test
  lm_test <- lm(as.formula(paste(y_col, "~", x_col)), data = test_data)
  r2_test <- summary(lm_test)$r.squared
  
  # Compute R² for test
  lm_all <- lm(as.formula(paste(y_col, "~", x_col)), data = al_preds)
  r2_all <- summary(lm_all)$r.squared
  
  # Save results
  r2_results[[i]] <- tibble(
    x_var = x_col,
    y_var = y_col,
    r2_train = r2_train,
    r2_test = r2_test, 
    r2_all = r2_all, 
  )
}

# Combine into one dataframe
r2_table <- bind_rows(r2_results)

wb <-loadWorkbook(file = "Output/Tables/Addhealth_Tables.xlsx")
addWorksheet(wb, sheetName = "train_test_r2")
writeData(wb, sheet = "train_test_r2", x = r2_table)
saveWorkbook(wb, "Output/Tables/Addhealth_Tables.xlsx", overwrite = T)

al_preds <- al_preds %>%
  mutate(train_test = factor(train_test, levels = c("train", "test")))  # ensures order


ggplot(al_preds, aes(x = AL_mean, y = AL_mean_al_msa_probes_pred, group = train_test)) +
  geom_point(alpha = 0.2, pch = 21, color = 'black', aes(fill = train_test)) +
  geom_smooth(method = "lm", se = FALSE, aes(color = train_test)) +
  labs(
    x = "Allostatic Load",
    y = "DNAmAL",
    color = "Data Source", fill = "Data Source"
  ) +
  scale_fill_manual(
    values = c('deepskyblue4', 'darkorange4'),
    labels = c("Training (R² = 0.762)", "Testing (R² = 0.544)")
  ) +
  scale_color_manual(
    values = c('deepskyblue4', 'darkorange4'),
    labels = c("Training (R² = 0.762)", "Testing (R² = 0.544)")
  ) +
  theme_bw() +
  theme(legend.position = c(0.80, 0.10))

ggsave(here::here("Output/Figures", "Allostatic_Load_v_DNAmAL_polished.pdf"), height = 5, width=5)








# Scatterplot of all (scaled) variables and the AL measure
al_preds %>%
  select(AL_mean_al_msa_probes_pred, AL_mean, age, sbp:hdl_inv, -race) %>%
  rename(Age = age, SBP = sbp, DBP = dbp, PR = pr, tcells = tc, BMI = bmi, 'Waist Circ.' = waist, 'HbA1c (log)' = hba1c_log, 'Trigs (log)' = tg_log, 'CRP (log)' = crp_log, 'IL6 (log)' = il6_log, 'TNFa (log)' = tnf_log, 'HDL (inverted)' = hdl_inv, 'Allostatic Load' = AL_mean, "DNAmAL" = AL_mean_al_msa_probes_pred) %>%
  ggpairs(., 
          lower = list(continuous = function(data, mapping, ...) {
            ggplot(data = data, mapping = mapping) +
              stat_binhex(bins = 28, na.rm = T, show.legend = F) +  
              geom_smooth(method = "lm", se = FALSE, color = "black")
          }),
          upper = list(continuous = wrap("cor", size = 5)),
          diag = list(continuous = wrap("barDiag", alpha = 0.5))) +
  theme_minimal()



 ggsave(here::here("Output/Figures", "AL_measures_scatterplot_matrix.pdf"), height = 14, width = 14)

