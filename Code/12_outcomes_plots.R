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
pacman::p_load(tidyverse, data.table, sjlabelled, caret, glmnet, openxlsx, GGally, patchwork, paletteer)


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
    .cols = c(AL_mean, AL_mean_al_msa_probes_pred, methALT),   # replace with your variable list
    .fns  = ~ as.numeric(scale(.)), 
    .names = "{.col}_scaled"
  ))




# 1) Select variables of interest in your specified order
corr_df <- df %>%
  dplyr::select(age, 
         AL_mean,
         AL_mean_al_msa_probes_pred,
         sbp:tc,
         bmi:hdl_inv)

names(corr_df) <-c("Age", "Allostatic Load", "DNAmAL", "SBP", "DBP", "Pulse Rate", "Tcells", "BMI", "Waist Circ.", "HbA1c", 'Triglycerides', "CRP", "IL6", "TNF-a", "HDL")



plot_corr_heatmap <- function(df, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  
  vars <- c("Age", "Allostatic Load", "DNAmAL", "SBP", "DBP", "Pulse Rate", "Tcells", "BMI", "Waist Circ.", "HbA1c", 'Triglycerides', "CRP", "IL6", "TNF-a", "HDL")
  
  # 2) Correlation matrix
  R <- cor(corr_df, use = "pairwise.complete.obs", method = method)
  
  # 3) Long format + helper columns
  df_long <- as.data.frame(R) %>%
    rownames_to_column("var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "r") %>%
    mutate(
      var1 = factor(var1, levels = vars),
      var2 = factor(var2, levels = vars),
      r2   = r^2,
      i    = as.integer(var1),
      j    = as.integer(var2),
      upper = i < j,          # strictly upper triangle
      diag  = i == j
    )
  
  # 4) Heatmap of r, with r² labels on the upper triangle
  ggplot(df_long, aes(x = var1, y = var2, fill = r)) +
    geom_tile(color = "white", linewidth = 0.3) +
    # r² labels only on the upper triangle (omit diagonal)
    geom_text(
      data = dplyr::filter(df_long, upper),
      aes(label = sprintf("%.2f", r2)),
      size = 3
    ) +
    scale_fill_gradient2(
      name = "r",
      limits = c(-1, 1),
      midpoint = 0,
      low = "dodgerblue3", mid = "white", high = "firebrick3"
    ) +
    scale_y_discrete(limits = rev(levels(df_long$var2))) +
    coord_fixed() +
    labs(x = NULL, y = NULL,
         title = "Correlation heatmap (tiles: r; upper triangle labels: r²)") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right"
    )
}

# Use it:
plot_corr_heatmap(df, method = "pearson")

ggsave("Output/Figures/Correlation_Matrix_AL_measures_individual_components.pdf", width = 7)



corr_df %>%
  ggpairs(., 
          lower = list(continuous = function(data, mapping, ...) {
            ggplot(data = data, mapping = mapping) +
              geom_point(alpha = 0.2, size = 1, color = 'deepskyblue4') +  # Explicit transparency
              geom_smooth(method = "lm", se = FALSE, color = "black")
          }),
          upper = list(continuous = wrap("cor", size = 5)),
          diag = list(continuous = wrap("barDiag", alpha = 0.5))) +
  theme_minimal()


ggsave("Output/Figures/Correlation_Matrix_AL_measures_individual_components.pdf", width = 12, height = 12)




###################
# Basic descriptive plots
###################
## Check summary
df %>% dplyr::select(contains('sle'), -contains('cat'))%>% summary()


# Set colors
my_colors <- c(
  "sle_earlylife" = "#1f78b4",  # blue
  "sle_adult"     = "#ff7f00",  # orange
  "sle"           = "#33a02c"   # teal
)

# Distribution plots of Stress variables
p1 <-ggplot(df, aes(x = sle_earlylife)) + 
  geom_histogram(binwidth = 1, color = "#1f78b4", fill = "#1f78b4", alpha =0.2) +
  theme_minimal() +
  labs(x = "SLE Count Early Life", y = "Frequency")

p2 <-ggplot(df, aes(x = sle_adult)) +
  geom_histogram(binwidth = 1, color = "#ff7f00", fill = "#ff7f00", alpha = 0.2) +
  theme_minimal() +
  labs(x = "SLE Count Adulthood", y = "Frequency")

p3 <-ggplot(df, aes(x = sle)) + 
  geom_histogram(binwidth = 1, color = "#33a02c", fill = "#33a02c", alpha=0.2) +
  theme_minimal() +
  labs(x = "SLE Count (All)", y = "Frequency")

p1|p2|p3

ggsave(here::here("Output/Figures", "Distribution_of_SLEs.pdf"), width = 10, height = 6)

###########
# GGpairs
###########

p <-ggpairs(df[, c('sex', "sle", "sle_adult", "sle_earlylife", "AL_mean_scaled", "AL_mean_al_msa_probes_pred_scaled", "methALT_scaled")],
mapping = aes(color = sex, fill = sex),
lower = list(continuous = function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.2, size = 1) +
    geom_smooth(method = "lm", se = FALSE)
  }),
upper = list(continuous = wrap("cor", size = 3, method = "pearson")),
diag = list(continuous = wrap("barDiag", alpha = 0.5))
)

# Add manual color + fill scales globally
p + 
  scale_color_manual(values = c('firebrick', 'steelblue3')) +
  scale_fill_manual(values = c('firebrick', 'steelblue3')) +
  theme_minimal()

ggsave(here::here("Output/Figures", "Correlation_of_AL_measures_and_SLE.pdf"), width = 12, height = 12)


##############
## Boxpplots of Adult SLE
##############

a1 <-df %>% 
  filter(!is.na(sle_adult)) %>%
  ggplot(., aes(x = as_factor(sle_adult), y = AL_mean_scaled, fill = sle_adult)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, color = 'black', pch = 21, alpha = 0.2, aes(fill = sle_adult))+
  theme_minimal() +
  labs(x = "Adult SLEs", y = "Allostatic Load (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) 

a2 <-df %>% 
  filter(!is.na(sle_adult)) %>%
  ggplot(., aes(x = as_factor(sle_adult), y = AL_mean_al_msa_probes_pred_scaled, fill = sle_adult)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, color = 'black', pch = 21, alpha = 0.2, aes(fill = sle_adult))+
  theme_minimal() +
  labs(x = "Adult SLEs", y = "AddHealth DNAmAL (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) 

  
a3 <-df %>% 
  filter(!is.na(sle_adult)) %>%
  ggplot(., aes(x = as_factor(sle_adult), y = methALT_scaled, fill = sle_adult)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, color = 'black', pch = 21, alpha = 0.2, aes(fill = sle_adult))+
  theme_minimal() +
  labs(x = "Adult SLEs", y = "MethALT (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) 

a1|a2|a3

ggsave(here::here("Output/Figures", "AL_measures_by_Adult_SLE_boxplots.pdf"), width = 12, height = 6)

#################
# Boxplots of early life SLE
#################

b1 <-df %>% 
  filter(!is.na(sle_earlylife)) %>%
  ggplot(., aes(x = as_factor(sle_earlylife), y = AL_mean_scaled, fill = sle_earlylife)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, color = 'black', pch = 21, alpha = 0.2, aes(fill = sle_earlylife))+
  theme_minimal() +
  labs(x = "Early Life SLE", y = "Allostatic Load (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) 

b2 <-df %>% 
  filter(!is.na(sle_earlylife)) %>%
  ggplot(., aes(x = as_factor(sle_earlylife), y = AL_mean_al_msa_probes_pred_scaled, fill = sle_earlylife)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, color = 'black', pch = 21, alpha = 0.2, aes(fill = sle_earlylife))+
  theme_minimal() +
  labs(x = "Early Life SLE", y = "AddHealth DNAmAL (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) 


b3 <-df %>% 
  filter(!is.na(sle_earlylife)) %>%
  ggplot(., aes(x = as_factor(sle_earlylife), y = methALT_scaled, fill = sle_earlylife)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, color = 'black', pch = 21, alpha = 0.2, aes(fill = sle_earlylife))+
  theme_minimal() +
  labs(x = "Early Life SLE", y = "MethALT (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) 

b1|b2|b3

ggsave(here::here("Output/Figures", "AL_measures_by_EarlyLife_SLE_boxplots.pdf"), width = 12, height = 6)

###############
# All SLE boxplots
###############

c1 <-df %>% 
  filter(!is.na(sle)) %>%
  ggplot(., aes(x = as_factor(sle), y = AL_mean_scaled, fill = sle)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, color = 'black', pch = 21, alpha = 0.2, aes(fill = sle))+
  theme_minimal() +
  labs(x = "Number of Stressful Life Events", y = "Allostatic Load (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) 

c2 <-df %>% 
  filter(!is.na(sle)) %>%
  ggplot(., aes(x = as_factor(sle), y = AL_mean_al_msa_probes_pred_scaled, fill = sle)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, color = 'black', pch = 21, alpha = 0.2, aes(fill = sle))+
  theme_minimal() +
  labs(x = "Number of Stressful Life Events", y = "AddHealth DNAmAL (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) 


c3 <-df %>% 
  filter(!is.na(sle)) %>%
  ggplot(., aes(x = as_factor(sle), y = methALT_scaled, fill = sle)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, color = 'black', pch = 21, alpha = 0.2, aes(fill = sle))+
  theme_minimal() +
  labs(x = "Combined SLE", y = "MethALT (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) 

c1|c2

ggsave(here::here("Output/Figures", "AL_measures_by_ALL_SLE_boxplots.pdf"), width = 9, height = 4)


##############
# Adult SLE continuous
##############
d1 <-df %>% 
  filter(!is.na(sle_adult)) %>%
  ggplot(., aes(x = sle_adult, y = AL_mean_scaled)) +
  geom_jitter(width = 0.2, color = 'black', pch = 21, alpha = 0.2, aes(fill = AL_mean_scaled))+
  geom_smooth(method = 'lm', color = 'darkorange4', linetype = 'dashed')+
  theme_minimal() +
  labs(x = "Adult SLE", y = "Allostatic Load (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) +
  scale_x_continuous(breaks = 0:12)

d2 <-df %>% 
  filter(!is.na(sle_adult)) %>%
  ggplot(., aes(x = sle_adult, y = AL_mean_al_msa_probes_pred_scaled)) +
  geom_jitter(width = 0.2, color = 'black', pch = 21, alpha = 0.2, aes(fill = AL_mean_al_msa_probes_pred_scaled))+
  geom_smooth(method = 'lm', color = 'darkorange4', linetype = 'dashed')+
  theme_minimal() +
  labs(x = "Adult SLE", y = "AddHealth DNAmAL (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) +
  scale_x_continuous(breaks = 0:12)


d3 <-df %>% 
  filter(!is.na(sle_adult)) %>%
  ggplot(., aes(x = sle_adult, y = methALT_scaled)) +
  geom_jitter(width = 0.2, color = 'black', pch = 21, alpha = 0.2, aes(fill = methALT_scaled))+
  geom_smooth(method = 'lm', color = 'darkorange4', linetype = 'dashed')+
  theme_minimal() +
  labs(x = "Adult SLE", y = "MethALT (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) +
  scale_x_continuous(breaks = 0:12)

d1|d2|d3

ggsave(here::here("Output/Figures", "AL_measures_by_Adult_SLE_continuous.pdf"), width = 12, height = 6)


##############
# Early life SLE continuous
##############
e1 <-df %>% 
  filter(!is.na(sle_earlylife)) %>%
  ggplot(., aes(x = sle_earlylife, y = AL_mean_scaled)) +
  geom_jitter(width = 0.2, color = 'black', pch = 21, alpha = 0.2, aes(fill = AL_mean_scaled))+
  geom_smooth(method = 'lm', color = 'darkorange4', linetype = 'dashed')+
  theme_minimal() +
  labs(x = "Early life SLE", y = "Allostatic Load (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) +
  scale_x_continuous(breaks = 0:9)

e2 <-df %>% 
  filter(!is.na(sle_earlylife)) %>%
  ggplot(., aes(x = sle_earlylife, y = AL_mean_al_msa_probes_pred_scaled)) +
  geom_jitter(width = 0.2, color = 'black', pch = 21, alpha = 0.2, aes(fill = AL_mean_al_msa_probes_pred_scaled))+
  geom_smooth(method = 'lm', color = 'darkorange4', linetype = 'dashed')+
  theme_minimal() +
  labs(x = "Early life SLE", y = "AddHealth DNAmAL (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) +
  scale_x_continuous(breaks = 0:9)


e3 <-df %>% 
  filter(!is.na(sle_earlylife)) %>%
  ggplot(., aes(x = sle_earlylife, y = methALT_scaled)) +
  geom_jitter(width = 0.2, color = 'black', pch = 21, alpha = 0.2, aes(fill = methALT_scaled))+
  geom_smooth(method = 'lm', color = 'darkorange4', linetype = 'dashed')+
  theme_minimal() +
  labs(x = "Early life SLE", y = "MethALT (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) +
  scale_x_continuous(breaks = 0:9)

e1|e2|e3

ggsave(here::here("Output/Figures", "AL_measures_by_EarlyLife_SLE_continuous.pdf"), width = 12, height = 6)


##############
# ALL SLE continuous
##############
f1 <-df %>% 
  filter(!is.na(sle)) %>%
  ggplot(., aes(x = sle, y = AL_mean_scaled)) +
  geom_jitter(width = 0.1, color = 'black', pch = 21, alpha = 0.2, aes(fill = AL_mean_scaled))+
  geom_smooth(method = 'lm', color = 'darkorange4', linetype = 'dashed')+
  theme_minimal() +
  labs(x = "Combined SLE", y = "Allostatic Load (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) +
  scale_x_continuous(breaks = 1:18)

f2 <-df %>% 
  filter(!is.na(sle)) %>%
  ggplot(., aes(x = sle, y = AL_mean_al_msa_probes_pred_scaled)) +
  geom_jitter(width = 0.1, color = 'black', pch = 21, alpha = 0.2, aes(fill = AL_mean_al_msa_probes_pred_scaled))+
  geom_smooth(method = 'lm', color = 'darkorange4', linetype = 'dashed')+
  theme_minimal() +
  labs(x = "Combined SLE", y = "AddHealth DNAmAL (scaled)")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) +
  scale_x_continuous(breaks = 1:18)


f3 <-df %>% 
  filter(!is.na(sle)) %>%
  ggplot(., aes(x = sle, y = methALT_scaled)) +
  geom_jitter(width = 0.1, color = 'black', pch = 21, alpha = 0.2, aes(fill = methALT_scaled))+
  geom_smooth(method = 'lm', color = 'darkorange4', linetype = 'dashed')+
  theme_minimal() +
  labs(x = "Combined SLE", y = "MethALT_scaled")+
  theme(legend.position = 'none')+
  paletteer::scale_fill_paletteer_c("grDevices::Sunset", direction = -1) +
  scale_x_continuous(breaks = 1:18)

f1|f2|f3

ggsave(here::here("Output/Figures", "AL_measures_by_combined_SLE_continuous.pdf"), width = 12, height = 6)
