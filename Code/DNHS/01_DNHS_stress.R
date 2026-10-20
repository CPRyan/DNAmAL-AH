####################
# Load packages
####################
pacman::p_load(tidyverse, haven, sjlabelled, GGally, sjPlot, broom, gt, gtsummary, purrr, patchwork, openxlsx, geepack)

####################
# Load Data
####################

# Load DNHS demographic information
demo <-read_sas("Data/DNHS/dnhs_demo.sas7bdat")
stress <-read_sas("Data/DNHS/dnhs_stress.sas7bdat")
smoke <-read_sas("Data/DNHS/dnhs_smoke.sas7bdat")

# Load DHNS DNAm data and measures
dnhs <-read_csv("Data/DNHS/DNAmAL_measures/DNHS_full_pheno_dataset.csv")
dnamAL <-read_csv("Data/DNHS/DNAmAL_measures/DNHS_DNAmAL.csv") %>% rename(DNAmAL = lambda.min)
methAL <-read_csv("Data/DNHS/DNAmAL_measures/DNHS_MethAL_ALT.csv")
pace <-read_csv("Data/DNHS/DNHS_CPR_noob_DunedinPACE.csv")
pc_clocks <-read_csv("Data/DNHS/DNHS_CPR_noob_PC_clocks.csv")



dnhs_dnam <-left_join(dnhs, dnamAL) %>% 
  left_join(., methAL) %>% 
  left_join(., pace) %>% 
  left_join(., pc_clocks) %>% 
  # Keep only those samples that passed CU QC. 
  filter(passed_cu_core_QC == "pass") 

my_df <-left_join(dnhs_dnam, demo) %>% 
  left_join(., smoke) %>% 
  left_join(., stress)

####################
# Clean data
####################


# Source code for matching age with wave, taking summed wave data, or if there are multiple matching waves to age, taking the mean of all waves. 
source("Code/DNHS/01b_clean_waves.R")


# Overview of how many people have more than one wave with data
# One person has no waves, most have 1, and a few have 2 waves that match DNAm age.
df_focal %>%
  count(n_age_matches)
#   n_age_matches     n
#             <dbl> <int>
# 1             0     1
# 2             1   795
# 3             2    51


# Subset to key columns

df_focal <-df_focal %>% 
  select(RESP, Female, Age, n_age_matches, focal_wave, focal_age, focal_smoking, DNAmAL, methAL, methALT, DunedinPACE, PCGrimAge, PCPhenoAge, focal_trauma, focal_stress) %>% 
  # One sample is missing age - 847 becomes 846
  filter(!is.na(Age))


# Age residualized clocks
clock_vars <- c(
  "DunedinPACE",
  "PCGrimAge",
  "PCPhenoAge",
  "DNAmAL",
  "methAL",
  "methALT"
)

df_focal <- df_focal %>%
  mutate(across(
    all_of(clock_vars),
    ~ resid(lm(.x ~ Age)),
    .names = "{.col}_r"
  ))


# Scale
df_focal <- df_focal %>%
  mutate(across(
    .cols = c(DNAmAL_r, methAL_r, methALT_r, DunedinPACE_r, PCGrimAge_r, PCPhenoAge_r, focal_stress, focal_trauma),   # replace with your variable list
    .fns  = ~ as.numeric(scale(.)), 
    .names = "{.col}_scaled"
  ))


df_focal_nodups <-df_focal %>% 
  # Remove Duplicates
  group_by(RESP) %>%
  slice_max(Age, n = 1, with_ties = FALSE) %>%
  ungroup() 

############################
####### EDA
############################
table_all <-df_focal %>% 
  tbl_summary(include = c(Age, Female, focal_smoking, DNAmAL, methAL, methALT, DunedinPACE, PCGrimAge, PCPhenoAge, focal_trauma, focal_stress), 
              statistic = all_continuous() ~ "{mean} ({sd})",
              label = list(
                Age ~ "Age (Years)",
                Female ~ "Sex (F)", 
                focal_smoking ~ "Smoking Status", 
                focal_trauma ~ "Trauma Score", 
                focal_stress ~ "SLE Score")) %>% 
  add_n()


table_nodups <-df_focal_nodups %>% 
  tbl_summary(include = c(Age, Female, focal_smoking, DNAmAL, methAL, methALT, DunedinPACE, PCGrimAge, PCPhenoAge, focal_trauma, focal_stress), 
              statistic = all_continuous() ~ "{mean} ({sd})",
              label = list(
                Age ~ "Age (Years)",
                Female ~ "Sex (F)", 
                focal_smoking ~ "Smoking Status", 
                focal_trauma ~ "Trauma Score", 
                focal_stress ~ "SLE Score")) %>% 
  add_n()



tables <- list(
  "All Observations" = table_all,
  "Last Observation (no dups)" = table_nodups
)

output_dir <- "Output/DNHS/Tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save each table as HTML
iwalk(tables, \(tbl, name) {
  filename <- paste0(str_replace_all(name, " ", "_"), "_summary.html")
  gtsave(as_gt(tbl), file.path(output_dir, filename))
})

# Save all tables as separate Excel sheets
wb <- createWorkbook()
header_style <- createStyle(textDecoration = "bold", fgFill = "#D9E2F3",
                            halign = "center", border = "bottom")

iwalk(tables, \(tbl, name) {
  x <- as_tibble(tbl, col_labels = TRUE)
  
  # Remove Markdown formatting from the Excel column names
  names(x) <- names(x) %>%
    str_remove_all("\\*\\*") %>%
    str_replace_all("\\s*\\n\\s*", " ")
  
  addWorksheet(wb, name, gridLines = FALSE)
  writeData(wb, name, x, headerStyle = header_style)
  setColWidths(wb, name, cols = 1:ncol(x), widths = "auto")
  freezePane(wb, name, firstActiveRow = 2)
})

saveWorkbook(
  wb,
  file.path(output_dir, "DNHS_summary_tables.xlsx"),
  overwrite = TRUE
)


############################
# Histograms
############################


histogramz <- df_focal_nodups %>% 
  select(-c(contains('scaled'))) %>% 
  pivot_longer(c(DNAmAL, methAL, methALT, focal_trauma, focal_stress), values_to = 'stress_trauma_score_value', names_to = 'stress_trauma_score') %>% 
  mutate(stress_trauma_score = fct_relevel(stress_trauma_score, "DNAmAL", "methAL", "methALT")) %>% 
  ggplot(., aes(x = stress_trauma_score_value, fill = stress_trauma_score))+geom_histogram()+facet_wrap(~stress_trauma_score, scales = 'free_x')+theme_bw()+theme(legend.position = "none") +
  paletteer::scale_fill_paletteer_d("MetBrewer::Tiepolo") 

ggsave(histogramz, filename = "Output/DNHS/Figures/DNHS_SLE_histograms_nodups.pdf")


# Correlation matrix of all AL scores and SLE scores
ggp <-ggpairs(df_focal_nodups[, c("DNAmAL_r_scaled", "methAL_r_scaled", "methALT_r_scaled", "DunedinPACE_r_scaled", "PCGrimAge_r_scaled", "PCPhenoAge_r_scaled", "focal_stress", "focal_trauma")])+theme_bw()

# Plots to match Add health
ggsave(ggp, filename = "Output/DNHS/Figures/DNHS_correlation_matrix_nodups.pdf")




##############
## Boxpplots of Adult SLE
##############
source("Code/HRS/f_plot_sle_boxplot.R")

sle_labels <- c(
  focal_trauma = "Trauma Score",
  focal_stress = "Stress Score"
)

dnam_labels <- c(
  DNAmAL  = "DNAmAL",
  methAL  = "methAL",
  methALT = "methALT"
)


plot_sle_boxplot <- function(data, sle_var, y_var, sle_label, y_label) {
  data %>%
    filter(!is.na(.data[[sle_var]]), !is.na(.data[[y_var]])) %>%
    ggplot(aes(x = factor(.data[[sle_var]]), y = .data[[y_var]],
               fill = .data[[sle_var]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, shape = 21, colour = "black", alpha = 0.2) +
    scale_fill_gradientn(colours = grDevices::hcl.colors(10, "Viridis", rev = TRUE)) +
    labs(x = sle_label, y = y_label) +
    theme_minimal() +
    theme(legend.position = "none")
}

output_dir <- "Output/DNHS/Figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sle_plots <- list()

for (sle_var in names(sle_labels)) {
  for (y_var in names(dnam_labels)) {
    plot_name <- paste0(y_var, "_by_", sle_var)
    
    sle_plots[[plot_name]] <- plot_sle_boxplot(
      df_focal_nodups, sle_var, y_var, sle_labels[[sle_var]], dnam_labels[[y_var]]
    )
    
    ggsave(
      file.path(output_dir, paste0(plot_name, "_boxplot.pdf")),
      sle_plots[[plot_name]], width = 6, height = 5
    )
  }
}













