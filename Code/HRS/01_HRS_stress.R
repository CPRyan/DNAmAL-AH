# Packages
pacman::p_load(tidyverse, haven, sjlabelled, GGally, sjPlot, broom, gt, gtsummary, purrr, patchwork, openxlsx)



############################
###### Load Data
############################

# Outcomes
outcomes <-read_csv("Data/HRS/HRS_outcome_data.csv") %>% 
  mutate(smoker1213 = if_else(smoker1213 == 0, "Never", 
                              if_else(smoker1213 == "1", "Former", 
                                      if_else(smoker1213 == "2", "Current", NA)))) 
names(outcomes)

# sle
sle <- read_dta("Data/HRS/SLE.dta") %>% 
  as_numeric(hhid, pn) 

names(sle)

# key
key <-read_csv("Data/HRS/HRS_key_ids_barcode.csv") %>% 
  as_numeric(pn)



############################
####### Merge files
############################
hrs <-key %>% 
  left_join(., outcomes, by = "IDATid") %>% 
  left_join(., sle)

 
hrs <-hrs %>% select(hhid, pn, IDATid, Female, PAGE, smoker1213, DNAmAL, methAL, methALT, DunedinPACE, PCGrimAge, PCPhenoAge, sle_score1, sle_score2, sle_score3) %>% 
  mutate(smoker1213 = fct_relevel(smoker1213, "Never", "Former", "Current"))

# Age residualized clocks
clock_vars <- c(
  "DunedinPACE",
  "PCGrimAge",
  "PCPhenoAge",
  "DNAmAL",
  "methAL",
  "methALT"
)

hrs <- hrs %>%
  mutate(across(
    all_of(clock_vars),
    ~ resid(lm(.x ~ PAGE)),
    .names = "{.col}_r"
  ))


# Scale
hrs <- hrs %>%
  mutate(across(
    .cols = c(DNAmAL_r, methAL_r, methALT_r, DunedinPACE_r, PCGrimAge_r, PCPhenoAge_r, sle_score1, sle_score2, sle_score3
              ),   # replace with your variable list
    .fns  = ~ as.numeric(scale(.)), 
    .names = "{.col}_scaled"
  ))

############################
####### EDA
############################
table_all <-hrs %>% 
  tbl_summary(include = c(PAGE, Female, smoker1213, DNAmAL, methAL, methALT, DunedinPACE, PCGrimAge, PCPhenoAge, sle_score2 
                          #, sle_score2, sle_score3
                          ), 
              statistic = all_continuous() ~ "{mean} ({sd}) [{min}-{max}]",
              label = list(
                PAGE ~ "Age (Years)",
                Female ~ "Sex (F)", 
                smoker1213 ~ "Smoking Status")) %>% 
  add_n()



table1 <-hrs %>% 
  filter(!is.na(sle_score1)) %>% 
  tbl_summary(include = c(PAGE, Female, smoker1213, DNAmAL, methAL, methALT, DunedinPACE, PCGrimAge, PCPhenoAge, sle_score1), 
              statistic = all_continuous() ~ "{mean} ({sd}) [{min}-{max}]",
              label = list(
                DNAmAL ~ "DNAmAL-AH",
                PAGE ~ "Age (Years)",
                Female ~ "Sex (F)", 
                smoker1213 ~ "Smoking Status",
                sle_score1 ~ "SLE Score")) %>% 
  add_n()

table2 <-hrs %>% 
  filter(!is.na(sle_score2)) %>% 
  tbl_summary(include = c(PAGE, Female, smoker1213, DNAmAL, methAL, methALT, DunedinPACE, PCGrimAge, PCPhenoAge, sle_score2), 
              statistic = all_continuous() ~ "{mean} ({sd}) [{min}-{max}]",
              label = list(
                DNAmAL ~ "DNAmAL-AH",
                PAGE ~ "Age (Years)",
                Female ~ "Sex (F)", 
                smoker1213 ~ "Smoking Status",
                sle_score2 ~ "SLE Score")) %>% 
  add_n()

table3 <-hrs %>% 
  filter(!is.na(sle_score3)) %>% 
  tbl_summary(include = c(PAGE, Female, smoker1213, DNAmAL, methAL, methALT, DunedinPACE, PCGrimAge, PCPhenoAge, sle_score3), 
              statistic = all_continuous() ~ "{mean} ({sd}) [{min}-{max}]",
              label = list(
                DNAmAL ~ "DNAmAL-AH",
                PAGE ~ "Age (Years)",
                Female ~ "Sex (F)", 
                smoker1213 ~ "Smoking Status",
                sle_score3 ~ "SLE Score")) %>% 
  add_n()

table4 <-hrs %>% 
  filter(!is.na(sle_score2)) %>% 
  filter(!is.na(smoker1213)) %>% 
  tbl_summary(include = c(PAGE, Female, smoker1213, DNAmAL, methAL, methALT, DunedinPACE, PCGrimAge, PCPhenoAge, sle_score2), 
              statistic = all_continuous() ~ "{mean} ({sd}) [{min}-{max}]",
              label = list(
                DNAmAL ~ "DNAmAL-AH",
                PAGE ~ "Age (Years)",
                Female ~ "Sex (F)", 
                smoker1213 ~ "Smoking Status",
                sle_score2 ~ "SLE Score")) %>% 
  add_n()



tables <- list(
  "SLE score 1" = table1,
  "SLE score 2" = table2,
  "SLE score 3" = table3,
  "SLE score 2 - no missing" = table4
)

output_dir <- "Output/HRS/Tables"
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
  file.path(output_dir, "HRS_SLE_summary_tables.xlsx"),
  overwrite = TRUE
)




############################
# Histograms
############################

histogramz <-hrs %>% 
  select(-c(contains('scaled'))) %>% 
  pivot_longer(c(DNAmAL, methAL, methALT, sle_score2 #1, sle_score2, sle_score3
                 ), values_to = 'sle_score_value', names_to = 'sle_score') %>% 
  ggplot(., aes(x = sle_score_value, fill = sle_score))+geom_histogram()+facet_wrap(~sle_score, scales = 'free_x')+theme_bw()+theme(legend.position = "none") +
  paletteer::scale_fill_paletteer_d("MetBrewer::Tiepolo") 

ggsave(histogramz, filename = "Output/HRS/Figures/HRS_SLE_histograms.pdf")


# Correlation matrix of all AL scores and SLE scores
ggp <-ggpairs(hrs[, c("DNAmAL_r_scaled", "methAL_r_scaled", "methALT_r_scaled", "DunedinPACE_r_scaled", "PCGrimAge_r_scaled", "PCPhenoAge_r_scaled", "sle_score2_scaled" #,"sle_score2_scaled", "sle_score3_scaled"
                      )])+theme_bw()

# Plots to match Add health
ggsave(ggp, filename = "Output/HRS/Figures/HRS_correlation_matrix.pdf")

##############
## Boxpplots of Adult SLE
##############
source("Code/HRS/f_plot_sle_boxplot.R")

sle_labels <- c(
  sle_score1 = "SLE score 1",
  sle_score2 = "SLE score 2",
  sle_score3 = "SLE score 3"
)

dnam_labels <- c(
  DNAmAL  = "DNAmAL-AH",
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

output_dir <- "Output/HRS/Figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sle_plots <- list()

for (sle_var in names(sle_labels)) {
  for (y_var in names(dnam_labels)) {
    plot_name <- paste0(y_var, "_by_", sle_var)
    
    sle_plots[[plot_name]] <- plot_sle_boxplot(
      hrs, sle_var, y_var, sle_labels[[sle_var]], dnam_labels[[y_var]]
    )
    
    ggsave(
      file.path(output_dir, paste0(plot_name, "_boxplot.pdf")),
      sle_plots[[plot_name]], width = 6, height = 5
    )
  }
}



hrs %>%
  filter(!is.na(sle_score2) & !is.na(DNAmAL)) %>% 
  mutate(sle2_cluster = if_else(sle_score2 >= 8, 8, sle_score2)) %>% 
  ggplot(., aes(x = factor(sle_score2), y = DNAmAL,
             fill = sle_score2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, shape = 21, colour = "black", alpha = 0.2) +
  scale_fill_gradientn(colours = grDevices::hcl.colors(8, "Viridis", rev = TRUE)) +
  labs(x = "Number of stressful life events", y = "DNAmAL-AH") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("Output/HRS/Figures/DNAmAL_by_sle_score2_boxplot.pdf", width =5, height =4.5)


# Ready for models. 
source("Code/HRS/02_HRS_models.R")

