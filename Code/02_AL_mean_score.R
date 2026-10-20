
####################
# Load packages
####################
pacman::p_load(tidyverse, sjlabelled, psych, car, gtsummary, openxlsx)

# Scale (mean=0, sd =1) and Sum
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

####################
# Load raw AL data
####################
dat_outlier <-read_csv(here::here("Output/Data", "AL_outlierflag.csv"))


# Remove outliers and scale - reason is using scaled values before losing nearly 900 samples seems like a bad ideas
dat_al <-dat_outlier %>%
  # remove missing
  na.omit() %>%
  # Scale all variables (mean 0, sd 1)
  mutate_at(c('bmi', 'crp_log', 'dbp', 'hba1c_log', 'hdl_inv', 'il6_log', 'pr', 'sbp', 'tc', 'tg_log', 'tnf_log', 'waist'), scale2, na.rm = TRUE) %>%
  # Take the mean of all variables to make AL score
  mutate(AL_mean = rowMeans(select(., sbp:tc, bmi:hdl_inv), na.rm = TRUE))

dat_al_mean <-dat_al %>%
  select(aid, AL_mean)


####################
# Write score
####################
write_csv(dat_al_mean, here::here("Output/Data", "AL_mean_score.csv"))

write_csv(dat_al, here::here("Output/Data", "AL_dat_scaled_narm.csv"))



####################

corr_matrix <-dat_al %>%
  select(-c(aid, sex, race, outlier_status)) |>
  cor() |>
  round(3)
  
wb <-createWorkbook()
addWorksheet(wb, sheetName = "biomarker_corrmat")
writeData(wb, sheet = "biomarker_corrmat", x = corr_matrix)
saveWorkbook(wb, "Output/Tables/Addhealth_Tables.xlsx", overwrite = T)



