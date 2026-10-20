pacman::p_load(tidyverse, DataExplorer, sjlabelled, UpSetR, GGally, patchwork, openxlsx, gtsummary)


df_full <-read_csv(here::here("Output/Data", "al_scores_cells_manifest_for_ewas.csv")) %>% select(aid, age, sex, AL_mean)

dat <-read_csv(here::here("Output/Data", "unscaled_partly_cleaned_biomarkers.csv")) %>% select(-age, -sex)

tab_dat <-left_join(df_full, dat, by = 'aid')

# Having issue recreating this table since they updated R in Longleaf and fs package will not install and so DataExplorer will not install. Sigh. 
tab_dat %>%
  select(-c(aid))%>%
 tbl_summary() %>%
  as_gt()%>%
  gt::gtsave(filename = here::here("Output/Tables", "AddHealth_DNAm_Sample_Characteristics.html"))





