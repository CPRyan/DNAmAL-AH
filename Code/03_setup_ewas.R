
####################
# Load packages
####################
pacman::p_load(tidyverse, DataExplorer, sjlabelled, naniar, GGally)

####################
# Load Data
####################
epig_home <-"/ifs/sec/cpc/addhealth/epig"

##########
# Load scaled data
dat_scaled <-read_csv(here::here("Output/Data", "AL_dat_scaled_narm.csv"))  %>% sjlabelled::as_character(aid)

##########
# Load cell counts
cells <-read_csv(paste0(epig_home, "/clocks/clocks.nobatchcorrection.postdetp.fordissemination.092525.csv"))%>% 
  dplyr::select(aid, cell_CD4T, cell_CD8T, cell_Bcell, cell_NK, cell_Mono, cell_Gran) %>% 
  dplyr::rename(CD4T = cell_CD4T, 
                CD8T = cell_CD8T, 
                NK =  cell_NK, 
                Bcell = cell_Bcell, 
                Mono = cell_Mono, 
                Gran = cell_Gran) %>%
  sjlabelled::as_character(aid)


##########
# Load manifest
manifest <-read_csv(paste0(epig_home, "/manifest/", 'manifest.allBatches.csv'))
# All samples with data = 4798
# This includes 53 technical replicates, 204 biological replicates, and 11 QC failures. 

manifest %>% 
  group_by(FinalQCStatus) %>% count()

#   FinalQCStatus            n
#    <chr>                <int>
# 1 cross_contamination      6 (did not used to be present)
# fail:gender_mismatch       1
# 2 flag:bc                  6
# 3 flag:mu                  4
# 4 pass                  4728 (used to be 4734)
# 5 pass:replicate          53

manifest_clean <-manifest %>% 
  filter(FinalQCStatus == "pass") %>% # keep only samples that passed QC - no replicates
  mutate(Barcode = paste0(ChipBarcode, "_", ChipPos)) %>%
  as_character(aid)

####################
# Filter to data in both
####################
df_full <-reduce(list(dat_scaled, cells, manifest_clean), inner_join)


write_csv(df_full, here::here("Output/Data", "al_scores_cells_manifest_for_ewas.csv"))




