####################
# Set up library and cores:
####################
.libPaths("/ifs/sec/cpc/addhealth/users/cpr/R/x86_64-pc-linux-gnu-library/4.3")
# Where to find my custom package library


# Setup paths
cpr_home <-"/ifs/sec/cpc/addhealth/users/cpr/AL_clock"
epig_home <-"/ifs/sec/cpc/addhealth/epig"

####################
# Load packages
####################
pacman::p_load(tidyverse, data.table, sjlabelled, SuperLearner, caret, glmnet, RhpcBLASctl, randomForest, biglasso, ranger)

####################
# Setup Data
####################

###########
# Phenotype Data
###########
pheno <-read_csv(paste0(cpr_home, "/Output/Data/al_scores_cells_manifest_for_ewas.csv"))

###########
# EWASa: Version 1 has platform compatible AND Sugden Probes
###########
# Add Sugden probes
sugden_epic_450k <-read_csv(here::here("Data", "probelist_450k_epicv1_epicv2_sugden_ying.csv"))


epic450k_probes <-sugden_epic_450k %>%
  filter(four_50_k == 1 & epicv1 == 1 & epicv2 == 1 & sugden == 1) %>%
  pull(probe)

###########
# EWASb: Version 2 has MSA
###########

# Add MSA probes
msa_probes <-read_csv(here::here("Output/Data", "platform_msa_compatible.csv")) %>% pull(probes)
# 96129

msa_compatible <-epic450k_probes[epic450k_probes %in% msa_probes]

###########
# EWAS and reliable450k/EPIC/EPICv2/Sugden probes
###########
al <-read_csv(here::here("Output/Data", "AL_mean_ewas_toptable_full.csv"))
sig_al <-al %>% filter(adj.P.Val < 0.05) 
al_legacy_msa <-read_csv(here::here("Output/Data", "AL_mean_ewas_toptable_legacy_sugden_msa.csv"))
sig_al_legacy_msa <-al_legacy_msa %>% filter(adj.P.Val < 0.05) 





# Only use MSA-compatible legacy probes
key_probes <- sig_al_legacy_msa %>% pull(cg)

probes_list <-list(al_msa_probes = sig_al_legacy_msa %>% pull(cg))

save(probes_list, file = here::here("Output/Data", "ewas_sugden_probe_list.Rdata"))

print('saved probe list')

###########
# Beta matrix
###########
chunk1  <-fread(paste0(epig_home, "/betas/", "postqc.beta.batch1.csv")) %>% column_to_rownames("cg")
chunk2a <-fread(paste0(epig_home, "/betas/", "postqc.beta.batch2a.csv")) %>% column_to_rownames("cg")
chunk2b <-fread(paste0(epig_home, "/betas/", "postqc.beta.batch2b.csv")) %>% column_to_rownames("cg")
chunk3  <-fread(paste0(epig_home, "/betas/", "postqc.beta.batch3.csv")) %>% column_to_rownames("cg")
chunk4  <-fread(paste0(epig_home, "/betas/", "postqc.beta.batch4.csv")) %>% column_to_rownames("cg")

print("all betas are loaded separately")

# Bind and select only those in pheno
betas <-bind_cols(chunk1, chunk2a, chunk2b, chunk3, chunk4) %>% dplyr::select(any_of(c(pheno$Barcode)))

print('betas are merged into one')

# remove chunks and empty trash
rm(chunk1, chunk2a, chunk2b, chunk3, chunk4);gc()

# Save betas in EWAS and Sugden
write_csv(as_tibble(betas, rownames = 'probe'), here::here("Output/Data", "betas_from_ewas_and_sugden_probes.csv"))

print('mission accomplished')

