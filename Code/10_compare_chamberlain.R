####################
# Load packages
####################
pacman::p_load(tidyverse, data.table, sjlabelled, caret, glmnet, purrr, readxl, writexl, patchwork, rlang)

source(here::here("Code/methAL_calculator.R"))

# Setup directories
# Setup paths
cpr_home <-"/ifs/sec/cpc/addhealth/users/cpr/AL_clock"
epig_home <-"/ifs/sec/cpc/addhealth/epig"

####################
# Load Data
####################

# Read the EPIC version table
meth_al_table <- read_excel("Data/Supplementary Tables - CpGs.xlsx",
sheet = "Table 2.", skip = 1)

# Read the 450k compatible version
meth_alt_table <- read_excel("Data/Supplementary Tables - CpGs.xlsx",
                             sheet = "Table 4.", skip = 1)

# Note, that in the paper it says there are 32 (EPIC) and 23 (450k compatible) CpGs in the clocks. 
# For some reason, both have 1 fewer than expected coefficients (31 and 22). 

key_probes <- setdiff(
  unique(c(meth_al_table$`CpG name`, meth_alt_table$`CpG name`)),
  "(Intercept)"
)


########
# Load Betas
########
chunk1  <-fread(paste0(epig_home, "/betas/", "beta.batch1.csv")) %>% filter(cg %in% key_probes) %>% column_to_rownames("cg")
chunk2a <-fread(paste0(epig_home, "/betas/", "beta.batch2a.csv")) %>% filter(cg %in% key_probes) %>% column_to_rownames("cg")
chunk2b <-fread(paste0(epig_home, "/betas/", "beta.batch2b.csv")) %>% filter(cg %in% key_probes) %>% column_to_rownames("cg")
chunk3  <-fread(paste0(epig_home, "/betas/", "beta.batch3.csv")) %>% filter(cg %in% key_probes) %>% column_to_rownames("cg")
chunk4  <-fread(paste0(epig_home, "/betas/", "beta.batch4.csv")) %>% filter(V1 %in% key_probes) %>% column_to_rownames("V1")

print("loaded beta value chunks")

# Bind and select only those in pheno
betas <-bind_cols(chunk1, chunk2a, chunk2b, chunk3, chunk4) 

print('merged beta value chunks')

# remove chunks and empty trash
rm(chunk1, chunk2a, chunk2b, chunk3, chunk4);gc()


########
# Predict MethAL
########

chamberlain_values <-predict_cpg_scores(betas = betas, 
                                        model_tables = list(meth_al_table, meth_alt_table), 
                                        model_names = c("methAL", "methALT"), 
                                        id_colname = "Barcode")

write_csv(chamberlain_values, here::here("Output/Data", "AddHealth_MethAL_ALT.csv"))
