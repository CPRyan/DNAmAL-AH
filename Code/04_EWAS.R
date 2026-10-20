####################
# Setup Packages
####################
.libPaths("/ifs/sec/cpc/addhealth/users/cpr/R/x86_64-pc-linux-gnu-library/4.3")
# Where to find my custom package library

pacman::p_load(tidyverse, data.table, limma, vroom)

####################
# Setup Data
####################
cpr_home <-"/ifs/sec/cpc/addhealth/users/cpr/AL_clock"
epig_home <-"/ifs/sec/cpc/addhealth/epig"

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
length(msa_compatible)
## 38050
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

####################
# Build annotation 
####################
annot <-fread(paste0(epig_home, '/annotation/MethylationEPIC_v-1-0_B4.csv'),
              skip=7, # Skip unused Illumina rows at top
              fill = TRUE,
              select = c("Name", "Genome_Build", "CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "GencodeCompV12_Accession", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Group", "Methyl450_Loci"))

annot <-annot %>%
  filter(Name %in% rownames(betas)) %>%
  filter(CHR != "X" & CHR!="Y") %>%
  rename(cg = Name)

print("annotation is loaded")

####################
# Match pheno rows with betas columns
####################
# Ensure all rows in pheno match column names in betas
pheno_filtered <- pheno %>% 
  filter(Barcode %in% colnames(betas)) %>%
  arrange(match(Barcode, colnames(betas)))  # Ensure correct order

print("pheno filtered to betas")

# Reorder betas columns to match pheno$Barcode order
betas_reordered <- betas[, pheno_filtered$Barcode]
rm(betas)

print("betas reordered. Original Betas removed.")

# Check if the order is now correctly aligned
alignment_check <-all(colnames(betas_reordered) == pheno_filtered$Barcode)
# Must be TRUE


####################
# Build Design Matrix
####################
design_mean = model.matrix(~AL_mean + age + sex + CD4T + CD8T + NK + Bcell + Mono + Gran, data = pheno_filtered)

print("design matrices constructed")
####################
# Run EWAS
####################

###########
# AL Mean
###########
fit1 = lmFit(betas_reordered, design_mean)
ebfit1 = eBayes(fit1)

output1 = topTable(ebfit1,
                   number=Inf,
                   coef='AL_mean') %>%
  
  rownames_to_column('cg') %>%
  inner_join(annot, by='cg') %>%
  mutate(adj.P.Val = as.numeric(as.character(ifelse(adj.P.Val==0,10^-300,adj.P.Val)))) %>%
  mutate(chr = as.numeric(CHR))

print("fit1 created, output created, next to save as csv")

write_csv(output1, paste0(cpr_home, "/Output/Data/", "AL_mean_ewas_toptable_full.csv"))

######## Repeat for platform compatible
rn <-rownames(ebfit1$coefficients)
keep <-intersect(rn, msa_compatible)
ebfit2 <-ebfit1[keep, ]

output2 = topTable(ebfit2, 
                   number=Inf,
                   coef='AL_mean') %>%
  
  rownames_to_column('cg') %>%
  inner_join(annot, by='cg') %>%
  mutate(adj.P.Val = as.numeric(as.character(ifelse(adj.P.Val==0,10^-300,adj.P.Val)))) %>%
  mutate(chr = as.numeric(CHR))

print("fit1 MSA platform compatible created, output created, next to save as csv")

write_csv(output2, paste0(cpr_home, "/Output/Data/", "AL_mean_ewas_toptable_legacy_sugden_msa.csv"))



rm(fit1, ebfit1, output1, output2); gc()

####################
# Save Log utput
####################
# Save output to a log file
log_file <- paste0(cpr_home, "/Output/Out/ewas_setup_check.log")

log_conn <- file(log_file, "w")
writeLines(c(
  "EWAS Setup Check:",
  paste("Number of matching Barcodes:", nrow(pheno_filtered)),
  paste("Total columns in betas:", ncol(betas_reordered)),
  paste("Alignment correct:", alignment_check)
), con = log_conn)
close(log_conn)

print("write lines worked")
# Run garbage collection again to free up memory
gc()

# Capture and save session info safely
sink(log_file, append = TRUE)
print(sessionInfo())
sink()

# Final memory cleanup
rm(list = ls())
gc()

