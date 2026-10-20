# Load packages
pacman::p_load(tidyverse, DataExplorer, sjlabelled, naniar, GGally)


# Load functions:
# Scale (mean=0, sd =1) and Sum
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

# Setup directories
home<- "/ifs/sec/cpc/addhealth/users/cpr/AL_clock"
in.folder<- "/ifs/sec/cpc/addhealth/epig/allostatic_load/01_Prep_Datasets/input"
out.folder<- "/ifs/sec/cpc/addhealth/users/cpr/AL_clock/Output"
setwd(home)

# Load raw AL data
dat <- read_csv(file= file.path(in.folder,"w5bio.radhika.062821"))

# Clean up names (remove h5)
names(dat) <- sub("^h5", "", names(dat))

summary(dat)
# many with values of 9999, clearly a code for missing 

# Convert 9990 and over to NA
dat <- dat %>% 
  as_character(aid)%>%
  mutate(across(where(is.numeric), ~ if_else(. > 9990, NA_real_, .)))
summary(dat)

write_csv(dat, here::here("Output/Data", "unscaled_partly_cleaned_biomarkers.csv"))




# Transformations:
dat <-dat %>%
  mutate(hba1c_log = log(hba1c), 
         tg_log = log(tg), 
         crp_log = log(crp), 
         il6_log = log(il6), 
         tnf_log = log(tnf), 
         hdl_inv = hdl * -1)%>%
  select(-c(hba1c, tg, crp, il6, tnf, hdl))





# Check at missingness overlaps
gg_miss_upset(dat, nsets = 20, set_size.show = TRUE)
# All are lab measures. 
# Must save plot manually for some reason

gg_miss_fct(x = dat %>%
              select(-c(aid, age, sex)), fct = race) + labs(title = "NA in AL measures relative to race")
ggsave(here::here("Output/Figures", "NA_in_AL_measures_by_race.pdf"))

gg_miss_fct(x = dat %>%
              select(-c(aid, age)), fct = sex) + labs(title = "NA in AL measures relative to race")
ggsave(here::here("Output/Figures", "NA_in_AL_measures_by_sex.pdf"))

ggplot(dat, aes(x = bmi, y = crp_log))+
  geom_miss_point()+
  facet_wrap(~sex)+theme_bw()+scale_color_manual(values = c("darkgoldenrod3", 'mediumpurple4'))
ggsave(here::here("Output/Figures", "NA_in_crp_by_bmi_sex.pdf"))

dat %>% na.omit() %>% nrow()
# 4245 with complete data


# Plot with Age, stratified by sex
dat %>%
  pivot_longer(-c(aid, age, sex, race), values_to = "Values", names_to = 'AL Component') %>%
  ggplot(., aes(x = age, y = Values, color = sex))+
  geom_point(aes(fill = `AL Component`), pch = 21, color = 'black', alpha = 0.4)+
  geom_smooth(method = 'lm', alpha = 0.5)+
  facet_wrap(~`AL Component`, scales = 'free_y')+
  paletteer::scale_fill_paletteer_d("ggthemes::calc")+
  scale_color_manual(values = c("red", 'gray'))+
  theme_bw()+
  xlab("Age") + ylab("Allostatic Load Value")+
  theme(legend.position = 'none') 

ggsave(here::here("Output/Figures", "AL_measures_by_Age.pdf"), width = 10, height = 8)


dat_scaled <-dat %>%
  mutate_at(c('bmi', 'crp_log', 'dbp', 'hba1c_log', 'hdl_inv', 'il6_log', 'pr', 'sbp', 'tc', 'tg_log', 'tnf_log', 'waist'), scale2, na.rm = TRUE)
  

# Any > 5 SDs from mean? 
dat_scaled %>%
  select(-c(age))%>%
  summarise(across(where(is.numeric), ~ sum(abs(.) > 5, na.rm = TRUE))) 
#  sbp   dbp    pr    tc   bmi  waist hba1c_log tg_log  crp_log il6_log tnf_log hdl_inv
# <int> <int> <int> <int> <int> <int>     <int>   <int>   <int>   <int>   <int>   <int>
#   4     0     2     7     1     0        58      0       0      10      16       3

outliers <-dat_scaled %>%
  filter(if_any(-c(aid, age, sex, race), ~ . > 5 | . < -5))




# Make a column for "outliers" and save the dat results (in case anything looks weird later)
dat_w_outlier <-dat %>%
  mutate(outlier_status = if_else(aid %in% outliers$aid, "outlier", 'not_outlier'))
        

create_report(dat_w_outlier %>%
                select(-c(aid)), output_dir = here::here("Output/Out"), output_file = "al_key_vars_report_outlier.html")

dat_w_outlier %>%
  write_csv(here::here("Output/Data", "AL_outlierflag.csv"))



