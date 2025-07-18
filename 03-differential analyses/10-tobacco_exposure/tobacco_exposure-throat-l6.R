library(tidyverse)
library(Maaslin2)

getwd()

setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/")

# 
read_tsv("metadata-CL04-NOMMS-v5.txt") %>%
  filter(Anatomical_Location=="T") %>%
  rename(sample_ID = sample_name) %>%
  select(sample_ID, SubjectID, smoking, everything()) %>%
  mutate(tobacco_exposure = case_when(str_detect(smoking, "Ciga") ~ "tobacco exposed",
                                      str_detect(smoking, "Toba") ~ "tobacco exposed")) %>%
  select(sample_ID, SubjectID, smoking, tobacco_exposure, everything()) %>%
  mutate(tobacco_exposure = if_else(is.na(tobacco_exposure), "not exposed", tobacco_exposure)) %>% 
  # view
  write_tsv("09-MaAsLin/01-input/metadata-v5-throat-tobacco_exposure.tsv")

df_input_metadata = read.table(file             = "09-MaAsLin/01-input/metadata-v5-throat-tobacco_exposure.tsv", 
                               header           = TRUE, 
                               sep              = "\t", 
                               row.names        = 1,
                               stringsAsFactors = FALSE)%>%
  mutate(smoking_YN = if_else(smoking_YN == 0, "no", "yes"))
df_input_metadata[1:5, ]

df_input_data = read.table(file             = "09-MaAsLin/01-input/input_genus_rel_abund.tsv",
                           header           = TRUE,
                           sep              = "\t", 
                           row.names        = 1,
                           stringsAsFactors = FALSE)

df_input_data[1:5, 1:5]

# ASV relative abundance as data has already been normalized.

# Fit 4: CPLM-NONE

dir.create("09-MaAsLin/03-output/03-throat/l6taxa/tobacco_exposure + [age + sex + bmi ]")
fit_data_4 = Maaslin2::Maaslin2(output = "09-MaAsLin/03-output/03-throat/l6taxa/tobacco_exposure + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-diseasestat/", 
                                input_data     = df_input_data, 
                                input_metadata = df_input_metadata, 
                                analysis_method = "CPLM",
                                normalization  = "NONE",
                                transform = "NONE",
                                min_prevalence = 0.1,
                                standardize = TRUE,
                                cores = 4,
                                plot_heatmap = F,
                                plot_scatter = F,
                                fixed_effects  = c("age","sex","bmi",
                                                   "tobacco_exposure"),
                                reference = c("tobacco_exposure","not exposed"))









































