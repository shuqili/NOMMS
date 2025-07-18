library(Maaslin2)
library(tidyverse)
library(glue)

setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/")

# Set up variables:

Anatomical_Location = "throat"
anatomicallocation = "Throat"
folder.output <- "03-throat"

sex = "Males"


read_tsv("09-MaAsLin/01-input/metadata-v5-2yr-EDSS-50subjectPool.tsv") %>%
  filter(Anatomical_Location==glue("{anatomicallocation}")) %>%
  rename(sample_ID = `#Sample ID` ) %>%
  mutate(edss_change = EDSS_yr2 - EDSS_yr0) %>%   #select(edss_change)
  filter(sex == "Males") %>%
  write_tsv("../01-NOMMS_Final_qiime_v202309/09-MaAsLin/01-input/metadata-v5-2yr-EDSS-12MalesPool-throat.tsv")





df_input_metadata = read.table(file             = glue("09-MaAsLin/01-input/metadata-v5-2yr-EDSS-12MalesPool-{Anatomical_Location}.tsv"), 
                               header           = T, 
                               sep              = "\t", 
                               row.names        = 1,
                               stringsAsFactors = FALSE)
df_input_metadata[1:5, ]

df_input_data = read.table(file             = "09-MaAsLin/01-input/input_genus_rel_abund.tsv",
                           header           = TRUE,
                           sep              = "\t", 
                           row.names        = 1,
                           stringsAsFactors = FALSE)

df_input_data[1:5, 1:5]

# ASV relative abundance as data has already been normalized.
# Fit 1: LM-AST
# Fit 2: LM-LOG
# Fit 3: LM-LOGIT
# Fit 4: CPLM-NONE

fixed_effects = "edss_change"

dir.create("09-MaAsLin/03-output/03-throat/l6taxa/male-edss_change/")
# fit_data_1 = Maaslin2(output = glue("09-MaAsLin/03-output/{folder.output}/l6taxa/female-edss_change/lmfit-ast-transformed-relab-asv-to-edss_change"),
#                       input_data     = df_input_data, 
#                       input_metadata = df_input_metadata, 
#                       
#                       analysis_method = "LM",
#                       normalization  = "NONE",
#                       transform = "AST",
#                       min_prevalence = 0.1,
#                       standardize = TRUE,
#                       cores = 4,
#                       plot_heatmap = F,
#                       plot_scatter = F,
#                       fixed_effects  = c(glue("{fixed_effects}"), "age","sex","bmi"))

fit_data_2 = Maaslin2(output = glue("09-MaAsLin/03-output/{folder.output}/l6taxa/male-edss_change/lmfit-log-transformed-relab-asv-to-edss_change"), 
                      input_data     = df_input_data, 
                      input_metadata = df_input_metadata, 
                      analysis_method = "LM",
                      normalization  = "NONE",
                      transform = "LOG",
                      min_prevalence = 0.1,
                      standardize = TRUE,
                      cores = 4,
                      plot_heatmap = F,
                      plot_scatter = F,
                      fixed_effects  = c(glue("{fixed_effects}"), "age","bmi"))

fit_data_3 = Maaslin2(output = glue("09-MaAsLin/03-output/{folder.output}/l6taxa/male-edss_change/lmfit-logit-transformed-relab-asv-to-edss_change"), 
                      input_data     = df_input_data, 
                      input_metadata = df_input_metadata, 
                      analysis_method = "LM",
                      normalization  = "NONE",
                      transform = "LOGIT",
                      min_prevalence = 0.1,
                      standardize = TRUE,
                      cores = 4,
                      plot_heatmap = F,
                      plot_scatter = F,
                      fixed_effects  = c(glue("{fixed_effects}"), "age","bmi"))





read_tsv("09-MaAsLin/01-input/input_genus_rel_abund.tsv") %>%
  filter(!str_detect(genus, "p__Bacteroidota.c__Bacteroidia.o__Flavobacteriales.f__Weeksellaceae.g__Bergeyella")) %>%
  filter(!str_detect(genus,"p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Tannerellaceae.g__Tannerella")) %>%
  filter(!str_detect(genus,"p__Proteobacteria.c__Gammaproteobacteria.o__Burkholderiales.f__Burkholderiaceae.g__Lautropia")) %>%
  filter(!str_detect(genus,"p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Pasteurellaceae.__")) %>%
  filter(!str_detect(genus,"p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Paludibacteraceae.g__F0058")) %>%
  filter(!str_detect(genus,"p__Firmicutes.c__Clostridia.o__Peptostreptococcales.Tissierellales.f__Anaerovoracaceae.g__.Eubacterium._brachy_group")) %>%
  write_tsv("../03b-EbioMedicine-Revise/00-github repository/03-differential analyses/09-males-EDSS-change/filtered.genus.rel.abund.12males.throat.tsv")
  
  
  
df_input_data_trimmed = read.table(file  = "../03b-EbioMedicine-Revise/00-github repository/03-differential analyses/09-males-EDSS-change/filtered.genus.rel.abund.12males.throat.tsv",
                           header           = TRUE,
                           sep              = "\t",
                           row.names        = 1,
                           stringsAsFactors = FALSE)




fit_data_4 = Maaslin2(output = glue("09-MaAsLin/03-output/{folder.output}/l6taxa/male-edss_change/cplmfit-no-transformation-relab-asv-to-edss_change"), 
                      input_data     = df_input_data_trimmed, 
                      input_metadata = df_input_metadata, 
                      analysis_method = "CPLM",
                      normalization  = "NONE",
                      transform = "NONE",
                      min_prevalence = 0.1,
                      standardize = TRUE,
                      cores = 1,
                      plot_heatmap = F,
                      plot_scatter = F,
                      fixed_effects  = c(glue("{fixed_effects}"), "age","bmi"))

