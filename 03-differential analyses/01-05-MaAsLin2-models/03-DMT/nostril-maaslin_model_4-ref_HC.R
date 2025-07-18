library(Maaslin2)
library(tidyverse)

setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/")


# read_tsv("subgroup-metadata-DMT_N_over_5-CL04-NOMMS-v5-maaslin.txt") %>% 
#   filter(Anatomical_Location=="N") %>%
#   rename(sample_ID = `#Sample ID`) %>%
#   write_tsv("09-MaAsLin/01-input/subgroup-metadata-v5-DMT_N_over_5-nostril.tsv")
list.files("09-MaAsLin/01-input/")


df_input_metadata = read.table(file             = "09-MaAsLin/01-input/subgroup-metadata-v5-DMT_N_over_5-nostril.tsv", 
                               header           = TRUE, 
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
dir.create("09-MaAsLin/03-output/01-nostril/l6taxa/")
dir.create("09-MaAsLin/03-output/01-nostril/l6taxa/DMT-N_over_5 vs HC/")
fit_data_1 = Maaslin2(output = "09-MaAsLin/03-output/01-nostril/l6taxa/DMT-N_over_5 vs HC/lmfit-ast-transformed-relab-asv-to-DMT",
                      input_data     = df_input_data, 
                      input_metadata = df_input_metadata, 
                       
                      analysis_method = "LM",
                      normalization  = "NONE",
                      transform = "AST",
                      min_prevalence = 0.1,
                      standardize = TRUE,
                      cores = 4,
                      plot_heatmap = F,
                      plot_scatter = F,
                      fixed_effects  = c("DMT_v2",
                                         "age","sex","bmi",
                                         "smoking",
                                         "dental_procedure_simplified",
                                         "oral_disease_simplified",
                                         "oral_hygiene_cumsum_score"),
                      reference      = c("DMT_v2,HC;smoking,NONE;dental_procedure_simplified,none;oral_disease_simplified,none"))


fit_data_2 = Maaslin2(output = "09-MaAsLin/03-output/01-nostril/l6taxa/DMT-N_over_5 vs HC/lmfit-log-transformed-relab-asv-to-DMT", 
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
                      fixed_effects  = c("DMT_v2",
                                         "age","sex","bmi",
                                         "smoking",
                                         "dental_procedure_simplified",
                                         "oral_disease_simplified",
                                         "oral_hygiene_cumsum_score"),
                      reference      = c("DMT_v2,HC;smoking,NONE;dental_procedure_simplified,none;oral_disease_simplified,none"))

fit_data_3 = Maaslin2(output = "09-MaAsLin/03-output/01-nostril/l6taxa/DMT-N_over_5 vs HC/lmfit-logit-transformed-relab-asv-to-DMT", 
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
                      fixed_effects  = c("DMT_v2",
                                         "age","sex","bmi",
                                         "smoking",
                                         "dental_procedure_simplified",
                                         "oral_disease_simplified",
                                         "oral_hygiene_cumsum_score"),
                      reference      = c("DMT_v2,HC;smoking,NONE;dental_procedure_simplified,none;oral_disease_simplified,none"))

fit_data_4 = Maaslin2(output = "09-MaAsLin/03-output/01-nostril/l6taxa/DMT-N_over_5 vs HC/cplmfit-no-transformation-relab-asv-to-DMT", 
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
                      fixed_effects  = c("DMT_v2",
                                         "age","sex","bmi",
                                         "smoking",
                                         "dental_procedure_simplified",
                                         "oral_disease_simplified",
                                         "oral_hygiene_cumsum_score"),
                      reference      = c("DMT_v2,HC;smoking,NONE;dental_procedure_simplified,none;oral_disease_simplified,none"))
