library(Maaslin2)
library(tidyverse)

setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/")
list.files()

read_tsv("metadata-CL04-NOMMS-v5.txt") %>%
  rename(sample_ID = sample_name) 
# write_tsv("09-MaAsLin/01-input/metadata-v5-gums.tsv")

df_input_metadata = read.table(file             = "metadata-CL04-NOMMS-v5.txt", 
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
getwd()
dir.create(path = "09-MaAsLin/03-output/00-anatomical_location/l6-ref-oropharyngeal/")
fit_data_1 = Maaslin2(output = "09-MaAsLin/03-output/00-anatomical_location/l6-ref-oropharyngeal/lmfit-ast-transformed-relab-asv-to-diseasestat",
                      input_data     = df_input_data, 
                      input_metadata = df_input_metadata, 
                      
                      analysis_method = "LM",
                      normalization  = "NONE",
                      transform = "AST",
                      min_prevalence = 0.1,
                      standardize = TRUE,
                      cores = 4,
                      fixed_effects  = c("Anatomical_Location",
                                         "age","sex","bmi",
                                         "smoking",
                                         "dental_procedure_simplified",
                                         "oral_disease_simplified",
                                         "oral_hygiene_cumsum_score"),
                      reference      = c("Anatomical_Location,T;smoking,NONE;dental_procedure_simplified,none;oral_disease_simplified,none"))

# ERROR::AST transform is only valid for values between -1 and 1. Please select an appropriate normalization option or normalize your data prior to running.



fit_data_2 = Maaslin2(output = "09-MaAsLin/03-output/00-anatomical_location/l6-ref-oropharyngeal/lmfit-log-transformed-relab-asv-to-diseasestat", 
                      input_data     = df_input_data, 
                      input_metadata = df_input_metadata, 
                      analysis_method = "LM",
                      normalization  = "NONE",
                      transform = "LOG",
                      min_prevalence = 0.1,
                      standardize = TRUE,
                      cores = 4,
                      fixed_effects  = c("Anatomical_Location",
                                         "age","sex","bmi",
                                         "smoking",
                                         "dental_procedure_simplified",
                                         "oral_disease_simplified",
                                         "oral_hygiene_cumsum_score"),
                      reference      = c("Anatomical_Location,T;smoking,NONE;dental_procedure_simplified,none;oral_disease_simplified,none"))

fit_data_3 = Maaslin2(output = "09-MaAsLin/03-output/00-anatomical_location/l6-ref-oropharyngeal/lmfit-logit-transformed-relab-asv-to-diseasestat", 
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
                      fixed_effects  = c("Anatomical_Location",
                                         "age","sex","bmi",
                                         "smoking",
                                         "dental_procedure_simplified",
                                         "oral_disease_simplified",
                                         "oral_hygiene_cumsum_score"),
                      reference      = c("Anatomical_Location,T;smoking,NONE;dental_procedure_simplified,none;oral_disease_simplified,none"))

fit_data_4 = Maaslin2(output = "09-MaAsLin/03-output/00-anatomical_location/l6-ref-oropharyngeal/cplmfit-no-transformation-relab-asv-to-diseasestat", 
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
                      fixed_effects  = c("Anatomical_Location",
                                         "age","sex","bmi",
                                         "smoking",
                                         "dental_procedure_simplified",
                                         "oral_disease_simplified",
                                         "oral_hygiene_cumsum_score"),
                      reference      = c("Anatomical_Location,T;smoking,NONE;dental_procedure_simplified,none;oral_disease_simplified,none"))
